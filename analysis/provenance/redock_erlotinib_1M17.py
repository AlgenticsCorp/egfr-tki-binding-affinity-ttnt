#!/usr/bin/env python3
"""
Fix 1 (R-43), Neurosnap variant — re-dock erlotinib WT reference against the ACTIVE-conformation
crystal 1M17, replacing the inactive-conformation 4HJO reference, via Neurosnap's hosted GNINA
tool (chosen over local Docker GNINA: local Docker was unavailable — GNINA v1.3.2 doesn't exist
as a public tag, and a v1.3.1 pull filled the local disk to 195MB free, crashing Docker Desktop).

IMPORTANT METHODS DEVIATION (must be disclosed in the manuscript, not silently absorbed):
Neurosnap's hosted GNINA tool (`scripts/submit_gnina_job.py`) only accepts a receptor PDB +
ligand SDF — it does NOT expose --exhaustiveness, --seed, or --autobox_ligand controls the way
the study's local Docker GNINA calls do (see step4_docking.py: exhaustiveness=32, seeds=42/43/44,
autobox padding 6A around the reference ligand). This run is therefore:
  - a SINGLE run, not a 3-seed mean (Neurosnap's controls are opaque/undocumented; repeated
    identical submissions would not constitute genuine replicates without seed control)
  - dependent on Neurosnap's own internal pocket-detection, not the study's explicit autobox
This must be reported as a limitation/methods note if these numbers reach the manuscript, and the
value should be treated as a single-estimate reference, not directly poolable with the 3-seed local
values elsewhere in ddg_summary.csv without flagging the difference.

USAGE:
  cd study2_affinity_pfs_correlation
  /usr/local/bin/python3.11 v2_prespecified/src/redock_erlotinib_1M17_neurosnap.py            # env check
  /usr/local/bin/python3.11 v2_prespecified/src/redock_erlotinib_1M17_neurosnap.py --run       # submit + poll + recompute
"""

import sys, os, csv, json, shutil, statistics, argparse, re, time
from pathlib import Path
from datetime import datetime

STUDY2_ROOT = Path(__file__).resolve().parents[2]
V1_SRC      = STUDY2_ROOT / "v1_exploratory" / "src"
SCRIPTS_DIR = STUDY2_ROOT.parent / "scripts"
NEUROSNAP_UTILS_DIR = STUDY2_ROOT.parent / "simulation_experiment" / "src" / "utils"
sys.path.insert(0, str(V1_SRC))
sys.path.insert(0, str(SCRIPTS_DIR))
sys.path.insert(0, str(NEUROSNAP_UTILS_DIR))

from step3b_prep_structures import prepare_receptor   # PDBFixer pH7.4 + OpenMM AMBER14 (1000 steps)

DATA        = STUDY2_ROOT / "v1_exploratory" / "data"
STRUCT_DIR  = DATA / "structures" / "wildtype_1M17"
LIGAND_SDF  = DATA / "ligands" / "erlotinib_prepared.sdf"
MERGED_CSV  = DATA / "correlation" / "merged_ddg_ttnt.csv"
DDG_SUMMARY = DATA / "docking" / "ddg_summary.csv"
ARCHIVE     = STUDY2_ROOT / "v2_prespecified" / "archive"
V2_ANALYSIS = STUDY2_ROOT / "v2_prespecified" / "src" / "run_v2_analysis.py"
JOB_DIR     = DATA / "docking" / "wildtype_1M17_erlotinib_neurosnap"

NEW_REF_PDB = "1M17"

def _stamp(): return datetime.now().strftime("%Y%m%d_%H%M%S")

def _archive(p: Path):
    if p.exists():
        ARCHIVE.mkdir(parents=True, exist_ok=True)
        dst = ARCHIVE / f"{p.stem}.{_stamp()}{p.suffix}"
        shutil.copy(p, dst)
        print(f"  archived {p.name} -> {dst}")

def env_check() -> bool:
    ok = True
    print("Environment check:")
    if not os.getenv("NEUROSNAP_API_KEY"):
        # try loading from repo-root .env
        from dotenv import load_dotenv
        load_dotenv(STUDY2_ROOT.parent / ".env")
    if os.getenv("NEUROSNAP_API_KEY"):
        print("  NEUROSNAP_API_KEY: present")
    else:
        print("  NEUROSNAP_API_KEY: MISSING"); ok = False
    import urllib.request
    try:
        urllib.request.urlopen(f"https://files.rcsb.org/download/{NEW_REF_PDB}.pdb", timeout=15)
        print("  RCSB: reachable")
    except Exception as e:
        print(f"  RCSB: NOT reachable ({e})"); ok = False
    try:
        import pdbfixer, openmm, rdkit  # noqa
        print("  pdbfixer/openmm/rdkit: importable")
    except Exception as e:
        print(f"  Python deps: MISSING ({e})"); ok = False
    if not LIGAND_SDF.exists():
        print(f"  erlotinib ligand SDF: MISSING ({LIGAND_SDF})"); ok = False
    else:
        print(f"  erlotinib ligand SDF: present ({LIGAND_SDF})")
    return ok

def _parse_score(text: str):
    for pattern in [r"<minimizedAffinity>\s*([-\d.]+)", r"<CNNaffinity>\s*([-\d.]+)"]:
        m = re.search(pattern, text)
        if m:
            return float(m.group(1))
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    args = ap.parse_args()

    ok = env_check()
    if not args.run:
        print("\nDry run only. Re-run with --run to execute.")
        return
    if not ok:
        print("\nEnvironment not ready. Aborting."); sys.exit(1)

    # 1) Fetch + prepare the 1M17 WT receptor
    STRUCT_DIR.mkdir(parents=True, exist_ok=True)
    raw = STRUCT_DIR / "structure.pdb"
    if not raw.exists():
        import urllib.request
        print(f"\nDownloading {NEW_REF_PDB} from RCSB...")
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{NEW_REF_PDB}.pdb", raw)
    receptor = STRUCT_DIR / "prepared.pdb"
    if not receptor.exists():
        print("Preparing 1M17 receptor (PDBFixer pH7.4 + OpenMM AMBER14, 1000-step min)...")
        prepare_receptor("wildtype_1M17", raw, is_predicted=False)
    else:
        print("1M17 receptor already prepared, reusing.")

    # 2) Submit to Neurosnap GNINA — receptor as a raw multipart file part (FFProteinStructureList),
    #    ligand as a JSON-encoded SMILES entry (FFMolList). Verified working format via manual API
    #    probing (job/submit/GNINA returns the job_id as a bare JSON string, not {"job_id": ...} —
    #    the repo's existing submit_gnina_job.py helper assumed the wrong response shape and also
    #    sent the receptor as a JSON-wrapped string, which the server rejects with "not provided").
    #    No autobox/seed control available beyond the optional Center/Size X/Y/Z search-box fields
    #    (left at Neurosnap's auto full-complex box here — see module docstring on the deviation).
    import requests

    NEUROSNAP_API_KEY = os.environ["NEUROSNAP_API_KEY"]
    ERLOTINIB_SMILES = "C(#N)c1ccc2c(ncnc2c1)Nc1cc(OCC)c(OCC)cc1"

    JOB_DIR.mkdir(parents=True, exist_ok=True)
    print("\nSubmitting erlotinib + 1M17 WT to Neurosnap GNINA...")
    ligand_field = json.dumps([{"name": "erlotinib", "type": "smiles", "data": ERLOTINIB_SMILES}])
    with open(receptor, "rb") as fh:
        resp = requests.post(
            "https://neurosnap.ai/api/job/submit/GNINA",
            headers={"X-API-KEY": NEUROSNAP_API_KEY},
            params={"note": "R43_erlotinib_1M17_WT_redock"},
            files={"Input Receptor": (receptor.name, fh, "chemical/x-pdb")},
            data={"Input Ligand": ligand_field},
        )
    (JOB_DIR / "submit_response.json").write_text(f"HTTP {resp.status_code}\n{resp.text}")
    if resp.status_code != 200:
        print(f"[ABORT] submission failed: HTTP {resp.status_code} {resp.text}"); sys.exit(2)
    job_id = resp.json()  # bare string
    print(f"  job_id = {job_id}")

    # 3) Poll until complete
    print("Polling for completion (this can take several minutes)...")
    session = requests.Session(); session.headers.update({"X-API-KEY": NEUROSNAP_API_KEY})
    deadline = time.time() + 2 * 3600
    status = None
    while time.time() < deadline:
        r = session.get(f"https://neurosnap.ai/api/job/status/{job_id}")
        r.raise_for_status()
        status = r.json()
        print(f"    status: {status}")
        if status in ("completed", "done", "success"):
            break
        if status in ("failed", "error", "cancelled"):
            print(f"[ABORT] job ended with status: {status}"); sys.exit(2)
        time.sleep(20)
    else:
        print("[ABORT] timed out after 2h waiting for job completion"); sys.exit(2)
    print(f"  final status: {status}")

    # 4) Download output.sdf (verified filename via manual probing of /job/file/{id}/out/{name} —
    #    Neurosnap's GNINA output is a single multi-pose SDF; GNINA lists poses best-first, so the
    #    FIRST <minimizedAffinity> tag is the top pose, same convention as the local Docker parser).
    r = session.get(f"https://neurosnap.ai/api/job/file/{job_id}/out/output.sdf")
    r.raise_for_status()
    out_path = JOB_DIR / "output.sdf"
    out_path.write_bytes(r.content)
    text = r.content.decode(errors="ignore")
    all_scores = [float(m) for m in re.findall(r"<minimizedAffinity>\s*\n?([-\d.]+)", text)]
    (JOB_DIR / "all_pose_scores.json").write_text(json.dumps(all_scores, indent=2))
    print(f"  all pose scores (unsorted, as returned): {all_scores}")
    # IMPORTANT: unlike the local GNINA CLI (which sorts poses best-first), Neurosnap's output.sdf
    # pose order is NOT reliably sorted by affinity (verified empirically: scores were non-monotonic
    # across two separate submissions). Take the minimum (most negative = strongest binding) instead
    # of assuming pose order — do not silently reuse the local-CLI "first tag" convention here.
    score = min(all_scores) if all_scores else None
    print(f"  best (minimum) pose score = {score}")

    if score is None:
        print("[ABORT] could not parse an affinity score from any output file; "
              f"raw files saved under {JOB_DIR} for manual inspection. Nothing changed.")
        sys.exit(3)

    new_wt_dg = round(float(score), 3)
    print(f"\nNew erlotinib WT dG on 1M17 (Neurosnap GNINA, single run) = {new_wt_dg} kcal/mol; was -6.95 on 4HJO")
    (JOB_DIR / "parsed_result.json").write_text(json.dumps({
        "job_id": job_id, "wt_dg_kcal_mol": new_wt_dg, "n_runs": 1,
        "method": "neurosnap_gnina_hosted", "note": "no seed/exhaustiveness/autobox control available"
    }, indent=2))

    # 5) Recompute ΔΔG for erlotinib rows. Archive first.
    _archive(DDG_SUMMARY); _archive(MERGED_CSV)

    rows = list(csv.DictReader(open(DDG_SUMMARY)))
    changed = []
    for r in rows:
        if r["drug"].lower() == "erlotinib":
            old = float(r["ddg"]); r["wt_dg_mean"] = f"{new_wt_dg}"
            r["ddg"] = f"{round(float(r['dg_mean']) - new_wt_dg, 3)}"
            changed.append((r["profile"], old, float(r["ddg"])))
    with open(DDG_SUMMARY, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)

    mrows = list(csv.DictReader(open(MERGED_CSV)))
    ddg_lookup = {(r["profile"], r["drug"].lower()): r["ddg"] for r in rows}
    for r in mrows:
        k = (r["profile"], r["drug"].lower())
        if r["drug"].lower() == "erlotinib" and k in ddg_lookup:
            r["ddg"] = ddg_lookup[k]
    with open(MERGED_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=mrows[0].keys()); w.writeheader(); w.writerows(mrows)

    print("\nΔΔG updated (erlotinib rows):")
    for prof, old, new in changed:
        print(f"  {prof:22s} {old:+.3f} -> {new:+.3f}")

    _archive(STUDY2_ROOT / "v2_prespecified" / "results" / "v2_all_results.json")
    print("\nRe-running group-level analysis...")
    import subprocess
    subprocess.run([sys.executable, str(V2_ANALYSIS)], check=True)

    print("\nDONE (Neurosnap single-run variant). Report the printed new v1.0 rho/p/CI, the erlotinib "
          "ΔΔG above, AND the methods deviation (single run, no autobox/seed control) so the manuscript "
          "can flag this honestly and R-43 can be closed.")

if __name__ == "__main__":
    main()
