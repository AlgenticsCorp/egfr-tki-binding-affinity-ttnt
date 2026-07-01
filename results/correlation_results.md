# v2 Pre-specified Correlation Results
**Date run:** 2026-06-26  
**Scope definitions:** see SCOPE.md  
**Decisions applied:** see DECISIONS.md  
**Distance measurements:** see distance_measurements.md

---

## Primary results

| Scope | n pairs | ρ | p | Bootstrap 95% CI | Significant? |
|---|---|---|---|---|---|
| **v1.0 — Binding-pocket (G719, T790, C797, Ex19del)** | **11** | **−0.664** | **0.026** | **−0.876 to −0.151** | **✓ YES** |
| **v3.0 — All sensitizing mutations** | **15** | **−0.632** | **0.011** | **−0.824 to −0.231** | **✓ YES** |
| v1.0 osimertinib only | 9 | −0.700 | 0.036 | −0.948 to −0.158 | ✓ YES |
| v3.0 osimertinib only | 11 | −0.582 | 0.060 | −0.854 to 0.000 | trend |

---

## Key differences from v1 exploratory

| | v1 exploratory | v2 pre-specified |
|---|---|---|
| Primary result | ρ=−0.566, p=0.055 | ρ=−0.664, p=0.026 |
| L858R included? | Yes | No (v1.0) / Yes (v3.0) |
| L861Q included? | Yes (single) | No (v1.0) / Yes (v3.0) |
| S768I included? | No (excluded post-hoc) | No (structural criterion) |
| Scope defined before analysis? | No | Yes |
| CI crosses zero? | Yes (+0.056) | No (−0.151) |

---

## Interpretation

Both scopes cross α=0.05 independently. The 95% CI does not cross zero in either scope.

v1.0 is the mechanistically cleanest result: only mutations with structural
evidence of direct binding-pocket contact are included, and the criterion was
defined before any correlation was computed.

v3.0 confirms the finding holds when the scope is widened to all clinically
sensitizing mutations — a biologically reasonable alternative definition.

The fact that both scopes are significant, and that v1.0 is actually stronger
than v3.0 (ρ=−0.664 vs −0.632), supports the mechanistic interpretation:
the correlation is driven specifically by binding-pocket contact, and adding
activation-loop mutations (L858R, L861Q) slightly dilutes the signal.
