# RESULTS — Probe 03: cross-family transport and closure defect — **LOSSY**

Run 2026-07-28, one round, terminal. **A validation-and-measurement probe,
not a discovery probe** (named so in advance; PREREG committed before any
scan). Substrate: LMFDB via beta API, pinned as committed snapshots + SHA256
manifest (1,961 requests; the API itself is unversioned — disclosed).
Readings R1–R4 in `run.py`; notably R1: dim>1 newforms have no integer a_p
and enter as all-marker rows (3,488 of them). Artifacts in [`output/`](output/).

## Data (before any measurement)

Open rows: EC/Q **2,463** classes · newforms **5,092** (1,604 dim-1) ·
EC/K **575 + 1,045** · HMF **575 + 1,045**. Sealed (top decile, untouched):
859 / 57 / 105 / 57 / 105. All entries within Hasse alphabets: **0 violations.**

## The matrix

| cell | gate | result | verdict |
|---|---|---|---|
| H1 join/parse/recompute | 0 mismatches | 0 Hasse violations; 200-row recomputation vs Cremona: **0** | PASS |
| H2 energy reuse | byte-identical, ≤0.005 | all block streams identical; max rel 1.1×10⁻³ | PASS |
| H3 meter validity | known-commuting ≤ floor | **2,241 instances (C1+C2), every δ = 0** | PASS |
| N1 shuffled null | 0 cross matches | 0 | PASS |
| N2 wobble floor | fixed before K5 | deterministic meter ⇒ floor = 0 per type (R3) | set |
| N3 power/foils | δ > floor in ≥95% | **2,345 / 2,345 detected (power 1.00)** | PASS |
| K1 modularity, blind | prec ≥99%, recall ≥95% | **precision 1.0000, recall 1.0000** — 2,137 blind matches = 2,137 curated links; 0 residuals, 0 missed | PASS |
| K2 bc-square | expected at floor | all 104 instances δ = 0 | at floor |
| K5 closure census | fixed grid | 2,343 instances; **79 above floor (3.4%), all C3, all δ_frac = 0**; total 543 bits; concentration(top 1%) = **0.41** | **LOSSY** |

**K3 headroom (report-only).** Curve families sit at the probe-02 ceiling:
F1 0.48%, EC/K 1.83% / 1.29% (probe 02: 1.67%). Form families **halve**:
F2 46.5%, HMF 49.5% / 48.4% — every dim-1 form row is one cross-family DUP
away from free.

**K4 double dissociation.** (a) more rows, same family: entries 2 → 2,
fraction 0.52% → 0.09% — **flat**. (b) more families, same rows: entries
2 → 3,766, bits saved 1.2k → 244.2k — **growing**. Per the pre-fixed grid:
**compounding is real and scales with families, not data.** The K3 ceiling
caveat does not bite: the added form families have ≫3% headroom.

## Verdict (mechanical)

H pass, N pass, K1 pass; K5 above floor with concentration 0.41 ∈ (0.10,
0.50] ⇒ **LOSSY** — the defect is the machine's, and the instance list says
exactly where: all 79 above-floor cycles are C3 squares with δ_frac = 0 (the
two paths never disagree where both predict); the entire defect lives in
columns only the form-side path predicts — primes where the Q-curve has bad
reduction, where the good-reduction base-change formula is applied to U_p
eigenvalues outside its domain. δ quantizes at 6/12 bits (one or two escape
flags). The meter measured the solver, not the objects, and said so.
Full list: `output/cycles.csv`. Reported in the same sentence as required:
the K5 reading is interpretable because N3's power was 1.00, not blindness.

## What is NOT claimed

No new mathematics, no conjectures, no discovery: every blind match was
already a link on the corresponding LMFDB page (rediscovery, labelled as
such). The 100%/100% K1 result validates the cross-family substrate, not any
new correspondence. The closure defect is transport lossiness, not structure.
Nothing is claimed about the sealed bands or the paramodular zone; no minting
move exists here, and the transport meter added nothing to any codebook.
