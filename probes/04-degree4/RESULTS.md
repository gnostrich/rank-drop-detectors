# RESULTS — Probe 04: degree 4, incomplete correspondence — **UNDERPOWERED**

Run 2026-07-29, one round, terminal. Territory-only change from probe 03.
Pinned as probe 03 pinned: committed snapshots + SHA256 manifest (2,839
requests); the API is unversioned — the committed snapshot is the substrate
of record. Expected outcome was recorded in the prereg as UNDERPOWERED or
NULL; UNDERPOWERED obtained, and per the prereg it is **a verdict about the
sample, not the territory**.

## Data

Open rows: genus-2 curves **2,789** / classes **2,538** (cond ≤ 9000; 280 +
269 sealed) · **F7 = 0**: LMFDB contains **zero weight-[2,0]
Siegel/paramodular newforms** (the K family starts at weight [3,0], whose
λ_p demonstrably do not match surface traces — checked on the known 277
instance before commit 1). The F7 emptiness is the probe's acquisition
finding. EC/Q reference: 2,463 classes from probe 03's snapshots.

## The matrix

| cell | gate | result | verdict |
|---|---|---|---|
| H1 recomputation | any mismatch ⇒ VOID | 200 rows, **4,521 Euler coefficients (c₁ and c₂) recomputed over F_p and F_p²: 0 mismatches**, 0 functional-equation violations | PASS |
| H2 energy reuse | byte-identical, ≤0.005 | all streams identical (incl. degree-4 alphabets) | PASS |
| H3 known-link recovery | recall ≥95%; size reported | **curated subset size = 0.** This cell is a smoke test and cannot support a claim either way | size-0 |
| H4 decomposability sieve | implemented before K1 | **399 decomposable classes** (E₁×E₂ products); agreement with LMFDB endomorphism curation **399/399 both directions** (0 sieve hits outside product-type curation, 0 product-type classes missed) | PASS |
| N1 shuffled null | zero matches | 0 | PASS |
| N2 foil power | ≥95% | **3,000/3,000 (power 1.00)** | PASS |
| K1 residual census | count + list | 2,789 cross-matches, all curve↔class label bookkeeping; after curation and sieve: **0 residuals** | empty |
| K2 power/sample | ≥10,000 pairs | comparable degree-4 × F7 pairs = **0** | **SAMPLE-SIZE** |

**N3 depth curve** (the deliverable's most informative object; d = 50 and
d = all coincide with 25, only 25 primes exist):

| depth d | F6×F5 matches | {F5,F6}×F7 | sieve decomposables |
|---|---|---|---|
| 10 | 254 | 0 | 1 |
| 25 (= 50 = all) | 2,789 | 0 | 399 |

Counts *rise* with depth here because the ≥10-shared-primes eligibility
floor leaves no slack at d = 10 (any row with a bad prime below 29 is
ineligible); the curve measures eligibility at shallow depth, agreement at
full depth. No cross-family match existed at any depth to test for decay.

**K3 headroom** (billing stated as required): F5 = 0.695 — **duplication
artifact**: F5 curve rows share class L-functions with F6 by construction;
this fraction is a copy, NOT comparable to probe 02's 1.67%. F6 = 0.000
(no PATTERN fired on degree-4 traces). F7 = empty family.

## Verdict (mechanical)

H pass, N pass, K1 empty, K2 comparable pairs 0 < 10,000 ⇒ **UNDERPOWERED**,
terminal. The verdict is about the sample: the form side of the paramodular
correspondence has no weight-2 data in LMFDB to compare against — nothing
here is evidence about the territory's mathematics. Reported in the same
sentence as required: the K1 null sits next to N2 power = 1.00, so it is a
true sample fact, not blindness.

## What is NOT claimed

Nothing about the paramodular conjecture, in either direction. Nothing about
the territory beyond its data infrastructure: the experiment this probe was
built for cannot currently be run at scale from LMFDB because weight-2
Siegel/paramodular eigendata does not exist there — that infrastructure gap
is the finding. The H4 sieve's 399/399 agreement validates the sieve, not
any correspondence. H3 was a size-0 smoke test and supports no claim either
way. The sealed bands were never read.
