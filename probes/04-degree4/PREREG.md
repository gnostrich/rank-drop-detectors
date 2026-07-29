# PREREG — Probe 04: degree 4, incomplete correspondence

Committed before any scan or measurement. One round, terminal. **Exactly one
variable changes from probe 03: the territory.** Energy, coder, vocabulary
{DUP, PATTERN, TEST}, normalization and pipeline are inherited unchanged.
BRIDGE and HUB remain registered kills; **there is no minting move here**.
Expected outcome, recorded now as the directive requires: **UNDERPOWERED or
NULL** — paramodular data is sparse and the curated link set is small; a thin
result reads as the sample-size fact it probably is.

## 0. Acquisition scoping facts (metadata-level, gathered before this commit)

- `smf_newforms` (LMFDB Siegel forms, degree 2) contains **zero rows of
  weight [2, 0]** — the weight relevant to abelian surfaces. The K
  (paramodular) family starts at weight [3, 0]; on the known level-277
  instance, the weight-[3,0] general-type form's λ_p admit no affine match
  to the genus-2 curve 277.a's Frobenius traces (checked at p = 3…29), so
  those rows are not surface L-functions and are **not** substituted for F7.
  F7 is therefore whatever weight-[2,0] rows exist at download time —
  scoping says none, and if none, every F7-dependent cell degrades exactly
  as specified below. That emptiness is itself an acquisition finding.
- `lfunc_lfunctions.euler_factors` (by origin `Genus2Curve/Q/N/c`) stores
  exactly 25 good/bad Euler factors (p = 2…97), arithmetically normalized
  (verified against independent point counts on 277.a at p = 3, 5, 7).

## 1. Substrate — families, bounds, pinning

Pinned as probe 03 pinned: committed snapshots + SHA256 manifest from the
beta.lmfdb.org API, unversioned-API caveat disclosed in the same words: the
API is unversioned; the committed snapshot is the substrate of record.

| block | family | source | band |
|---|---|---|---|
| F5 | genus-2 curves / Q, one row per curve | `g2c_curves` + class L `euler_factors` | cond 1–10000 |
| F6 | abelian-surface isogeny classes (= g2c classes) | class L `euler_factors` | cond 1–10000 |
| F7 | Siegel/paramodular newforms, **weight [2,0]**, any family, available levels | `smf_newforms` | all available |
| F1, F2 | degree-2 reference rows | **probe 03's committed snapshots, reused byte-identical** | as probe 03 |

Entries: integer Frobenius traces a_p = −c₁ of the stored degree-4 Euler
factor; 25 rational-prime columns; degree-4 column-norm vector = 4p, so the
alphabet bound is ⌊2√(4p)⌋ = ⌊4√p⌋ (Hasse–Weil for degree 4) and degree-4
fingerprints sit in the same column space as degree-2 ones with the same
motivic-weight-1 normalization exponent. Markers at p | cond (bad Euler
factors are skipped as comparison positions and encoded as markers, never a
reason to drop a row — probe 02's reading, kept). F5 and F6 rows of the same
class share the class L-function by construction; every headroom number this
produces is billed as duplication in K3 and is **not comparable to 1.67%**.

Seal: F5/F6 conductor 9001–10000 (top decile), downloaded, flagged, never
scanned. F7: top decile by level of whatever exists (vacuous if F7 is
empty; stated either way). Probe 02/03 seals bind unchanged.

## 2. Codebook and scans — deterministic, inherited

DUP cross-family pairs allowed: F5×F6 (bookkeeping, curated by label),
{F5, F6}×F7. Same exact-agreement rule as probe 03 (equal integer values at
every shared good column, ≥ 10 shared, exhaustive), pointer costs
log₂ N_total. PATTERN as probe 02 (Π = {p+1, 0, 1, p}, l ∈ {5,7,11,13},
membership deterministic, ΔE < 0). TEST: probe-02 Φ predicates are EC/Q
predicates with no degree-4 extension; TEST is therefore evaluable on no
degree-4 family and this is disclosed, not repaired.

**H4 decomposability sieve** (runs before any residual count): for every
degree-4 row of conductor N, enumerate all factorizations N = N₁·N₂ with
F1 classes at both conductors (the minimal EC conductor 11 forces
N₁, N₂ ∈ [11, 909] ⊆ F1's band, so the sieve's coverage of product
decompositions is complete); the row is decomposable iff its trace vector
equals a_p(E₁) + a_p(E₂) at every shared good column for some candidate
pair (including N₁ = N₂ squares). "Symmetric square" of a degree-2 object
is degree 3 and cannot appear in this column space; the product-with-square
case is the implementable reading, disclosed. Sieve hits are cross-checked
(report-only) against `is_gl2_type` / `end_alg` curation.

## 3. Cells

### H (any failure ⇒ VOID)
- **H1 recomputation.** 200 random F5 rows (`default_rng(20260728)`);
  independent point counting of #C(F_p) and #C(F_p²) (χ-sums over F_p and
  F_p², complete-the-square, F_p² built as F_p[x]/(x²−r), r the least
  non-residue) at odd good primes p ∤ disc; recompute c₁ and c₂ and compare
  with the stored Euler factors; also verify the functional-equation shape
  c₃ = p·c₁, c₄ = p² at good p on the same sample. Any mismatch ⇒ VOID.
- **H2 energy reuse.** Probe-02 coder imported unchanged (probe-03's
  column-norm substitution, norms = 4p); decoder byte-identical on the new
  tables under ∅ and under the final codebook; accounting ≤ 0.005 ⇒ else VOID.
- **H3 known-link recovery.** Curated subset = LMFDB-recorded links between
  {F5, F6} and F7. Gate: blind DUP recall ≥ 95% on that subset. **The subset
  size is reported explicitly; if it is under 50 instances, this cell is a
  smoke test and cannot support a claim either way, and RESULTS.md says so
  in those words.** If the subset is empty the recall is undefined, the cell
  is a size-0 smoke test, and VOID is not triggered by vacuity.
- **H4 sieve implemented and run before K1**, as §2. Not implementable ⇒ VOID.

### N
- **N1 shuffled null.** Column shuffle within family, marginals preserved,
  `default_rng(20260728 + family_id)` with family ids F5 = 7, F6 = 8,
  F7 = 9. Gate: zero cross-family matches.
- **N2 foil power.** 1000 foils per type, seed `default_rng(20260728)`:
  (a) one good column's value shifted +1 (wrapped in-alphabet); (b) wrong
  level: fingerprint replaced by that of the next row of different
  conductor; (c) wrong weight: values doubled, clipped values marked.
  Detection = the foil fails exact-DUP against its unperturbed source row.
  Gate: ≥ 95% pooled. Reported in the same sentence as any K1 null.
- **N3 depth curve.** Cross-family match counts at comparison depth
  d ∈ {10, 25, 50, all} primes; with 25 columns available, d = 50 and
  d = all coincide with 25 and are reported as such. Computed for every
  cross-family pair class and for H4 sieve candidates. **Reported whatever
  the verdict.**

### K
- **K1 residual census.** Cross-family DUP matches at full depth, minus
  curated LMFDB links, minus H4-decomposables. Count; if non-zero, the list,
  with no interpretation (probe 01's five breaking pairs got one sentence;
  that line holds).
- **K2 power/sample.** Rows per family, primes per row, comparable
  cross-family pairs actually compared ({F5 ∪ F6} × F7), curated-link count
  for H3. **If comparable pairs < 10,000, a K1 null is a SAMPLE-SIZE
  verdict, not evidence about the territory.** Which case obtains is stated
  mechanically.
- **K3 headroom.** Compression fraction per degree-4 family beside probe
  02's 1.67% and probe 03's figures, with the duplication artifact billed
  explicitly in the table (F5 vs F6 share class L-functions by
  construction; that compression is a copy, not structure).

## 4. Ladder — mechanical, in order

1. **VOID** — any H fails. 2. **BLIND** — N2 < 95%. 3. **UNDERPOWERED** —
H, N pass, K1 empty, K2 comparable pairs < 10,000; the verdict is about the
sample and is written that way. 4. **NULL** — H, N pass, K1 empty, K2 ≥
10,000. 5. **RESIDUALS** — K1 non-empty after the sieve; list published
with per-instance N3 depth curves and no interpretation; sealed-band
confrontations are separate registered actions, not part of this probe.

No parameter changes after this commit. No follow-up round. Sealed bands
untouched. Residuals are judged against LMFDB links AND the H4 sieve.

## 5. What no verdict here establishes

UNDERPOWERED/NULL says nothing about the paramodular conjecture or the
territory's mathematics — only about this sample, this depth, this source.
RESIDUALS, if any, are candidate coincidences until their one sealed-band
confrontation, and are published without interpretation.
