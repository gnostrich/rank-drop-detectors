# PREREG — Probe 01: congruence census over elliptic curves

Committed before any run. Every analysis decision, threshold, and seed below is
fixed here; `run.py` implements this document mechanically and the verdict is
read off the ladder in §7 with no discretion.

## 1. Question

Which pairs of non-isogenous elliptic curves over **Q** (conductor < 10000) have
congruent mod-`l` Galois representations, as detected by
`a_p(E1) ≡ a_p(E2) (mod l)` at every available good prime, for
`l ∈ {5, 7, 11, 13}`?

This is a harness run on real mathematics: congruences between modular forms
are well studied (level-lowering/raising, Eisenstein congruences), so the
**expected outcome is VALIDATED** — the detector should re-find known
structure. It is not a discovery run.

## 2. Data — pinned

- Source: Cremona's elliptic curve database, `https://github.com/JohnCremona/ecdata`
- **Pinned commit: `25cec5ecfec8b9f016eb1631ac633194c2bed39f`**
- Files used:
  - `aplist/aplist.00000-09999` — 38,042 lines, one per **isogeny class**:
    `N class a_2 a_3 ... a_97` (25 primes, p = 2..97). Entries at primes of bad
    reduction appear as `+`, `-`, or `?`.
  - `allcurves/allcurves.00000-09999` — 64,687 lines, one per **curve**:
    `N class number [a1,a2,a3,a4,a6] rank torsion`.
- Scope: conductor 1–9999 only. No other files are read.

## 3. Parsing rules

- A token that parses as an integer is a good-reduction `a_p` value. Tokens
  `+`, `-`, `?` mark bad-reduction positions; they are **skipped, never used,
  and never a reason to discard the class**. (A previous run discarded any
  curve with a bad prime below 97 and lost 95% of the data. Matching is on
  good positions only.)
- `a_p` is an isogeny invariant, so one `aplist` row serves every curve in the
  class. Curve coefficients for point counting come from the class's curve
  number 1 in `allcurves`.
- The number of primes parsed is 25; the number used per modulus is 24 (§4).
  Both are reported. The count of usable classes (§4 eligibility) is reported
  before any hit rate.

## 4. The congruence relation — fixed definition

For modulus `l`, classes `C1, C2` are **directly congruent** iff:

- comparison columns = the 25 primes minus the column `p = l` (the mod-`l`
  representation is unreliable at `p = l`), i.e. 24 columns;
- shared good positions = columns where both classes have integer entries;
- `a_p(C1) ≡ a_p(C2) (mod l)` at **every** shared good position;
- the number of shared good positions is **≥ 10** (pairs below this are
  ineligible; the count of such skipped pairs is reported).

Congruence **groups** are the connected components of the graph on all 38,042
classes with an edge per directly-congruent eligible pair (the direct relation
is not transitive in the presence of skipped positions; components are the
preregistered grouping). The pair scan is exhaustive over all pairs — no
blocking or sampling.

## 5. Cells — all in one batch, one round

### H1 — harness (must pass or the run is VOID)

1. Join integrity: every `allcurves` row joins to exactly one `aplist` row by
   `(N, class)`, and every `aplist` row has ≥ 1 curve.
2. Grouping consistency: grouping curves by fingerprint (exact integer values,
   same relation as §4 but on all 25 columns and with no modulus), every
   isogeny class must land inside a single group. Isogenous curves have equal
   `a_p` at every good prime, so failure is a data/indexing bug and is itself
   the published verdict.
3. Independent recomputation: for a uniform random sample of **200 classes
   (numpy `default_rng(20260728)`)**, recompute `a_p` at all good `p ≤ 97` by
   direct point counting over `F_p` from the curve-1 Weierstrass model and
   compare with the parsed file values. **Any mismatch → VOID.**
4. Cross-class *exact* agreement pairs (all shared good positions equal as
   integers) are counted and reported (prior run: 0 at 25 primes). This count
   is informational, not pass/fail.

### H2 — null calibration

For each `l`: independently permute each of the 25 columns of
`(value, good/bad flag)` pairs across the 38,042 rows — marginals preserved,
row structure destroyed — with numpy `default_rng(20260728 + l)`, then rerun
the §4 grouping. Report the number of non-singleton components per `l`.
Expected ≈ 0; the `l^-22`-style estimate is replaced by this empirical rate.

### A1 — census

For each `l ∈ {5, 7, 11, 13}`: run the §4 grouping on the real data. Because
rows are isogeny classes, **every non-singleton component is a non-trivial
congruence** (isogeny is already collapsed; H1.2 guarantees no class is
split). Report per `l`: number of non-singleton components, classes involved,
component size distribution.

### A2 — conductor structure

Per surviving congruence, record conductors and ratios. Pairs are classified
by the first matching rule:

- **EIS**: both members satisfy `a_p ≡ p + 1 (mod l)` at every tested good
  prime (aplist columns and the A3 extension). This is the signature of a
  reducible (Eisenstein-type) mod-`l` representation, whose conductor is
  supported away from both levels — coprime conductors are expected here.
- **P1**: `N1 = N2`.
- **P2**: `min(N1,N2) | max(N1,N2)`.
- **P3**: `gcd(N1,N2) > 1`.
- **P4**: `gcd(N1,N2) = 1` and not EIS.

Level-lowering/raising phenomenology predicts irreducible-case congruent pairs
to sit in P1–P3 (the prime-to-`l` conductor of the shared representation
divides both levels). **Pattern-consistent = EIS ∪ P1 ∪ P2 ∪ P3.** The
observed pattern is reported without further assumption; P4 pairs are recorded
individually. EIS–EIS pair counts inside a component may be computed
combinatorially (the classification is member-level); all other pairs are
enumerated explicitly. If explicit pair enumeration would exceed 5·10^6 pairs,
pattern statistics fall back to a uniform random sample of 10^6 pairs
(`default_rng(20260728)`), and this is logged.

### A3 — depth

For every class in a surviving component, compute `a_p` by point counting for
all primes `97 < p < 1000` (143 primes), skipping `p | N`. A pair is
**persistent** if congruence holds at every tested good prime for both members
(aplist columns + extension). Otherwise its statistic is the smallest breaking
prime. Pattern statistics (§ A2) are computed over persistent pairs only;
breaking pairs are reported separately as finite-depth coincidences with their
breaking depths.

## 6. Reported numbers (fixed list)

Usable class count and curve count; primes parsed/used; per-`l`: A1 component
count and sizes, H2 count, persistent-pair count and consistency fraction,
breaking-pair count and depths; H1 pass/fail with the three sub-checks;
skipped-pair count from the ≥ 10 shared-positions rule.

## 7. Interpretation ladder — applied mechanically, in this order

1. **VOID** — any of: H1.1–H1.3 fails; or for some `l`, H2 non-singleton
   count ≥ 3; or for some `l`, H2 ≥ 1 and `4·H2_l ≥ A1_l` (null comparable to
   signal).
2. **KNOWN-ONLY** — `Σ_l A1_l = 0`: no congruence beyond isogeny. Publish as a
   clean negative.
3. *(guard)* — if `Σ_l A1_l > 0` but zero pairs are persistent, the result is
   "inconclusive at this scale with this method", terminal, recorded CLOSED.
4. **VALIDATED** — ≥ 90% of persistent pairs (pooled over `l`) are
   pattern-consistent (EIS ∪ P1 ∪ P2 ∪ P3).
5. **OPEN** — otherwise: < 90% consistent. Inconsistent persistent pairs are
   recorded individually, with no interpretation beyond the data.

One round, terminal. No parameter changes after this commit; a KNOWN-ONLY or
CLOSED verdict is a completed result, not a reason to re-run.

## 8. What a VALIDATED verdict does NOT establish

Not a conjecture engine, not new mathematics. It establishes that a rank-drop
detector, run as a systematic census over an arithmetic family, recovers known
structure — the precondition for pointing it anywhere unknown, and nothing
more. In particular: trace congruence at finitely many primes is evidence for,
not proof of, an isomorphism of semisimplified mod-`l` representations; EIS
classification is a trace-pattern signature, not a computed isogeny/torsion
certificate; and conductor-pattern consistency is descriptive, not a verified
instance of Ribet's theorem.
