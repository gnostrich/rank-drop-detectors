# PREREG — Probe 03: cross-family transport and closure defect

Committed before any scan, transport, or measurement runs. One round, terminal.
**This is a validation-and-measurement probe, not a discovery probe.** Near-zero
net compression gain is the expected outcome (probe 02's K7 ceiling: the whole
known taxonomy bought 1.67% of L₀), and is written down now so a thin result
reads as the ceiling it is.

Inherits from probe 02 (PARTIAL): vocabulary **{DUP, PATTERN, TEST}**; the
energy of `probes/02-encoder-energy/PREREG.md` §4–§5 imported unchanged and not
re-derived. `BRIDGE` and `HUB` are registered kills: not used, not repaired,
not re-proposed under another name. **This probe has no minting move**: nothing
here creates an object absent from the table. Registered kills from
`mathlib-structure-probes` also bind (no Mathlib object, no MDL judge on a
formal library, no Gibbs move in statement space).

## 1. Substrate — families, bounds, pinning

Source: **LMFDB via the `beta.lmfdb.org` API.** LMFDB does not version its
API; the achievable pin, declared here, is: (a) the exact query URLs in the
committed downloader; (b) the committed snapshot files themselves; (c) a
SHA256 manifest of every raw response, with retrieval timestamps. The
committed snapshot is the substrate of record. If any family cannot be
acquired this way, that is the finding and the probe closes there.

Family blocks, in table order (block id in parentheses):

| block | family | source | band |
|---|---|---|---|
| F1 (1) | EC/Q, one row per isogeny class (curve 1) | `ec_curvedata` (labels, ainvs); a_p **computed by point counting** | conductor 1–1000 |
| F2 (2) | classical newforms, weight 2, trivial character, **all dimensions** | `mf_newforms.traces` | level 1–1000 |
| F3a (3) / F3b (4) | EC over K₁ = Q(√5) `2.2.5.1` / K₂ = Q(√2) `2.2.8.1`, one row per isogeny class | `ec_nfcurves` (ainvs, `base_change`, `non_min_p`, `bad_primes`); a_P **computed by point counting** over residue fields | conductor norm 1–1000 |
| F4a (5) / F4b (6) | Hilbert modular forms over K₁ / K₂, parallel weight 2, **dimension 1** | `hmf_forms` + `hmf_hecke` (eigenvalues), prime order from `hmf_fields.primes` | level norm 1–1000 |

Genus-2 / abelian surfaces / paramodular forms are declared **out of scope**.

Columns: 25 per block. Q-blocks: rational primes 2…97. K-blocks: the first 25
prime ideals of the field in `hmf_fields.primes` order (norm-ascending, LMFDB
tie order), **shared between F3 and F4 of the same field**. Entries: raw
integer a_p / a_P. All six blocks have motivic weight 1, so the analytic
normalization exponent is identical across blocks and raw integers are already
commensurable; H1 verifies |entry| ≤ 2√(column norm) per cell. Markers (`?`)
at: bad-reduction primes, primes where the stored model is non-minimal
(`non_min_p`), ramified primes with no stored/derivable value, and any source
gap. Markers are never a reason to drop a row.

The combined table is the concatenation of blocks; pointer cost log₂ N_total;
per-(block, column) alphabets and coder contexts. This is an instantiation of
probe 02 §4–§5 (same KT baseline, same cost formulas, same cell code), not a
re-derivation; H2 verifies it against the probe-02 implementation directly.

**Seal.** EC/Q: probe 02's bands bind unchanged (this probe's F1 band and the
K4(a) extension stay inside the open 1–9999 train band). Each new block: rows
with conductor/level norm **901–1000 are sealed** (top decile of the declared
band): downloaded and stored, excluded from every scan, codebook, transport,
and measurement. Recorded in `SEAL.md` in this same commit.

## 2. Codebook construction — deterministic, no search, no minting

- **DUP within block**: exact-duplicate rows (expected none; one row per class).
- **DUP cross-block (the blind scan)**: allowed pairs, parent block first:
  F1→F2, F3a→F4a, F3b→F4b. Match = equal integer values at every shared good
  column, ≥ 10 shared good columns, exhaustive pairwise scan. Each match is a
  DUP(parent, child), accepted iff ΔE < 0. A child matching several parents
  takes the lowest parent index; the multiplicity is reported.
- **PATTERN**: for each (block, π ∈ Π, l ∈ {5,7,11,13}) with Π = {p+1, 0, 1, p}
  as in probe 02: members = rows whose every good cell conforms to π mod l;
  entry accepted iff ΔE < 0.
- **TEST**: F1 only, probe-02 Φ predicates (5-isogeny from pinned ecdata
  `allisog`), accepted iff ΔE < 0 and ΔE < ΔE(PATTERN on the same members).

Acceptance is ΔE < 0 under the probe-02 energy. There is no other criterion,
no proposal loop, and no entry type beyond these three.

## 3. Transports — meter, never codebook

Declared operation vocabulary, closed, v0 (each op maps an a-vector to a
predicted a-vector; predictions are exact integers; unpredictable columns are
skipped):

- `bc_K` (K ∈ {K₁, K₂}): split p (Kronecker (disc|p) = +1): a_P := a_p for each
  P | p; inert p: a_P := a_p² − 2p; ramified p: skip.
- `tw_χ` (χ ∈ {χ₅, χ₈}, Kronecker symbols of disc 5 and 8): a_p := χ(p)·a_p;
  p | disc(χ): skip.
- `mod` / `modK`: identity on values across a curated curve↔form link (F1↔F2
  via the LMFDB label link; F3↔F4 via the blind DUP match, cross-checked
  against `is_base_change`-consistent curation where available).

**Transport meter**: given a predicted value vector and an endpoint row, the
description cost E_X is the §4.3 cell code with exact-conform semantics (flag
+ nothing if predicted value equals the cell; escape + literal otherwise;
unpredicted or marker columns at baseline). Transport predictions are never
added to any codebook — the meter measures, it does not mint.

**Cycles** (two paths, same endpoints; instances = all eligible unsealed
starts, eligibility mechanical as stated):

- **C1 bc-square** (per field; known-commuting by theorem — the twist by the
  field character is trivial over K): start = F1 row E whose class has a base
  change in F3 (via `ec_nfcurves.base_change`); endpoint = that F3 row.
  Path A: `bc_K`. Path B: `tw_χ_K` then `bc_K`.
- **C2 mod out-and-back** (known-commuting trivially): start/endpoint = F1 row
  with a K1 link; Path A: identity; Path B: `mod` then `mod⁻¹`.
- **C3 modularity square** (per field): start = F1 row E with (i) a K1 form
  link g, (ii) a base-change row E_K in F3, (iii) E_K DUP-matched to an F4 row
  f. Endpoint = f. Path A: `bc_K` on E then `modK`. Path B: `mod` to g then
  `bc_K` on g's coefficients.
- **C4 twist-commute** (algebraically identical pointwise; reported as the
  triviality it is): `tw_χ₅∘tw_χ₈` vs `tw_χ₈∘tw_χ₅` on F1 rows, endpoint =
  the predicted vector itself (no row lookup).

**Closure defect** per instance, both reported: δ_bits = |E_A − E_B|;
δ_frac = fraction of columns predicted by both paths where the predictions
differ. The mathematics says C1–C4 commute: **a nonzero δ is the machine's
transport being lossy by default**, which is why the floor (N2) precedes any
reading and the concentration statistic, not the magnitude, is informative.

## 4. Cells — full matrix, one batch

### H — harness (any failure ⇒ VOID, published as such)

- **H1 join/parse/normalization.** Every row joins its source record exactly
  once; every integer entry satisfies |a| ≤ 2√(column norm); independent a_p
  recomputation for **200 random F1 rows (seed `default_rng(20260728)`)
  against the pinned ecdata aplist** (independent source, Cremona), joined by
  Cremona label: any mismatch at a shared good prime ⇒ VOID.
- **H2 energy reuse.** Probe 02's coder/energy imported unchanged; decoder
  test (byte-identical) on the combined table under ∅ and under the final
  codebook; emitted-vs-analytic ≤ 0.005; L₀ per probe-02 formulas. Drift ⇒ VOID.
- **H3 meter validity — load-bearing.** Every C1 and C2 instance must read
  δ ≤ floor (N2). If a known-commuting cycle reads above floor, the meter is
  measuring the machine; report that and stop. K5 is then uninterpretable and
  is not read.

### N — nulls, floor, power

- **N1 shuffled-family null.** Column shuffle within each block
  (`default_rng(20260728 + block_id)`, block ids 1–6 as tabled), rerun the
  cross-block DUP scan. Gate: **zero** cross-block matches.
- **N2 wobble floor.** Every cycle instance recomputed under **20 replicates**
  (`default_rng(20260728 + r)`, r = 0…19): replicate-permuted codebook order,
  column evaluation order, and endpoint tie-breaking. floor = q95 of δ per
  cycle type, fixed before K5 is read. (The transports are deterministic; a
  zero floor is a legitimate reading and binds K5 to "any nonzero δ is above
  floor".)
- **N3 power/foils.** Per cycle type: bc to the wrong field; tw by the wrong
  character (χ₅↔χ₈ swapped); mod against the next-labelled mismatched form.
  Gate: δ > floor in **≥ 95%** of foils. N3 is reported alongside every K5
  null; a K5 null without N3 passing is blindness, not agreement.

### K — known answers and measurements

- **K1 modularity, blind.** F1×F2 DUP scan (§2), scored afterwards against
  LMFDB's own link (lmfdb isogeny-class label N.x ↔ newform label N.2.a.x),
  unsealed band. Gate: **precision ≥ 99% and recall ≥ 95%**.
- **K2 base-change square.** All C1 instances, δ per instance. Expected at
  floor; feeds H3.
- **K3 degree-ladder headroom (report-only).** Bits saved as a fraction of
  L₀, per block, alongside probe 02's 1.67%.
- **K4 compounding, double dissociation.**
  (a) more rows, fixed family: EC/Q alone from **pinned ecdata**, bands
  N ≤ 1000 vs N ≤ 9999 (row ratio reported); codebook entries and bits saved,
  vocabulary of §2 only.
  (b) more families, fixed rows: blocks added in table order (F1, +F2, +F3a,
  +F3b, +F4a, +F4b) at this probe's bands; entries and bits saved after each.
  Interpretation fixed now: (a) flat & (b) growing ⇒ compounding scales with
  families; both growing ⇒ codebook tracks data volume; both flat ⇒
  saturation. **Ceiling caveat, pre-registered:** if K3 shows < 3% available
  compression in every family, a null in (b) is a ceiling reading, not
  evidence against compounding; RESULTS.md says which, mechanically, from K3.
- **K5 closure-defect census.** All declared cycles, all instances:
  distribution of δ (against per-type floor) and the **concentration
  statistic** = fraction of total δ_bits carried by the top 1% of instances.
  Fixed interpretation: at-floor throughout ⇒ transport closes (CLOSED-
  SYMMETRIC; N3 reported in the same sentence); above floor + top 1% < 10% ⇒
  uniform machine lossiness (LOSSY); above floor + top 1% > 50% ⇒ list the
  instances, no interpretation (CONCENTRATED). Intermediate concentration
  (10–50%): reported as LOSSY with the concentration number stated; the
  instance list is still published.

## 5. Ladder — mechanical, in order

1. **VOID** — any H fails. 2. **BLIND** — N3 fails. 3. **SUBSTRATE-FAIL** —
K1 misses either gate; K2–K5 not reported. 4. **CLOSED-SYMMETRIC** — K5 at
floor throughout. 5. **LOSSY** — above floor, low concentration.
6. **CONCENTRATED** — above floor, top 1% > 50%; instance list published, no
claims. K3 and K4 are reported under every verdict from 4 onward. Residuals
are judged against LMFDB's link tables, not our codebook: anything already
linked on the object's LMFDB page is rediscovery and is labelled as such.

No parameter changes after this commit. No follow-up round. "Inconclusive at
this scale with this method" is terminal and recorded CLOSED.

## 6. What no verdict here establishes

No new mathematics, no conjecture, no discovery claim, nothing about the
sealed bands or the paramodular zone. VALIDATED-shaped outcomes establish only
that the meter, the substrate, and the known correspondences are mutually
consistent at this scale; CONCENTRATED publishes a list of instances and
nothing else.
