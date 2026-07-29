# PREREG — LOOP v0: the promotion cycle

Committed before any run. One round, terminal, same ladder discipline as
probes 02–05. Build per the LOOP design document: inner loop on held-out
primes (cheap, unlimited), outer confrontation on held-out rows (scarce) —
**the confronter is specified and NOT run in this build**; the sealed bands
are untouched by every cell. Validation substrate: S-A only (EC/Q,
conductor < 10000, pinned ecdata, probe-02 energy unchanged, probe-05
vocabulary). Expected outcome, recorded now: **the known-answer anchor
promotes, shuffled data promotes nothing, and chain depth stays at 0 — the
loop works and does not compound.** A depth-0 result is not a bug.

## 1. Constants (v0, locked)

| constant | value |
|---|---|
| primes | 25 (p = 2…97) |
| folds F | 5; each hides a uniform-random 10 of 25 (`default_rng(20260728 + fold)`) |
| P_fit / P_hold | 15 / 10 |
| θ | 10⁻⁶, Bonferroni over (anchor, fold) scores in the round; never moved |
| R_max | 5 rounds; stop early on a round with zero promotions |
| moduli | {5, 7, 11, 13} |
| candidate source | probe-05 clusters (probe-01 components, k ≥ 2), computed on **fit columns only** (masked array; held columns unreadable at fit time) |
| quarantine Λ | member-set hash per modulus, append-only (`output/quarantine.jsonl`); never re-proposed at that modulus |
| N2 orderings | 20 (`default_rng(20260728 + r)`) |

## 2. Fit / predict / score / promote

- **Fit.** Residue vector on fit columns from the lowest-indexed member good
  there. Anchor type: **closed-form** iff the fit-column residues equal
  π(p) mod l for some π ∈ Π = {p+1, 0, 1, p} at every fit column; else
  **free**. Promotion rates are reported per type (the design predicts
  closed-form is stronger).
- **Predict.** Held prime p: closed-form anchors predict π(p) mod l for
  every member. Free anchors predict the reference member's (lowest index)
  held value mod l for the **other** members; the reference's own held
  cells are excluded from the trial count (nothing else in the fit read
  them).
- **Score.** Trials n = held-good member cells (bad positions excluded,
  never imputed); hits ĉ; p_tail = upper tail of Binom(n, 1/l). Promote iff
  p_tail < θ / (#scores this round). Kill otherwise → Λ.
- **Promote.** The anchor becomes a row of T in the MINT family
  (probe-05 uniform coding, norms 43), fingerprint = residue vector over
  all 25 columns (reference values at previously held columns); total
  description bits recomputed; next round's descent sees the mint rows, so
  chained anchors (chain depth ≥ 1) are possible.
- An anchor promoted in any fold is promoted once; scored-and-never-passing
  member sets are quarantined.

## 3. Cells and ladder — mechanical, in order

1. **H (VOID if any fails):** H-a masked-fit enforcement (the fit function
   receives an array with held columns structurally absent; verified by
   construction test); H-b probe-02 coder byte-identical on S-A + MINT
   rows; H-c scorer calibration: 10⁶ simulated null trials per l give
   empirical tail within 2× of Binom at p = 10⁻⁴.
2. **FATAL:** the full loop on marginal-preserving shuffled S-A
   (`default_rng(20260728 + 1)`) promotes **anything**, in any round. The
   shuffled promotion curve (per round) is reported alongside the real one
   under every verdict.
3. **SCORER-DEAD:** the mod-5 Eisenstein anchor fails to promote on real
   data (it predicts held primes exactly; failure means the scorer is
   wrong, not the mathematics).
4. **NO-COMPOUND:** loop runs, promotes, but chain depth = 0 in all rounds.
   Terminal. (Expected.)
5. **SATURATED:** total description bits fail to decrease round over round
   while promotions continue. Terminal.
6. **COMPOUNDING:** chain depth ≥ 1 achieved off non-shuffled data with all
   H/N cells clean. Report depth distribution, promoted list, N2 spread; no
   interpretation.

Under every verdict from 3 on: promoted-set size per round, L₀ trajectory,
chain-depth histogram, per-type promotion rates, N2 order spread (Jaccard
over 20 orderings of the final promoted set + bits spread, floor = q95),
and `loop.jsonl` per the design document's schema.

## 4. Not claimed

Promotion is a statistical statement about held-out columns on fitted rows
— not mathematics. No sealed row is read; no conjecture is published; the
confronter exists in code and is not exercised. A working loop with depth 0
establishes exactly: the promotion cycle functions and does not compound on
this substrate at this scale.
