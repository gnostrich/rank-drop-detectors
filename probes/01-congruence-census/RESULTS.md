# RESULTS — Probe 01: congruence census — **VALIDATED**

Run 2026-07-28, one round, terminal. Data: ecdata `25cec5ec`, conductors
1–9999. Prereg: [PREREG.md](PREREG.md), committed before the run. Full
artifacts in [`output/`](output/).

## Usable data (reported before any hit rate)

**38,042 isogeny classes (64,687 curves) — 100% usable, 0 discarded.**
Matching used good positions only; bad-prime markers were skipped, never a
reason to drop a class. Only 315 classes have no bad prime below 97, so a
reader that discards on bad primes loses >99% of the corpus. Primes parsed:
25 (p = 2..97); used per modulus: 24 (column p = l excluded). Minimum shared
good positions on any accepted pair: 17 (preregistered floor: 10); pairs
ineligible under the floor: 0.

## The four cells

**H1 — harness: PASS.** Join integrity clean (every curve ↔ one aplist row,
no duplicates); 0 isogeny classes split across fingerprint groups; independent
point-count recomputation of a_p for 200 random classes at all good p ≤ 97:
**0 mismatches**. Cross-class exact-agreement pairs: 0 (prior run also 0).

**H2 — null calibration: clean.** Marginal-preserving column shuffle, rerun
identically: **0 non-singleton components for every l ∈ {5, 7, 11, 13}**.
Empirical false-positive rate at this depth: 0.

**A1 — census** (non-trivial congruence groups after isogeny is collapsed):

| l  | groups | classes involved | largest group sizes |
|----|--------|------------------|---------------------|
| 5  | 3,521  | 8,441            | 143, 64, 53, 37, 22 |
| 7  | 896    | 1,906            | 17, 11, 11, 6, 6    |
| 11 | 32     | 64               | all pairs           |
| 13 | 5      | 10               | all pairs           |

**A2 + A3 — structure and depth.** All group members' a_p extended by point
counting to every good p < 1000. Of 22,598 pairs evaluated, **22,593 are
persistent** (congruent at every tested good prime) and **100% are
pattern-consistent** (threshold: 90%): EIS 10,289, P1 (equal conductor)
2,219, P2 (one conductor divides the other) 3,472, P3 (common factor) 6,613,
**P4 (coprime, non-Eisenstein): 0**. The largest mod-5 group (143 classes,
containing 11a, 38b, 50b, …) is entirely Eisenstein-type (a_p ≡ p + 1 mod 5)
— the known reducible family, where coprime conductors are expected. Every
mod-11 and mod-13 hit is a same-level or divisibility pair with prime
conductor ratio (mod 11: ratios 1, 7, 13, 17, 19; mod 13: ratios 13, 19) —
the level-raising shape. **Depth:** 5 pairs broke in the extension, all in
one mod-5 group {5408a, 7200bg, 7200br, 9248e}: agreement through the aplist
columns, breaking at p = 13/101/113 — finite-depth coincidences, recorded in
`output/breaking_pairs.csv` with no further interpretation.

## Ladder verdict (mechanical)

Not VOID (H1 passed; H2 = 0 everywhere). Not KNOWN-ONLY (4,454 non-trivial
groups). Persistent pairs exist; consistency 1.00 ≥ 0.90 → **VALIDATED**:
the detector recovers known congruence structure — Eisenstein families and
level-raising/lowering conductor patterns — from raw a_p data alone.

## What is NOT claimed

No new mathematics and no conjectures. Trace congruence at finitely many
primes (24 columns + 143 extension primes) is evidence for, not proof of,
isomorphic semisimplified mod-l representations. EIS classification is a
trace-pattern signature, not a computed torsion/isogeny certificate.
Conductor-pattern consistency is descriptive; no instance of Ribet's theorem
was verified. The result establishes exactly one thing: a rank-drop detector,
run as a systematic census over an arithmetic family, recovers known
structure — the precondition for pointing it anywhere unknown, and nothing
more.
