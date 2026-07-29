# RESULTS — LOOP v0: the promotion cycle — **SATURATED**

Run 2026-07-29, one round, terminal. S-A only (EC/Q, 38,042 classes, pinned
ecdata); the confronter was **not** run; no sealed row was read. The
expected outcome ("promotes known answer, shuffled promotes nothing, depth
0") was **half right** — and the half that failed is the informative part.

## The matrix

| cell | gate | result | verdict |
|---|---|---|---|
| H-a masked fit | held columns unreadable | mutation canary: fit residues unchanged | PASS |
| H-b coder | byte-identical | S-A stream identical | PASS |
| H-c scorer calibration | within 2× at 10⁻³/10⁻⁴ | 10⁶ simulated trials per modulus: within band | PASS |
| FATAL (shuffled loop) | any promotion | 75 scored, **0 promoted**; curve flat at 0 | PASS |
| SCORER-DEAD (known answer) | EIS-5 must promote | promoted (closed-form, p_tail ≈ 0) | PASS |
| SATURATED (bits trajectory) | must decrease | **3,929,900 → 4,114,692 → 4,338,438 → 4,666,755 → 5,296,037** — rising every round, and above the 3,866,992 baseline from round 1 | **FIRES** |

## What the loop actually did

Promotions per round: **765, 824, 976, 1,388, 2,584** (6,537 total; scored
22,523 then 3,830–7,010/round; quarantined 4,298). Chain-depth histogram:
depth 0: 765 · 1: 824 · 2: 976 · 3: 1,388 · 4: 2,584. **Promotion rate by
type: closed-form 125/126 (99.2%) vs free 22,843/42,273 (54.0%)** — the
design document's prediction that closed-form patterns are stronger is
confirmed emphatically. **N2 order spread: Jaccard mean 0.074 (min 0.068)
across 20 orderings; bits-spread floor (q95) = 862 bits** — the first live
antisymmetric measurement in the repo, and it is large: the promoted set is
~93% order-noise.

## Mechanism, disclosed (registry defects for any v1)

1. **Mint-echo cascade.** A promoted anchor's row re-clusters with its own
   members next round under a fresh member-set hash; the echo promotes
   (it inherits the original's predictive success), spawning depth d+1
   echoes of depth-d anchors. The rising promotion counts and the entire
   depth ≥ 1 population are echoes, not compounding. v0 has no echo
   exclusion; the quarantine is keyed on exact member sets and cannot catch
   supersets-by-one-mint-row.
2. **Promotion is decoupled from the energy.** The promoter's criterion is
   predictive significance only; most promoted free anchors have positive
   full-column ΔE (pointer + row costs exceed savings), so total
   description **rises from round 1** even before echoes. The single
   authority of probes 02–05 (ΔE) and the loop's statistical filter
   genuinely disagree, and the bits trajectory is where that disagreement
   is visible.
3. **Order dependence at Jaccard 0.07** is the quantitative cost of 1+2:
   which echoes exist depends on processing order.

The ladder never reaches COMPOUNDING: depth ≥ 1 exists, but SATURATED
(bits non-decreasing while promotions continue) fires first, and per the
design document that is terminal. The depth histogram is reported above as
required, with the echo mechanism attached.

## What is NOT claimed

Promotion is a statistical statement about held-out columns on fitted rows
— not mathematics; no conjecture was generated, and no promoted object
earned a sealed-band confrontation. The 99%-vs-54% closed/free split and
the flat shuffled curve validate the scorer, not the objective. The loop as
specified **functions and does not compound in any meaningful sense**: what
compounds is its own bookkeeping, at increasing description cost, in an
order-dependent way. That is the honest v0 result, and any v1 needs echo
exclusion and an energy-coupled promotion rule before the sealed band is
worth touching.
