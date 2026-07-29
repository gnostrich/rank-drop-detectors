# RESULTS — Probe 05: minting on a non-duplicate substrate — **NO-MINT**

Run 2026-07-29, one round, terminal. Two ledgers, kept separate as the
directive requires:

**Objective ledger (the one that counts): FAILURE.** K3 = 0. No anchor was
minted that is accepted, reused, and not identifiable with an existing
object. Per the ladder, NO-MINT **closes the conjecture-engine claim on
this substrate**, and the branch closes with it. Nothing below changes
this.

**Instrument ledger (a precondition, not progress):** every harness and
null cell passed; the numbers are recorded because later preregs cite
them, not because they advance the objective.

## Data

S-A: 38,042 EC/Q classes (pinned ecdata). S-B: **954 usable Artin reps**
(of 2,781 enumerated at conductor ≤ 2000: **1,521 (55%) have no L-function
in LMFDB** — the acquisition fact of this probe; 306 sealed) + 2,463 EC/Q
classes + 5,092 newforms from probe-03 snapshots. Manifest: 3,885 requests.

## The matrix

| cell | gate | result | verdict |
|---|---|---|---|
| H1 recomputation | any mismatch ⇒ VOID | SA 0/200 rows; Artin bound 0 violations full-table; factor-membership 0/200 | PASS |
| H2 coder + minted rows | byte-identical | all five streams incl. the **MINT block** decode identically | PASS |
| H3 reduction to HUB | exact | reuse-stripped ANCHOR reproduces probe-02's k4_groups.csv **exactly**; in-band 0.351541 = 0.351541 | PASS |
| N1 shuffled, full path | any accept ⇒ FATAL | **0** (S-A), **0** (S-B) | PASS |
| N2 wobble floor | report always | 20 orderings: **Jaccard = 1.00 across all pairs; bits spread q95 = 2×10⁻¹² ⇒ floor ≈ 0**. The acceptance pass is order-insensitive in practice despite live ownership coupling | set |
| N3 foils, full path | any accept ⇒ BLIND | 0/2,000 | PASS |
| K1 known-answer mint | ΔE < 0, accepted | EIS-5 mints: k = 143, **ΔE = −2,752**; reuse = 0 (reported ungated — weak on a single-family substrate by design) | PASS |
| K2 census (S-B) | report | **111 anchors accepted**; every one reuse = 0; every one identifiable (`row:` shadows) | reported |
| K3 objective | count | **0** | empty |
| K4 compounding | report always | chain depth 0 for all 20 S-A and all 111 S-B anchors — **no anchor was built on** | depth 0 |
| K5 power | ≥500 on S-B | 4,454 (S-A) and **8,288 (S-B) candidates — ADEQUATE** | met |

## Disclosed defect (registry item for any successor prereg)

Identifiability rule (a) as committed is satisfied by an anchor's own
members (each member's good columns reduce to r_h by construction), so
**K3 = 0 is structurally guaranteed under this prereg**, independent of the
data. Implemented verbatim per the no-parameter-changes constraint; the
honest reading of this probe's K3 is therefore "empty by definition", and
the informative cells are K2 (111 anchors, all reuse 0), K4 (all depths 0),
and N2 (floor ≈ 0). A successor prereg needs an identifiability definition
that a member cannot satisfy — e.g. requiring the identifying object to be
good at every comparison column.

Also flagged as required: reuse is the repo's first non-single-bit
acceptance rule; it is a downstream ΔE in the same currency (authority
deferred, not split).

## What is NOT claimed

Instrument validation above is **not** progress toward the objective — the
probe validates cleanly and mints nothing, which the directive defines as
failure against the objective, and that is this probe's primary result.
Zero reuse and zero chain depth mean the compounding hypothesis received no
support here. The 55% Artin L-function gap is a data-infrastructure fact,
not mathematics. Sealed bands untouched; no sealed-band confrontation was
earned by any anchor.
