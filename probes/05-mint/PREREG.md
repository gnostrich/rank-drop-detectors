# PREREG — Probe 05: minting, on a non-duplicate substrate

Committed before any scan, mint, or measurement. One round, terminal. **This
probe is on the critical path to the objective (a conjecture engine); probes
01–04 validated the instrument and are a precondition, not progress.**
RESULTS.md keeps the two ledgers separate: instrument validation is not
objective progress, and a probe that validates cleanly and mints nothing is
a failure against the objective. Expected outcome, recorded now: **NO-MINT
or UNDERPOWERED**; the objective-relevant probability is concentrated in K3
and K3 is expected empty.

Two variables change (move vocabulary + substrate), decomposed by running
the move on both substrates in one batch: S-A isolates the move, S-B the
territory. BRIDGE and HUB stay dead. If the reuse condition is relaxed
mid-run the probe is VOID.

## 1. The move — ANCHOR(h, S, l)

Same energy arithmetic as the killed HUB (probe 02 §5: entry cost
H + Γ(k) + k·log₂N + C_l·log₂ l, member cells claimed at residues r_h),
**plus a minted row and a second acceptance condition**:

```
accept  ⟺  ΔE < 0   AND   reuse(h) ≥ 1
```

- The minted row: h enters the table as a row in a dedicated MINT block,
  fingerprint = its residue vector r_h (values 0…l−1), all 25 columns.
  MINT-block coding is the uniform code per column (alphabet from norm 43,
  so bound 13 ≥ l−1 for every declared modulus; model cost 0; cell cost
  log₂|A_c|). The row's full coding cost is charged inside ΔE. Pointer
  costs use log₂ N₀ (initial open rows), frozen at table build; disclosed.
- reuse(h) ≥ 1 iff, after all ΔE-accepted anchors are minted, h participates
  in at least one further accepted (ΔE < 0) entry from the fixed menu:
  DUP between MINT rows; PATTERN over the MINT block (Π, moduli as probe
  02); ANCHOR over clusters of MINT rows (mod-l congruence scan among MINT
  rows, same rules — these are chained anchors, K4's depth). No other reuse
  channel exists. Reuse is a downstream ΔE in the same currency: authority
  deferred, not split (flagged in RESULTS as the first non-single-bit
  acceptance rule in the repo).
- Candidate clusters: connected components of the probe-01 congruence
  relation (mod l ∈ {5,7,11,13}, shared good ≥ 10, p = l column excluded on
  rational-prime blocks) over the substrate's full row set, k ≥ 2. One
  ANCHOR candidate per (component, l); r_h[c] = residue of the lowest row
  good at c, else 0.
- Order sensitivity is live (cross-modulus cell ownership + reuse pass):
  the acceptance pass processes candidates in a declared order (by l, then
  −k, then lowest row); N2 permutes this order.

## 2. Substrates

- **S-A (known answer):** probe 02's EC/Q class table, conductor < 10000,
  from pinned ecdata, unchanged (38,042 rows; probe-02 seal bands bind).
- **S-B (non-duplicate):** three blocks — Artin representations (LMFDB
  `artin_reps`: Dim 2, CharacterField 1 (rational), conductor ≤ 2000;
  a_p from `lfunc_instances` → `lfunc_lfunctions.euler_factors`, the
  two-step route verified on 2.23.3t2.b against the known weight-1 form;
  reps without an L-function row are counted and excluded); probe 03's F1
  (EC/Q ≤ 1000) and F2 (weight-2 newforms, R1 reading) snapshots reused
  byte-identical. Artin block norms = 1 per column (|a_p| ≤ 2, motivic
  weight 0; commensurability across weights is NOT asserted — congruence
  mod l is the only relation used). Markers at bad/hard primes and any
  missing coefficient.
- Pin/seal: committed snapshots + SHA256 manifest, unversioned-API caveat
  disclosed in the same words as probes 03–04. Artin top decile (conductor
  1801–2000) sealed, never scanned. All prior seals bind.

## 3. Cells

### H (any failure ⇒ VOID)
- **H1.** 200 random rows per family (seed 20260728): S-A rows re-counted
  by point counting vs ecdata (probe-01 machinery); Artin rows' a_p
  cross-checked against character theory where determinable (the identity
  class bound |a_p| ≤ 2 at every good prime for every row — full-table
  check — plus 200-row recomputation of a_p parity/values against
  `GaloisConjugates.Character` where the local factor list pins the value:
  unramified p with all distinct factors sharing... if not mechanically
  determinable for a row, the row is skipped and counted). Any mismatch ⇒
  VOID.
- **H2.** Probe-02 coder unchanged; byte-identical decode of both
  substrates **including tables containing minted anchor rows** (MINT block
  streams with the uniform code); accounting ≤ 0.005.
- **H3.** With the reuse condition removed and the minted-row charge
  removed, ANCHOR's ΔE arithmetic must equal probe 02's HUB ΔE on the
  identical 714 k≥3 groups, reproducing A_HUB in-band rate = probe 02's
  measured 35.154% exactly (same data, same formula ⇒ exact equality is
  the gate; any deviation means the move is silently different ⇒ VOID).

### N
- **N1.** Full accept path (ΔE + reuse) on marginal-preserving shuffled
  tables (`default_rng(20260728 + fid)`, fid: S-A = 1; Artin/F1/F2 = 2/3/4).
  **Any anchor accepted on shuffled data ⇒ FATAL.**
- **N2 wobble floor, genuine.** 20 replicates on S-A
  (`default_rng(20260728 + r)`): candidate order permuted; full accept +
  reuse pass rerun. Report the Jaccard distribution over final accepted
  anchor sets and the bits spread; floor = q95 of bits spread, fixed before
  K cells are read. Reported under every verdict.
- **N3.** 1000 foils per substrate (seed 20260728): random member sets of
  matched size distribution, random modulus, random residue vector, run
  through the FULL accept path including reuse. **Any accepted ⇒ BLIND.**

### K
- **K1 known-answer mint (S-A).** The mod-5 Eisenstein cluster (the
  143-member component) must yield an accepted-on-ΔE anchor: gate
  ΔE < 0 and the anchor mints. reuse reported **without gating** (weak on a
  single-family substrate by design).
- **K2 mint census (S-B).** Every anchor accepted on S-B: count; per anchor
  k, ΔE, reuse, and LMFDB-identifiability. **Identifiable** ⟺ (a) some
  existing row in any S-B block satisfies a_p ≡ r_h[c] (mod l) at every one
  of its good columns (≥10 shared) — the anchor is that object's mod-l
  shadow; or (b) r_h equals a closed-form Π pattern (p+1, 0, 1, p mod l) at
  every column — a named Eisenstein-type object. Both labelled rediscovery.
- **K3 the objective cell.** Anchors accepted ∧ reuse ≥ 1 ∧ not identifiable.
  The count is the probe's answer; if non-zero, one sentence each, no
  interpretation.
- **K4 compounding.** Chain-depth distribution (anchors over MINT rows =
  depth ≥ 1). First honest test; probe 03's version is registered as
  confounded.
- **K5 power.** Rows, primes, candidate clusters per substrate.
  **< 500 candidates evaluated on S-B ⇒ a K3 zero is SAMPLE-SIZE.**

## 4. Ladder — mechanical, in order

1. **VOID** — any H fails (or reuse condition relaxed mid-run).
2. **FATAL** — N1 accepts any anchor on shuffled data.
3. **BLIND** — N3 accepts any foil.
4. **MOVE-DEAD** — K1 fails; closes minting on this energy.
5. **UNDERPOWERED** — K1 passes, K3 empty, K5 < 500.
6. **NO-MINT** — K1 passes, K3 empty, K5 ≥ 500. Closes the conjecture-engine
   claim on this substrate; the branch closes with it.
7. **MINTED** — K3 non-empty: list + N2 floor + K4 depths, no
   interpretation; sealed-band confrontations are separate registered
   actions.

N2 and K4 reported under every verdict from 4 onward. No parameter changes
after this commit; no follow-up round; sealed bands untouched; anything
identifiable with an existing LMFDB object is rediscovery, labelled.

## 5. What no verdict here establishes

MINTED items are statistical candidates, not mathematics — they earn one
registered sealed-band confrontation each, later, separately. NO-MINT
closes the conjecture-engine claim on this substrate only. And in every
case: instrument validation reported here is not progress toward the
objective.
