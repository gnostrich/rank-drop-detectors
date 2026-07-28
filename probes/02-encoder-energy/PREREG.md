# PREREG — Probe 02: the encoder and its energy

Committed before any encoder runs. Every definition, constant, threshold and seed
below is fixed here; `run.py` implements this document mechanically and the verdict
is read off the ladder in §9 with no discretion. One round, terminal.

Builds past: `probes/01-congruence-census` (VALIDATED, REGISTRY 2026-07-28) and the
kill ledger in `gnostrich/mathlib-structure-probes` — in particular **no MDL judge
is pointed at any Mathlib object here; the compression judge lives on arithmetic
data only**, and **no Gibbs move is made in formal-statement space**.

---

## 0. Where this sits in the machine

The machine is two parts (design notes, `NOTES-2026-07-28-generator-design.md`):
a **generator** (Gibbs process on fingerprint tables, all scientific content) and a
**notary** (stateless Lean exporter, adds trust not truth, one-way pipe).

The generator has five pieces: substrate, energy, proposal distribution, moves,
schedule. **Probe 02 tests pieces 1–2 only: substrate and energy.** No search, no
descent, no temperature, no proposals. Codebooks are hand-supplied; the only thing
under test is whether the energy is a real code length, is not gameable, and pays
the right number of bits for the moves we already know the answers to.

Probe 03 gets pieces 3–5 (proposals, moves, descent) and only inherits move types
that pass a K-cell here. A move type that fails K is deleted from probe 03's
vocabulary; it is not repaired and re-run.

Rationale for the split: if the energy is gameable, every later result is an
artifact of the energy, and no amount of descent quality can recover it. The energy
is the single authority; it gets audited alone, before anything is allowed to
optimise against it.

---

## 1. What each component does, in plain English

- **Substrate.** A fixed table of raw numbers. One row per arithmetic object, one
  column per prime, the entry being that object's `a_p`. Nothing derived, nothing
  computed by us, pinned to an external commit so we cannot have sculpted it.
- **Baseline model.** The dumbest honest description of the table: for each column
  separately, a histogram of how often each value occurs. This is what "no
  structure found" costs.
- **Codebook.** A list of claims about the table — "this row matches that one mod
  5", "these 143 rows all satisfy `a_p ≡ p+1` mod 5". Each claim is a *hint to the
  encoder*, never a constraint: the encoder stays a valid encoder even if the claim
  is false, because it can emit an escape symbol and the literal value. A false
  claim therefore costs bits rather than breaking anything. This is why
  falsification is not a module — it is the arithmetic of the escape symbol.
- **Energy.** The total number of bits: the codebook itself, plus the table given
  the codebook. A claim is accepted only if writing it down costs fewer bits than
  it saves. Nothing else adjudicates.
- **Decoder.** The thing that proves the energy is a code length and not a score.
  If the bits do not decode back to the exact table, we are optimising a
  scoreboard, and scoreboards are gameable.

---

## 2. Notation

| symbol | type | shape / range | rate | locus |
|---|---|---|---|---|
| `P` | prime list | `p = 2 … 97`, `C = 25` | fixed | substrate |
| `T_cls` | table | `38042 × 25` | fixed, append-only | substrate A (isogeny classes) |
| `T_crv` | table | `64687 × 25` | fixed, append-only | substrate B (curves; K1 only) |
| `x[r,c]` | cell | `ℤ ∪ {+,−,?}` | fixed | substrate |
| `A_c` | column alphabet | `{+,−,?} ∪ [−⌊2√p_c⌋, ⌊2√p_c⌋]` | fixed | encoder |
| `q_c` | column marginal | prob. vector on `A_c` | recomputed once, at baseline | encoder |
| `N` | rows | 38042 (A) / 64687 (B) | fixed | substrate |
| `𝒦` | codebook | list of entries | hand-supplied in this probe | state |
| `e` | entry | one of DUP / BRIDGE / PATTERN / HUB / TEST | — | state |
| `l` | modulus | `{5, 7, 11, 13}` | fixed | entry field |
| `r_h` | hub residue vector | `(ℤ/l)^C` | — | entry field |
| `S` | member set | subset of rows, `k = |S|` | — | entry field |
| `ε₀` | escape probability | `2^-6` (v0) | prereg-locked | encoder |
| `ρ` | rent per unproved assumption | `0` bits (v0, inert here) | declared | energy |
| `L₀` | baseline code length | bits | — | energy |
| `E(𝒦)` | energy | bits | — | energy |
| `ΔE(e)` | entry delta | bits, negative = accepted | — | energy |
| `S_i` | per-row saving | bits | measured from `q_c` | energy |

Units are bits throughout. All costs are real-valued (arithmetic coding); **no
ceilings anywhere** — a `⌈·⌉` in the implementation is a bug and H3 will catch it.

---

## 3. Factor graph

```
        column marginals  q_1 … q_25          (baseline factors, one per column)
              │    │    │        │
              ▼    ▼    ▼        ▼
   rows  ┌──────────────────────────────┐
    r=1  │  x[1,1] x[1,2] … x[1,25]     │
    r=2  │  x[2,1] x[2,2] … x[2,25]     │  ← cells: the only thing encoded
     ⋮   │             ⋮                │
    r=N  │  x[N,1] x[N,2] … x[N,25]     │
         └──────────────────────────────┘
              ▲              ▲
              │              │
        ┌─────┴────┐   ┌─────┴──────────────┐
        │ BRIDGE e │   │ PATTERN / HUB e    │   (codebook factors)
        │ i → j, l │   │ S, l, π  or  S,r_h │
        └──────────┘   └────────────────────┘
              ▲                   ▲
              │                   │
        (2 row pointers)   (k row pointers, one per member)
```

Every codebook factor attaches to a set of **cells** and does exactly one thing:
it narrows the alphabet at those cells from `A_c` to the residue class it names,
at the price of its own description plus an escape channel per attached cell.

The single structural fact that the whole entry-type hierarchy turns on: **a BRIDGE
names two rows; a PATTERN/HUB membership names one.** Everything in §5's arithmetic
is downstream of that.

---

## 4. The encoder, exactly

### 4.1 Canonical table

`T_canon` = the parse of the pinned `aplist` / `allcurves` files into a
`N × 25` array over `A_c`, serialised canonically (one row per line, entries
space-separated, no leading zeros, bad markers as the literal `+ − ?`).
`parse` is pure and total.

### 4.2 Baseline (empty codebook)

Column marginals `q_c` are plug-in estimates over `A_c` with Laplace smoothing
`α = 1/2` (Krichevsky–Trofimov), transmitted as part of the model:

```
L₀ = Σ_c [ model_c ]  +  Σ_{r,c} −log₂ q_c(x[r,c])
model_c = (|A_c| − 1)/2 · log₂ N          (KT / Rissanen parametric cost)
```

`L₀` is the number every codebook is scored against. `E(∅) = L₀`, `ρ` inert.

### 4.3 Cells under an entry

Let entry `e` claim, at cell `(r,c)`, residue `ρ_e(r,c) ∈ ℤ/l`. Write
`Q_c(ρ) = Σ_{y ∈ A_c, y ≡ ρ (mod l)} q_c(y)` (bad markers `+,−,?` are **never**
claimed by any entry; they always fall through to baseline).

The cell is coded in two parts:

```
conform:  cost = −log₂(1−ε₀)  −log₂ q_c(x)  + log₂ Q_c(ρ_e(r,c))     if x ≡ ρ_e
escape:   cost = −log₂ ε₀     −log₂ q_c(x)                            otherwise
```

with `ε₀ = 2^-6`, so conforming costs `0.0227` bits of flag and escaping costs
`6` bits of flag. Per-cell saving when conforming:
`s(r,c) = −log₂ Q_c(ρ_e(r,c)) − 0.0227`.
Per-row saving `S_r = Σ_{c good, c ≠ column(l)} s(r,c)`.

**A false claim is not an error. It is 6 bits.** That is the entire falsification
mechanism, and it is why no separate refutation module exists.

Overlapping entries: cells are claimed by at most one entry; if two entries claim a
cell, the one **earlier in codebook order** owns it. Codebook order is part of the
codebook and is transmitted implicitly by the list order. No cell is double-counted.

### 4.4 Energy

```
E(𝒦) = L₀_model + Σ_e cost(e) + Σ_{r,c} cell_cost(r,c | 𝒦) + ρ·(#unproved)
ΔE(e) = E(𝒦 ∪ {e}) − E(𝒦)
```

Entry accepted iff `ΔE(e) < 0`. There is no other criterion anywhere in the
machine. `ρ = 0` in v0 and every entry here is verified on the encoded range, so
the rent term is declared-and-inert; it becomes live only in probes that assert
beyond-range predictions.

---

## 5. Entry vocabulary and its exact arithmetic

Header for every entry: `type` (3 bits, 5 declared types incl. one reserved) +
`modulus` (2 bits, 4 declared moduli) = **`H = 5` bits**. Row pointer =
`log₂ N` bits (`15.215` for substrate A, `15.981` for B). Member-list length prefix
= Elias-gamma `Γ(k) = 2⌊log₂ k⌋ + 1`.

| type | fields | `cost(e)` | claims |
|---|---|---|---|
| `DUP(i→j)` | 2 pointers | `H + 2 log₂ N` | every cell of `j` equals `i` (exact, `l` field unused) |
| `BRIDGE(i→j, l)` | 2 pointers, `l` | `H + 2 log₂ N` | `a_p(j) ≡ a_p(i) mod l` |
| `PATTERN(S, l, π)` | member list, `π ∈ Π` | `H + 2 + Γ(k) + k log₂ N` | `a_p ≡ π(p) mod l` for all of `S` |
| `HUB(h, S, l)` | member list, free residue vector | `H + Γ(k) + k log₂ N + C_l·log₂ l` | `a_p ≡ r_h[c] mod l` for all of `S` |
| `TEST(φ, l, π, X)` | predicate id, exception list | `H + 2 + log₂|Φ| + Γ(|X|) + |X| log₂ N` | membership = `φ`, corrected by `X` |

`Π` (closed-form patterns, `|Π| = 4`, 2 bits): `π_EIS: p ↦ p+1`, `π_0: p ↦ 0`,
`π_1: p ↦ 1`, `π_id: p ↦ p`.
`Φ` (predicates, external, pinned): `has-rational-l-isogeny` (from ecdata isogeny
matrices at the pinned commit), `N ≡ 0 mod m` for `m ≤ 12`, `has-CM`. `|Φ| = 14`.
`C_l = 24` (the `p = l` column is excluded, as in probe 01 §4).

### 5.1 The three analytic predictions this probe exists to test

These are derived from §4–§5 and are what the K-cells check. They are stated with
their formulas, not just their numbers, and are evaluated at run time with the
measured `S_r`.

**(P-a) A bridge pays for itself iff a row's saving exceeds two pointers.**
`ΔE(BRIDGE) = H + 2 log₂ N − S_j`. With `l = 5` and equidistributed residues,
`S_j ≈ 24 log₂ 5 = 55.7` against `5 + 30.43 = 35.4` ⇒ `≈ −20.3` bits.
For `l = 13`: `S_j ≈ 24 log₂ 13 = 88.8` ⇒ `≈ −53.4` bits. Bridges are always
accepted at every declared modulus; a bridge that fails to compress means the
measured residues are far from equidistributed and the prediction is wrong.

**(P-b) PATTERN beats a spanning tree of BRIDGEs by the pointer economy alone.**
The honest comparator for a `k`-member group is **not** `k(k−1)/2` pairwise bridges
(that is raw detector output, which MDL prunes) but the minimum spanning tree of
`k−1` bridges. Against that:

```
A_PAT(k) = ΔE(PATTERN) − ΔE(MST) = 2 + Γ(k) + (2−k)·log₂ N − 5(k−1) − S_root + S_root
         = 2 + Γ(k) + (2−k)·log₂ N − 5(k−1)
```

(the root row's saving cancels: PATTERN compresses it, MST does not, and PATTERN's
pattern id is O(1) where MST would need it too). At `k = 143`, substrate A:
`A_PAT = 2 + 15 − 141·15.215 − 710 = **−2 838 bits**`, against a total group
`ΔE(PATTERN) ≈ −5 700` bits. So roughly **half** of the 143-family's compression is
the pointer economy of not naming a partner.

**(P-c) A HUB's own residue vector costs exactly what the spanning tree's root
fails to save, so hub-vs-tree reduces to pointer economy and nothing else.**

```
ΔE(HUB) = H + Γ(k) + k log₂ N + C_l log₂ l − Σ_{i∈S} S_i
ΔE(MST) = (k−1)(H + 2 log₂ N) − Σ_{i∈S∖root} S_i
A_HUB(k) = Γ(k) + (2−k)·log₂ N − 5(k−1)          ← C_l log₂ l and S_root cancel
```

The cancellation is exact and is the sharpest content in this document: the hub's
free residue vector `r_h` costs `C_l log₂ l`, and the root row of the spanning tree
forgoes exactly `S_root = C_l log₂ l` (up to escape overhead). **Therefore
`A_HUB(k) < 0` for every `k ≥ 3`, independent of `l`, independent of the data.**
`A_HUB(2) = 3 − 0 − 5 = −2` (marginal), `A_HUB(3) = 5 − 15.2 − 10 = −20.2`,
`A_HUB(17) = 9 − 228.2 − 80 = −299`, `A_HUB(143) = 15 − 2145.3 − 710 = −2 840`.

**Correction to the design notes, §3.** The note says a hub collapses `k²` edges to
`k` spokes and "the energy drop is the evidence the new object exists." The `k²`
figure is the unpruned detector output, not the codebook a proper MDL judge would
hold; the honest advantage is `(k−2) log₂ N + 5(k−1) − Γ(k)` bits over a spanning
tree, i.e. **linear in k, not quadratic**. The conclusion survives — hubs win from
`k = 3` — but for a different and cleaner reason, and the note should be amended
before probe 03 cites it. More importantly: this is a **pointer-economy argument,
not an existence argument**. Compression alone licenses hubs as bookkeeping. The
claim that a hub is a *real object* is carried entirely by the two properties that
are *not* tested here — that it predicts members outside the table, and that it pays
rent across a second table. Neither is a probe-02 fact and neither may be cited as
one.

---

## 6. Cells — full matrix, one batch, one round

All cells below are specified now and run in a single batch. No cell is added,
dropped, or re-parameterised after this commit.

### H — harness validity (any failure ⇒ VOID; the VOID is itself the published result)

- **H1 parse fidelity.** `serialise(parse(file)) == file` modulo declared
  whitespace normalisation, on 100% of rows of both files. Join integrity re-run
  from probe 01 H1.1 (every `allcurves` row joins exactly one `aplist` row).
- **H2 decoder test.** `decode(encode(T)) == T` **byte-identical**, for
  `T ∈ {T_cls, T_crv}` × `𝒦 ∈ {∅, K1-book, K2-book, K3-book, K4-book}`.
  Any single mismatched byte ⇒ VOID. This is the cell that makes the energy a
  description length rather than a score.
- **H3 accounting fidelity.** Actual bytes emitted by the range coder vs the
  analytic `Σ −log₂ q`: `|emitted − analytic| / analytic ≤ 0.005`. Failure means
  the reported energy is not the energy.
- **H4 baseline calibration.** `(L₀ − L₀_model)/(N·C_good)` equals the mean plug-in
  column entropy to within `0.01` bits/cell.

### N — nulls and anti-gaming (failure ⇒ terminal verdict, not a re-run)

- **N1 foil gate.** For each entry type in {BRIDGE, PATTERN, HUB}, generate
  `M = 1000` synthetic entries: member sets drawn uniformly at random with the size
  distribution matched to the real probe-01 group-size distribution, modulus drawn
  uniformly from `{5,7,11,13}`, reference row / residue vector drawn uniformly.
  Seed `default_rng(20260728)`. **Gate: `#{ΔE < 0} = 0` for each type, out of
  1000.** One compressing foil ⇒ verdict GAMEABLE.
- **N2 shuffled-table gate.** Marginal-preserving column shuffle (seed
  `20260728 + l`, identical procedure to probe 01 H2). The real probe-01 codebook
  entries, replayed on the shuffled table, must all have `ΔE > 0`.
  **Gate: `#{ΔE < 0} = 0`.**
- **N3 relabelling invariance.** `L₀` invariant to the bit under a random row
  permutation (seed `20260728`); `E(𝒦)` invariant under simultaneous relabelling of
  rows and entry pointers. **Gate: exact equality of the analytic total, and
  `≤ 1` byte difference in emitted stream length.**

### K — known answers (the scientific content; each licenses one move type for probe 03)

- **K1 DUP, substrate B.** Isogenous curves share `a_p` at every good prime, so
  `T_crv` contains exact duplicate rows by construction. Book: one DUP per
  non-representative curve in each class. **Gate: `ΔE < 0`, and the duplicated
  portion compresses to `≤ 5%` of its baseline cost.** If the energy cannot see an
  exact duplicate it cannot see anything; K1 failing ⇒ UNDERPOWERED, terminal.
- **K2 BRIDGE.** Book: MST of bridges over the 143-member mod-5 group from probe 01
  (`11a, 38b, 50b, …`). **Gate: `ΔE < 0`, and measured mean per-member saving
  within `±15%` of `−Σ_c log₂ Q_c(ρ) − 0.0227·C_l` evaluated on the measured
  baseline marginals.**
- **K3 PATTERN vs BRIDGE.** Same group, `π_EIS`, `l = 5`.
  **Gate: `A_PAT` measured within `±20%` of the P-b prediction `−2 838` bits, and
  `A_PAT < 0`.**
- **K4 HUB vs BRIDGE, distributional.** For **every** probe-01 group with `k ≥ 3`
  (all `l`), compute `A_HUB`. **Gate: `A_HUB < 0` in `≥ 99%` of groups, and the
  measured `A_HUB(k)` within `±20%` of `Γ(k) + (2−k) log₂ N − 5(k−1)` in `≥ 95%`
  of groups.** Also report, without gating, the `k = 2` groups (predicted
  marginal at `−2` bits) — the sign at `k = 2` is a diagnostic of coder overhead,
  not a verdict.
- **K5 TEST vs PATTERN, supplied predicate.** `φ = has-rational-l-isogeny` at
  `l = 5`, computed from the pinned ecdata isogeny matrices, exception list `X` =
  symmetric difference against the 143-member set. **Gate: `ΔE(TEST) < ΔE(PATTERN)`
  iff `|X| < k − (log₂|Φ| )/log₂ N`, i.e. iff the predicate is right about more
  members than the exception list costs.** Report `|X|` regardless — `|X| = 0` would
  mean the human definition and the trace signature coincide exactly on this range,
  which is the known answer probe 03's search must recover blind.
- **K6 sensitivity to `ε₀`.** Repeat K2–K5 at `ε₀ ∈ {2^-4, 2^-6, 2^-8}`.
  **Gate: no verdict flips.** A flip ⇒ SENSITIVE.
- **K7 headroom / ceiling (report-only, no gate).** Total `ΔE` of the full probe-01
  codebook (all 4 454 non-trivial groups, best entry type per group) as a fraction
  of `L₀`. This number exists so that a null in probe 03 can be distinguished from
  a ceiling: if the entire known taxonomy buys `< 1%` of `L₀`, then descent finding
  little is a headroom fact about the substrate, not evidence about the method.

---

## 7. Constants — v0, locked

| constant | value | why this value |
|---|---|---|
| `ε₀` | `2^-6` | escape costs 6 bits ≈ 3 residue classes' worth at `l = 5`; swept in K6 |
| `α` (KT smoothing) | `1/2` | standard; avoids zero-probability cells |
| `H` (entry header) | 5 bits | 3 type + 2 modulus, 5 declared types, 4 declared moduli |
| `ρ` (rent) | 0 bits | declared, inert in this probe; live from probe 04 |
| `log₂ N` | 15.215 (A) / 15.981 (B) | `N = 38042 / 64687` |
| `C` / `C_l` | 25 / 24 | `p = 2…97`; the `p = l` column excluded per probe 01 §4 |
| foil count `M` | 1000 per type | gives a `0/1000` gate with resolution `10^-3` |
| seed | `default_rng(20260728)` (+`l` where noted) | same as probe 01 |
| data pin | ecdata `25cec5ecfec8b9f016eb1631ac633194c2bed39f` | same commit as probe 01 |

---

## 8. Log schema

`output/energy.jsonl`, one record per `(cell, entry, ε₀)`:

```json
{"cell":"K3","eps0":-6,"entry_type":"PATTERN","l":5,"k":143,
 "members_hash":"sha256:…","cost_entry_bits":2192.7,
 "saving_bits":7894.1,"delta_E_bits":-5701.4,
 "predicted_delta_E_bits":-5698.0,"rel_err":0.0006,
 "escapes":0,"cells_claimed":3432,"gate":"PASS"}
```

Plus `output/baseline.json` (`L₀`, `L₀_model`, per-column entropies, alphabet
sizes), `output/foils.csv` (one row per foil: type, `l`, `k`, `ΔE`),
`output/decoder.json` (per `(T, 𝒦)`: byte-identical true/false, emitted bytes,
analytic bits, ratio), `output/headroom.json` (K7).

Every record carries `members_hash` so that probe 03's registry can auto-reject a
re-proposal of an entry already scored here.

---

## 9. Interpretation ladder — applied mechanically, in this order

1. **VOID** — any of H1–H4 fails. Published as-is; the harness failure is the
   result and the probe is closed.
2. **GAMEABLE** — N1 fails (any foil compresses). Terminal. The energy is
   withdrawn. A redesigned energy is a **new prereg**, not a re-run of this one.
3. **NULL-LEAK** — N2 fails. Terminal, same treatment as 2.
4. **INVARIANCE-FAIL** — N3 fails. Terminal, same treatment as 2.
5. **UNDERPOWERED** — H and N pass, K1 fails. The energy cannot see exact
   duplicates. Terminal; probe 03 does not start.
6. **SENSITIVE** — H, N, K1–K5 pass but K6 flips a verdict. The energy is
   tuning-dependent; recorded as such, and probe 03 inherits only the move types
   that were stable across all three `ε₀`.
7. **PARTIAL** — H, N, K1 pass; some of K2–K5 fail. Terminal. Each failing move
   type is **deleted from probe 03's vocabulary** and logged as a kill. This is a
   completed result, not a reason to repair the entry type.
8. **VALIDATED** — H, N pass; K1–K5 pass; K6 shows no flips.

"Inconclusive at this scale with this method" is a binding terminal verdict,
recorded in `REGISTRY.md` as CLOSED, never as pending. No parameter changes after
this commit. No follow-up round is authorised.

---

## 10. Pre-registered falsifiers, and the load-bearing risk

The prediction most likely to break, and the one worth breaking:

- **P-c is a bookkeeping theorem, not an existence theorem.** If K4 passes exactly
  as predicted — which it should, since the cancellation is analytic — that
  establishes only that hubs are cheaper bookkeeping than spanning trees. It does
  **not** establish that a hub names a real object. The design notes' claim that
  "the energy drop is the evidence the new object exists" is not supported by
  anything in probe 02 and must not be cited from it. The evidence for existence is
  out-of-table prediction and cross-table rent, both of which are probe 04+.
- **If K4 fails**, the escape-the-table move is not licensed by compression at all,
  and §3 of the design notes needs rewriting before probe 03 is specified.
- **If N1 fails at any type**, the whole generator design is suspended: a gameable
  energy under a Gibbs process manufactures structure at exactly the rate the
  temperature allows, and every downstream discovery would be an artifact.
- **If K7 shows the entire known taxonomy buys `< 1%` of `L₀`**, then the substrate
  is nearly incompressible by congruence structure at this conductor range, and
  probe 03's descent is being run near its ceiling. That is a headroom fact to
  publish alongside probe 03's verdict, not a reason to change probe 03.

## 11. What a VALIDATED verdict does NOT establish

No new mathematics, no conjectures, no discovery. VALIDATED here means exactly one
thing: **the energy is a real code length, is not gameable by the declared foils,
and pays the predicted number of bits for the four moves whose answers we already
know.** It says nothing about whether descent over this energy finds anything,
nothing about whether the move vocabulary is complete (that is probe 03's
completeness audit), and nothing about whether any entry is true beyond the encoded
range (that is the sealed range, §12).

## 12. Seal

Conductors `≥ 10000` are sealed in the same commit as this document
(`SEAL.md`). Nothing in probe 02 reads them. Every codebook entry minted on the
training range gets exactly one confrontation with the sealed range, ever, enforced
by the registry.
