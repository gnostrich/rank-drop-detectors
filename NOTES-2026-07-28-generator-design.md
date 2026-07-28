# Session notes, 2026-07-28 — from probe 01 to the generator design

Written for a cold reader (human or another session). Context: probe 01 shipped
this session (see `probes/01-congruence-census/`, verdict VALIDATED, merged as
PR #1). The rest of the session designed what gets built next. These notes are
design, not results — nothing below has been run.

## 1. The machine, split in two

**Part 1 — the generator.** Self-contained; never touches Lean. A Gibbs
process running on fingerprint tables. Data in, discoveries out. All the
scientific content lives here.

**Part 2 — the notary.** A small, stateless, downstream exporter that turns
finished discoveries into Lean artifacts. Adds trust (portability,
kernel-checked certificates), not truth. Can lag years behind Part 1 without
stalling it. The pipe is one-way: nothing flows from Lean back into the
generator (see §6 for why this is now a hard rule).

## 2. Part 1 concretely

- **Substrate:** fingerprint tables — one row per arithmetic object (elliptic
  curve, modular form, …), one column per prime, entries a_p. Currently:
  38,042 isogeny classes, conductor < 10000 (probe 01's corpus). Fixed,
  append-only, raw invariants only.
- **State ("the library" / codebook):** the current list of entries the
  encoder may use — bridges ("row X matches row Y mod 5"), rules ("these k
  rows obey a_p ≡ p+1 mod 5"), membership tests, hubs (see §3), plus a
  baseline value distribution. On day one this is a JSON file, nearly empty.
- **Energy:** E = bits(codebook) + bits(table | codebook) + rent per unproved
  assumption. Real code length — the encoding must actually decode back to
  the table (see §5, decoder test).
- **Process:** propose a small codebook edit; accept if E drops (with
  temperature for occasional uphill moves); cool. Detectors (probe-01-style
  scans) are the proposal distribution — legs, not judge. Falsification is
  not a module: it is the energy re-evaluated on data the entry hasn't seen.
- **Held-out seal:** conductors ≥ 10000 (ecdata has them to 500000) are
  sealed NOW. Every codebook entry minted on the training range gets exactly
  one confrontation with the sealed range. One look per hypothesis, enforced
  by the registry.

## 3. How discovery output escalates (bridge → rule → test → hub)

- A **bridge** is one match. A **rule** bundles many. Both still reference
  rows.
- A **membership test** replaces a rule's explicit member list with a
  computable predicate — that step is definition-formation in the MDL sense.
  Known-answer check: the 143-member mod-5 family from probe 01 must yield
  the predicate "has a rational 5-isogeny" (the human definition, recovered
  blind).
- A **hub** is the move that escapes the table: posit ONE new node that all
  k members of a cluster map to (probe 01's 143-family: the mod-5 Eisenstein
  representation — not an elliptic curve, not a row in the table).
  *Amended (probe 02, prereg §5.1 P-c):* the honest comparator is not the
  ~k² unpruned detector edges — an MDL judge already holds a spanning tree
  of k−1 bridges. Against that tree, the hub's own residue vector costs
  C·log₂ l and the tree's root row forgoes exactly S_root = C·log₂ l; these
  cancel identically, and what survives is pointer economy alone:
  `A_HUB(k) = Γ(k) + (2−k)·log₂ N − 5(k−1)` — linear in k, not quadratic;
  independent of l and of the data; negative for every k ≥ 3. The
  conclusion (hubs win from k = 3) survives; the mechanism differs.
  [Measured in probe 02: the sign held in 714/714 groups; the ±20%
  quantitative band held in only 35% — the cancellation is a clean-marginals
  idealization. See `probes/02-encoder-energy/RESULTS.md`.]
  *Amended (probe 02, §10):* the energy drop is evidence that hubs are
  cheaper **bookkeeping**, not evidence that the new object exists.
  Existence is carried entirely by the two properties compression does not
  test — out-of-table prediction and cross-table rent, below — and neither
  is a probe-02 fact; neither may be cited as one. A hub earns "real
  definition, not a column" only through those untested properties: it
  predicts members outside the table, it gets its own fingerprint row
  (becomes a first-class object), and it must pay rent across tables (a hub
  that only ever organizes one table was a coincidence dressed up).
- **Concept-gap case:** a cluster that compresses but matches no predicate in
  the declared vocabulary = a proven missing concept, with its extension
  enumerated. That is a deliverable, not a failure.

## 4. Why this can produce new mathematics (the functoriality story)

A fingerprint match IS a map: the Langlands programme predicts that matching
a_p sequences are always explained by a functorial relationship. Modularity —
the theorem behind FLT — began as exactly this species of noticed fingerprint
coincidence. The pipeline:

1. (done) Table = elliptic curves only; every bridge is known math → harness.
2. (next harness) Add modular forms to the table; the machine must rediscover
   modularity as the giant bridge family.
3. (live ammo) Add families where the map is genuinely incomplete: high
   conductors, genus-2 curves / abelian surfaces (paramodular conjecture
   territory).
4. Output = residuals: "X and Y match at hundreds of primes, survive the null
   and the sealed range, and no codebook entry explains it." That sentence is
   a conjecture of functoriality type.

Honest odds (stated in-session and held): modal outcome is a validated
instrument plus rediscoveries. The murmuration-class find (a real unnoticed
regularity) is a low-single-digit-percent event. The five breaking pairs from
probe 01 ({5408a, 7200bg, 7200br, 9248e}, mod 5, broke at p = 13/101/113) are
the first residual specimens on file.

## 5. Architecture with faithfulness/naturalness checks, piece by piece

Governing principle: every arbitrary-looking choice either passes an
invariance test (result independent of it) or is prereg-locked with a
sensitivity report.

1. **Substrate.** Faithfulness: independent a_p recomputation by point
   counting on random samples, zero mismatch tolerance (probe 01 H1.3, kept
   forever); join integrity. Naturalness: data is external and pinned by
   commit — we cannot have sculpted it.
2. **Energy.** Faithfulness: the **decoder test** — encode, decode, get the
   table back byte-identical; otherwise it is a score, not a description
   length, and scores are gameable. Naturalness: empty-codebook cost must
   match empirical entropy; energy invariant under row permutation; the
   **anti-gaming foil gate** (synthetic conjunction-style entries must reduce
   energy 0% of the time — transplanted from mathlib-structure-probes);
   encoder fixed before any run.
3. **Proposal distribution (detectors).** Faithfulness: completeness audit —
   on a small slice, exhaustive search and detector-guided proposals must
   accept the identical move set; declared blind spots only. Naturalness:
   **path-independence** — different seeds/orderings must converge to the
   same (or energy-tied) codebook.
4. **Moves.** Closed, pre-declared vocabulary; each move type licensed by a
   known answer: isogeny bridges compress near-perfectly; the 143-rule beats
   its member list; test-search recovers "5-isogeny"; hub-factoring collapses
   the 143-hairball with predicted-vs-measured energy drop. No move types
   invented mid-run — a needed-but-missing move is recorded as a gap and
   added in the NEXT prereg.
5. **Schedule.** Faithfulness: deterministic given seed; plateau criterion
   preregistered. Naturalness: final codebook invariant across the declared
   cooling-rate range.
6. **Machine-level nulls.** (a) Global null: entire machine on
   marginal-preserving shuffled tables → codebook must converge empty.
   (b) **Implant test:** plant a synthetic congruence family in shuffled
   data → machine must recover exactly the implant, nothing else (measures
   detection power instead of assuming it). (c) Held-out seal as in §2.
7. **Output format.** Every discovery carries: statement, bits saved, null
   rate, held-out verdict, and the codebook version it was not explained by —
   with the codebook's coverage of textbook phenomena published as a
   checklist (isogeny, twist, Eisenstein, CM, …). Registry of kills is
   append-only; re-proposals of registered corpses are auto-rejected.

Build order: **probe 02** = pieces 1–2 (substrate + energy; decoder test and
foil gate as H-cells; known-answer compression gates). **Probe 03** = pieces
3–5 (moves + descent; blind rediscovery of the known taxonomy;
path-independence). Piece 6 wraps both; piece 7 is standing discipline.

## 6. The mathlib failure ledger — binding constraints

Read this session from `gnostrich/mathlib-structure-probes` (see its README
and `harness/STATUS.md`): eleven probes on Mathlib's structure run to
terminal verdicts; eight killed outright, one "do not wire", two partial
salvages. Standing positives there: human importance ≈ foundational depth
(AUC 0.81) while reuse ≈ chance (0.52); de-gamed MDL compression adds nothing
beyond depth (0.569 < locked 0.58 line); AUC + permutation, never rank
correlation, under extreme class imbalance; a working deterministic
Lean-variant harness.

Constraints this imposes here (things we almost repeated in-session and must
not):

- **No MDL/compression judge pointed at the Mathlib graph.** Refuted there at
  a pre-committed line. The compression judge lives on arithmetic data only.
- **No Gibbs moves in formal-statement space.** ~95% of local perturbations
  to real Mathlib theorems break immediately (probe 16 ceiling): statement
  space is a cliff field, not a descent landscape.
- **No mining Mathlib structure for guidance/taste/interestingness.** Seven
  attempts, all reduced to degree/gradient/blind.
- Lean is notary and checker only. Also practical: Lean toolchain downloads
  were egress-blocked in these sessions — the notary needs a remote
  built-Mathlib runner.
- Meta-lesson (learned live this session): elegant dead ideas get re-derived
  by fresh contexts. The kill registry is load-bearing — proposals matching
  registered kills are rejected without re-evaluation, and any new prereg
  cites the kills it builds past.

## 7. Part 2 (notary) in four levels

1. **Statement:** the discovery as a Lean proposition (a_p via point counts is
   definable in mathlib today). Round-trip check: Lean recomputes sample a_p
   values against the raw data before the statement is trusted.
2. **Finite certificate:** kernel-checked verification up to the tested bound.
3. **Conditional reduction:** Sturm bound + modularity turn the finite check
   into the infinite statement — Lean certifies "IF those two, THEN the
   discovery." Buildable now; this is the standard export form.
4. **Unconditional:** blocked on the FLT-formalization effort landing
   modularity/Sturm in mathlib; our conditional entries flip automatically
   when it does. (Some Eisenstein-family facts can go unconditional early via
   explicit torsion points.)

## 8. Immediate next action

Write `probes/02-.../PREREG.md`: encoder spec, decoder test, foil gate,
known-answer compression gates, seeds, thresholds, ladder — committed before
any encoder runs; seal conductors ≥ 10000 in the same commit. One round,
terminal, as always.
