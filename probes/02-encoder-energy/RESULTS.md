# RESULTS — Probe 02: the encoder and its energy — **PARTIAL**

Run 2026-07-28, one round, terminal. Data: ecdata `25cec5ec`, conductors
1–9999; nothing in the sealed range was read (path guard enforced). Prereg
and seal committed before any code existed. Implementation readings R1–R5
(forced by the prereg text, e.g. marker cells take the escape branch — the
only decoder-realizable form) are documented in `run.py`'s docstring. Full
artifacts in [`output/`](output/).

## Usable data (reported before any energy number)

38,042 classes / 64,687 curves, 100% usable; 25 columns parsed, 24 used per
modulus. Baselines: **L₀(A) = 3,974,375 bits**, L₀(B) = 6,396,373 bits.
Probe-01 groups reconstructed exactly: 3521/896/32/5, the 143-family intact
with 11a. K5 predicate: 534 classes have a rational 5-isogeny; |X| = 391.

## H — harness: PASS

H1: 0 serialise∘parse failures on both files; join clean. H2: **7/7**
well-typed (table × book) streams decode **byte-identical**. H3: emitted
vs analytic ≤ 1.5×10⁻⁵ on every stream (bar 5×10⁻³) — and the first launch
was caught by H3 exactly as designed: a ledger bug (DUP delta missing its
baseline term) showed as a 37% gap, was fixed, rerun. H4: exact.

## N — nulls: PASS

N1: **0/1000 foils compress, per entry type**. N2: **0/31,243** real
entries compress when replayed on marginal-preserving shuffled tables.
N3: L₀ exactly invariant; book energies invariant to <10⁻⁹ bits; permuted
stream lengths differ by 0 bytes.

## K — known answers

| cell | predicted | measured | gate | verdict |
|---|---|---|---|---|
| K1 DUP | ΔE < 0, ≤5% of baseline | ΔE = −1,380,808 b; **0.64%**; 0 escapes | both | **PASS** |
| K2 BRIDGE | mean saving 54.9 b; ΔE < 0 | saving 47.4 b (−13.7%, in band); **ΔE = +52.96 b** | ΔE < 0 | **FAIL** |
| K3 PATTERN | A_PAT = −2,838 b (±20%) | **−3,185 b** (12.2% off); ΔE = −2,928 b | both | **PASS** |
| K4 HUB | <0 in ≥99%; in ±20% band ≥95% | **negative in 714/714**; in band **35.2%** | band | **FAIL** |
| K5 TEST | TEST beats PATTERN iff \|X\| < k−0.25 | 391 ≥ 142.75 and TEST loses (+851 vs −2,928) | iff holds | **PASS** |
| K6 | no verdict flips across ε₀ | K2 flips at ε₀ = 2⁻⁴ (recorded; K2 already failed) | — | noted |
| K7 | report-only | full taxonomy (5,275 entries) buys 66,523 b = **1.67% of L₀** | — | reported |

What the two failures are, mechanically. **K2:** a mod-5 bridge saves
~47 bits at real (non-equidistributed) marginals but pays 35.4 pointer bits
plus ~12 bits of marker escapes — net **+0.4 bits per bridge**. Pairwise
bridges sit at break-even and land positive; the compression lives in
family-level entries. **K4:** the sign claim of P-c is *universal* on this
data (hub beats spanning tree in every one of 714 groups) but the ±20%
quantitative band fails in 65% of groups: the S_root = C_l·log₂ l
cancellation assumes marker-free roots and equidistributed residues, and
real data has neither. k=2 diagnostic as predicted marginal: mean A = +7.3,
1.6% negative — coder overhead dominates at k=2.

## Ladder verdict (mechanical)

H pass, N pass, K1 pass, K2 and K4 fail → **PARTIAL** (§9.7), terminal.
Probe 03's move vocabulary inherits **DUP, PATTERN, TEST**. **BRIDGE and
HUB are logged as kills** — HUB's kill is the quantitative prediction band,
not the sign (714/714 negative), and per §9.7 it is deleted, not repaired.
K7's 1.67% headroom is the ceiling context probe 03 must be read against.

## What is NOT claimed

The energy is a real code length (it decodes), is not gameable by the
declared foils, leaks nothing under the null, and pays the predicted bits
for DUP/PATTERN/TEST — nothing more. No claim that hubs name real objects
(P-c was bookkeeping; even its band failed on real marginals). No claim
about BRIDGE at other moduli than tested. No discovery, no new mathematics,
no statement about the sealed range, which remains sealed and unread.
