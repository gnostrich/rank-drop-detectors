# rank-drop-detectors

Rank-drop detectors for arithmetic objects. Generic gives zero; exceptional
gives a number; the number is the finding.

## The idea

Arithmetic objects — elliptic curves, modular forms, Galois representations —
each carry local data: an integer `a_p` for every prime `p`. The sequence of
those integers is a fingerprint.

The programme is a single move, applied in different settings: establish what
generic looks like, measure, and where the measurement departs from generic,
something is there. The size of the departure is graded, not binary.

## Two validated detectors

**Detector A — instance level.** Combine two objects; a degenerate
combination shows up in its own numbers with no lookup required. For the
tensor product of two elliptic curves, the degree-4 Euler polynomial is

```
1 - a1*a2*T + (p*(a1^2 + a2^2) - 2*p^2)*T^2 - p^2*a1*a2*T^3 + p^4*T^4
```

and substituting `T = 1/p` gives `(a1 - a2)^2 / p`. The polynomial degenerates
at `p` exactly when `a1 = a2`. Measured on Cremona data (fraction of primes
where `a1 = a2`): isogenous **1.000**, quadratic twists **0.528** (theory
0.5), unrelated **0.068**; 5 of 5000 unrelated pairs exceeded 0.35.

**Detector B — structure level.** Instead of comparing two objects, compute
the dimension of the space of maps between two structures. `dim Hom(A, B)` is
the nullity of a linear system — solved, not searched. Measured: unrelated
matrices **0**, conjugate **n**, sharing one eigenvalue **1–2**. On nilpotent
Jordan types (all with identical spectrum, so spectrum comparison is blind) it
separates all 11 partitions of 6, from 6 to 36. On quiver representations the
self-Hom recovers representation type exactly: finite **4**, tame **2**,
wild **1**.

**The finding that dictates the design: random families are useless.** Across
800 random pairs in each of three quiver families, `dim Hom` was exactly 0,
every time, no exceptions. Genericity is precisely the absence of structure.
The family must have arithmetic origin.

## Current probe: congruences between elliptic curves

Two elliptic curves have a nonzero Hom between their mod-`l` Galois
representations exactly when

```
a_p(E1) = a_p(E2)   (mod l)   for all p
```

Congruences of this kind are the mechanism behind level-lowering and Serre's
conjecture, so a hit is not a curiosity; it is the species of fact that
carries weight. The family is arithmetic, not random; the null is effectively
exact (for `l >= 5`, an accidental congruence across 22 primes has probability
about `l^-22`); and the harness has a known answer: isogenous curves are
congruent modulo everything, so the detector must catch them first.

See [`probes/01-congruence-census/`](probes/01-congruence-census/):
[`PREREG.md`](probes/01-congruence-census/PREREG.md) (committed before any
run), [`run.py`](probes/01-congruence-census/run.py), and
[`RESULTS.md`](probes/01-congruence-census/RESULTS.md) (one page: numbers,
ladder verdict, what is not claimed). Closed probes are logged one line each
in [`REGISTRY.md`](REGISTRY.md).

## Reproducing

```
git clone --filter=blob:none --no-checkout https://github.com/JohnCremona/ecdata.git
cd ecdata && git sparse-checkout set aplist allcurves && git checkout 25cec5ecfec8b9f016eb1631ac633194c2bed39f
cd .. && python3 probes/01-congruence-census/run.py --ecdata ecdata --out probes/01-congruence-census/output
```

Requires Python 3 and numpy. MIT licence.
