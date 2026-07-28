# SEAL — held-out range, declared 2026-07-28

Committed in the same commit as `probes/02-encoder-energy/PREREG.md`, before any
encoder run.

## What is sealed

All elliptic curves of **conductor ≥ 10000** in ecdata at the pinned commit
`25cec5ecfec8b9f016eb1631ac633194c2bed39f` (the database extends to 500000).

Concretely: every `aplist/aplist.NNNNN-NNNNN` and `allcurves/allcurves.NNNNN-NNNNN`
file whose range starts at or above 10000.

## The rule

1. **No probe reads the sealed range except through a registered confrontation.**
2. A confrontation is registered *before* it happens: a line in `REGISTRY.md`
   naming the codebook version, the entry (by `members_hash`), the prediction, and
   the pass/fail criterion.
3. **One look per hypothesis, ever.** A codebook entry that has been confronted is
   marked spent. A revised entry is a new entry with a new hash and gets its own
   single look; a re-look at the same hash is auto-rejected.
4. A failed confrontation is a published result, not a reason to revise and retry
   on the same range. Revisions are confronted on the *next* unsealed band
   (conductor ≥ 100000), which is sealed under the same rule.
5. Nothing that touches the sealed range may feed back into the training-range
   encoder, the marginals `q_c`, or the move vocabulary. The pipe is one-way.

## Why it is declared now

The encoder's marginals and the move vocabulary are both fitted on the training
range. If the sealed range were consulted while either was being designed, the
held-out verdict would be a training number wearing a held-out label, and every
downstream discovery claim would inherit the contamination. Sealing before the
encoder exists is the only version of this that is worth anything.

## Bands

| band | conductors | status | released by |
|---|---|---|---|
| train | 1 – 9999 | open | — |
| seal-1 | 10000 – 99999 | **sealed** | registered confrontation, one per hypothesis |
| seal-2 | 100000 – 500000 | **sealed** | reserved for revised hypotheses only |

## Probe-03 bands (added 2026-07-28, same commit as probes/03-transport/PREREG.md)

Top decile of each newly added family's declared band, sealed before any scan.
Sealed rows are downloaded and stored but excluded from every scan, codebook,
transport, and measurement; one-look rule as above.

| family block | measure | open | sealed |
|---|---|---|---|
| F2 newforms wt 2 (Q) | level | 1 – 900 | **901 – 1000** |
| F3a EC / 2.2.5.1 | conductor norm | 1 – 900 | **901 – 1000** |
| F3b EC / 2.2.8.1 | conductor norm | 1 – 900 | **901 – 1000** |
| F4a HMF / 2.2.5.1 | level norm | 1 – 900 | **901 – 1000** |
| F4b HMF / 2.2.8.1 | level norm | 1 – 900 | **901 – 1000** |
