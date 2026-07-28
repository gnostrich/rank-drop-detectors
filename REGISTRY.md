# REGISTRY — closed probes, append-only, one line each

| probe | closed | substrate | verdict | numbers |
|-------|--------|-----------|---------|---------|
| 01-congruence-census | 2026-07-28 | mod-l congruences, elliptic curves N<10000, l∈{5,7,11,13} | VALIDATED | 38,042 classes; H1 pass; H2 = 0/0/0/0; A1 groups 3521/896/32/5; 22,593 persistent pairs, 100% pattern-consistent, P4 = 0; 5 finite-depth coincidences |
| 02-encoder-energy | 2026-07-28 | MDL energy over the probe-01 corpus; H/N/K matrix, codebooks hand-supplied | PARTIAL | H: 7/7 byte-identical, ≤1.5e-5; N: 0/3000 foils, 0/31,243 null replays, invariant; K1 −1.38 Mb (0.64%) PASS; K2 ΔE +53 b FAIL → BRIDGE killed; K3 −3,185 vs −2,838 PASS; K4 sign 714/714, band 35% FAIL → HUB killed; K5 \|X\|=391 PASS; K6 K2 flips at 2^-4; K7 headroom 1.67% of L₀ |
| 03-transport | 2026-07-28 | LMFDB cross-family (EC/Q, wt-2 newforms, EC+HMF over 2.2.5.1/2.2.8.1); blind modularity, transport cycles, closure defect | LOSSY (CLOSED) | H pass; N: 0 shuffled matches, foil power 1.00, floors 0; K1 blind modularity precision 1.0000 recall 1.0000 (2,137/2,137, 0 residuals); K2/C1+C2: 2,241 cycles all δ=0; K5: 79/2,343 above floor, all C3, δ_frac=0, conc 0.41 → transport lossiness at bad-reduction columns; K3: curves ~0.5–1.8%, forms 46–50%; K4: (a) flat, (b) growing → compounding scales with families |
