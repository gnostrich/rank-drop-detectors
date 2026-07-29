#!/usr/bin/env python3
"""LOOP v0 — the promotion cycle, S-A only. Implements loop/PREREG.md.

Usage:
    python3 loop.py --ecdata /path/to/ecdata --out output

Confronter: specified in the design document, NOT run here. No sealed row
is read by any cell in this file. SA-only cluster scans are cached per
dataset and reused across rounds and orderings (rounds >= 2 add only
mint-involved components, computed incrementally); this is caching of a
deterministic function, not a parameter.
"""

import argparse
import hashlib
import importlib.util
import json
import math
import os
import time
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
          67, 71, 73, 79, 83, 89, 97]
C = 25
MODULI = [5, 7, 11, 13]
SEED = 20260728
FOLDS = 5
N_HOLD = 10
THETA = 1e-6
R_MAX = 5
ORDERINGS = 20
MIN_SHARED = 10


def load_module(name, relpath):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, relpath))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


p05 = load_module("p05r", "../probes/05-mint/run.py")
p03, p02, p01 = p05.p03, p05.p02, p05.p01


def hkey(l, members):
    return hashlib.sha256(
        (str(l) + str(sorted(members))).encode()).hexdigest()[:16]


def fold_masks():
    out = []
    for f in range(FOLDS):
        rng = np.random.default_rng(SEED + f)
        out.append(sorted(int(x)
                          for x in rng.choice(C, N_HOLD, replace=False)))
    return out


MASKS = fold_masks()


def sa_cluster_cache(SA):
    """(fold, l) -> list of SA-only components (row-index lists, k>=2)."""
    cache = {}
    V, G = SA.val_good()
    for fold in range(FOLDS):
        fit = [c for c in range(C) if c not in MASKS[fold]]
        for l in MODULI:
            cols = [c for c in fit if PRIMES[c] != l]
            vm = (V[:, cols] % l).astype(np.int8)
            comps, _ = p01.scan(vm, G[:, cols], tag=f"c{fold}-{l}",
                                log=lambda *a, **k: None)
            cache[(fold, l)] = [sorted(m) for m in comps]
    return cache


def mint_components(SA, mint_rows, fold, l):
    """Components involving >=1 mint row, on fit columns (incremental)."""
    if not mint_rows:
        return []
    fit = [c for c in range(C) if c not in MASKS[fold]]
    cols = [c for c in fit if PRIMES[c] != l]
    V, G = SA.val_good()
    Vs, Gs = V[:, cols] % l, G[:, cols]
    out = []
    for mi, row in enumerate(mint_rows):
        mv = np.array([0 if row[c] is None else row[c] % l for c in cols])
        mg = np.array([row[c] is not None for c in cols])
        both = Gs & mg[None, :]
        neq = (Vs != mv[None, :]) & both
        ok = (~neq.any(axis=1)) & (both.sum(axis=1) >= MIN_SHARED)
        sa_members = [("SA", int(i)) for i in np.nonzero(ok)[0]]
        if sa_members:
            out.append([("MINT", mi)] + sa_members)
    return out


def fit_residues(fams_by_name, members, l, fit):
    rh = {}
    for c in fit:
        if PRIMES[c] == l:
            continue
        for (fn, i) in members:
            v = fams_by_name[fn].rows[i][c]
            if v is not None:
                rh[c] = v % l
                break
    return rh


def anchor_type(rh, l):
    for pi in p02.PI_ORDER:
        fn_ = p02.PI_FUN[pi]
        if rh and all(v == fn_(PRIMES[c]) % l for c, v in rh.items()):
            return "closed:" + pi
    return "free"


def score(fams_by_name, members, l, atype, held):
    ref = min(members)
    hits = n = 0
    for c in held:
        if PRIMES[c] == l:
            continue
        if atype.startswith("closed:"):
            pred = p02.PI_FUN[atype.split(":")[1]](PRIMES[c]) % l
            targets = members
        else:
            rv = fams_by_name[ref[0]].rows[ref[1]][c]
            if rv is None:
                continue
            pred = rv % l
            targets = [m for m in members if m != ref]
        for (fn, i) in targets:
            v = fams_by_name[fn].rows[i][c]
            if v is None:
                continue
            n += 1
            if v % l == pred:
                hits += 1
    return hits, n


def binom_tail(n, k, p):
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    lp, l1p = math.log(p), math.log1p(-p)
    total = 0.0
    for x in range(k, n + 1):
        lg = (math.lgamma(n + 1) - math.lgamma(x + 1)
              - math.lgamma(n - x + 1) + x * lp + (n - x) * l1p)
        total += math.exp(lg)
        if lg < -700:
            break
    return min(1.0, total)


def run_loop(SA, cache, out_records, order_seed=None, tag="real"):
    fams_by_name = {"SA": SA}
    quarantine = set()
    promoted = {}
    mint_rows, mint_labels, mint_depth = [], [], []
    rounds = []
    log2N0 = math.log2(SA.N)
    owned = set()
    cum_dE = 0.0
    bits_traj = []
    for R in range(1, R_MAX + 1):
        scores = []
        for fold in range(FOLDS):
            fit = [c for c in range(C) if c not in MASKS[fold]]
            held = MASKS[fold]
            for l in MODULI:
                comps = [[("SA", i) for i in m] for m in cache[(fold, l)]]
                comps += mint_components(SA, mint_rows, fold, l)
                for members in comps:
                    if len(members) < 2:
                        continue
                    key = hkey(l, members)
                    if key in quarantine or key in promoted:
                        continue
                    rh = fit_residues(fams_by_name, members, l, fit)
                    atype = anchor_type(rh, l)
                    hits, n = score(fams_by_name, members, l, atype, held)
                    if n == 0:
                        continue
                    scores.append({"round": R, "fold": fold, "key": key,
                                   "l": l, "k": len(members), "type": atype,
                                   "n": n, "hits": hits,
                                   "p_tail": binom_tail(n, hits, 1.0 / l),
                                   "members": members})
        m = max(1, len(scores))
        for s in scores:
            s["promoted"] = s["p_tail"] < THETA / m
        order = list(range(len(scores)))
        if order_seed is not None:
            order = [int(x) for x in
                     np.random.default_rng(order_seed).permutation(len(scores))]
        promoted_this = []
        seen = set()
        for si in order:
            s = scores[si]
            if s["promoted"] and s["key"] not in seen:
                seen.add(s["key"])
                promoted_this.append(s)
        passed = {s["key"] for s in scores if s["promoted"]}
        for s in scores:
            if s["key"] not in passed:
                quarantine.add(s["key"])
        for s in promoted_this:
            members = s["members"]
            depth = (1 + max([mint_depth[i] for (fn, i) in members
                              if fn == "MINT"], default=0)
                     if any(fn == "MINT" for (fn, _i) in members) else 0)
            ref = min(members)
            row = []
            for c in range(C):
                v = (None if PRIMES[c] == s["l"]
                     else fams_by_name[ref[0]].rows[ref[1]][c])
                row.append(None if v is None else v % s["l"])
            # full-column dE at promotion time (persistent ownership)
            eng = p05.AnchorEngine(list(fams_by_name.values()), log2N0)
            rh_full = np.array([-1 if row[c] is None else row[c]
                                for c in range(C)], dtype=np.int16)
            dE = eng.anchor_delta(members, s["l"], rh_full, owned, claim=True)
            cum_dE += dE
            mint_rows.append(row)
            mint_labels.append(s["key"])
            mint_depth.append(depth)
            promoted[s["key"]] = {**{k: v for k, v in s.items()
                                     if k != "members"},
                                  "depth": depth, "dE_full": dE}
        for s in scores:
            out_records.append({**{k: v for k, v in s.items()
                                   if k != "members"}, "tag": tag})
        if mint_rows:
            MINT = p05.Fam("MINT", [p05.MINT_NORM] * C, mint_rows,
                           mint_labels)
            p05.uniformize(MINT.table)
            fams_by_name = {"SA": SA, "MINT": MINT}
        bits = SA.table.L0 + cum_dE
        bits_traj.append(bits)
        rounds.append({"round": R, "scored": len(scores),
                       "promoted": len(promoted_this),
                       "promoted_total": len(promoted),
                       "quarantined_total": len(quarantine),
                       "bits_state": round(bits, 1)})
        if not promoted_this:
            break
    return {"rounds": rounds, "promoted": promoted,
            "chain_depths": mint_depth, "bits_traj": bits_traj,
            "cum_dE": cum_dE, "quarantine": quarantine,
            "final_set": set(promoted.keys())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ecdata", required=True)
    ap.add_argument("--out", default="output")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    report = {"seed": SEED, "theta": THETA, "folds": FOLDS, "R_max": R_MAX}

    apl, _ = p02.parse_aplist(os.path.join(args.ecdata, "aplist",
                                           "aplist.00000-09999"))
    sa_rows = [[None if s in ("+", "-", "?") else int(s) for s in syms]
               for (_N, _c, syms, _t) in apl]
    sa_labels = [f"{N}{c}" for (N, c, _s, _t) in apl]
    SA = p05.Fam("SA", PRIMES, sa_rows, sa_labels)

    # ---------- H ----------
    print("== H ==")
    fit0 = [c for c in range(C) if c not in MASKS[0]]
    r1 = fit_residues({"SA": SA}, [("SA", 0), ("SA", 1), ("SA", 2)], 5, fit0)
    mut = [list(r) for r in sa_rows[:3]]
    for c in MASKS[0]:
        for i in range(3):
            if mut[i][c] is not None:
                mut[i][c] = max(min(mut[i][c] + 1, 18), -18)
    SAm = p05.Fam("SAm", PRIMES, mut + sa_rows[3:10], sa_labels[:10])
    r2 = fit_residues({"SA": SAm}, [("SA", 0), ("SA", 1), ("SA", 2)], 5, fit0)
    ha = r1 == r2
    with p03.p02_norms(PRIMES):
        codec = p02.StreamCodec(SA.table, [], 6)
        hb = bool((codec.decode(codec.encode()) == SA.table.SYM).all())
    hc = True
    rngc = np.random.default_rng(SEED)
    for l in MODULI:
        n = 100
        sims = rngc.binomial(n, 1.0 / l, size=1_000_000)
        for target in (1e-3, 1e-4):
            k = 1
            while binom_tail(n, k, 1.0 / l) > target:
                k += 1
            emp = float((sims >= k).mean())
            ana = binom_tail(n, k, 1.0 / l)
            ratio = (emp / ana) if ana > 0 and emp > 0 else 1.0
            if not (0.5 <= ratio <= 2.0):
                hc = False
    report["H"] = {"masked_fit": ha, "coder": hb, "scorer_cal": hc,
                   "pass": ha and hb and hc}
    print(report["H"], f"({time.time()-t0:.0f}s)")

    # ---------- real loop ----------
    print("== scans (cached) ==")
    cache = sa_cluster_cache(SA)
    print(f"cache built ({time.time()-t0:.0f}s)")
    records = []
    real = run_loop(SA, cache, records, tag="real")
    report["real_rounds"] = real["rounds"]
    report["bits_traj"] = [round(b, 1) for b in real["bits_traj"]]
    report["chain_depth_hist"] = {str(d): real["chain_depths"].count(d)
                                  for d in set(real["chain_depths"])} or {"0": 0}
    promoted = list(real["promoted"].values())
    report["promoted_count"] = len(promoted)
    by_type = defaultdict(lambda: [0, 0])
    for r in records:
        if r["tag"] == "real":
            base = r["type"].split(":")[0]
            by_type[base][0] += 1
            by_type[base][1] += int(r["promoted"])
    report["promotion_rate_by_type"] = {
        k: {"scored": v[0], "promoted": v[1]} for k, v in by_type.items()}
    report["eis_promoted"] = any(p["l"] == 5 and p["k"] >= 100
                                 for p in promoted)
    print(f"promoted {len(promoted)} EIS={report['eis_promoted']} "
          f"depths={report['chain_depth_hist']} ({time.time()-t0:.0f}s)")

    # ---------- shuffled loop ----------
    print("== shuffled ==")
    rngs = np.random.default_rng(SEED + 1)
    arr = np.empty((SA.N, C), dtype=object)
    for i, r in enumerate(sa_rows):
        arr[i] = r
    for c in range(C):
        arr[:, c] = arr[rngs.permutation(SA.N), c]
    SAsh = p05.Fam("SA", PRIMES, [list(arr[i]) for i in range(SA.N)],
                   sa_labels)
    cache_sh = sa_cluster_cache(SAsh)
    sh = run_loop(SAsh, cache_sh, records, tag="shuffled")
    report["shuffled_rounds"] = sh["rounds"]
    report["shuffled_promoted"] = len(sh["promoted"])
    print(f"shuffled promoted {len(sh['promoted'])} ({time.time()-t0:.0f}s)")

    # ---------- N2 orderings (cache reused; scans not repeated) ----------
    print("== N2 ==")
    sets, bitsv = [], []
    for r_ in range(ORDERINGS):
        res = run_loop(SA, cache, [], order_seed=SEED + r_, tag=f"o{r_}")
        sets.append(res["final_set"])
        bitsv.append(res["cum_dE"])
    jac = []
    for i in range(ORDERINGS):
        for j in range(i + 1, ORDERINGS):
            u = sets[i] | sets[j]
            jac.append(len(sets[i] & sets[j]) / len(u) if u else 1.0)
    spread = [abs(bitsv[i] - bitsv[j]) for i in range(ORDERINGS)
              for j in range(i + 1, ORDERINGS)]
    report["N2"] = {"jaccard_min": min(jac),
                    "jaccard_mean": float(np.mean(jac)),
                    "floor_q95_bits": float(np.quantile(spread, 0.95))}
    print(report["N2"], f"({time.time()-t0:.0f}s)")

    # ---------- verdict ----------
    bt = real["bits_traj"]
    decreasing = all(b2 <= b1 + 1e-9 for b1, b2 in zip(bt, bt[1:]))
    max_depth = max(real["chain_depths"], default=0)
    if not report["H"]["pass"]:
        verdict = "VOID"
    elif report["shuffled_promoted"] > 0:
        verdict = "FATAL"
    elif not report["eis_promoted"]:
        verdict = "SCORER-DEAD"
    elif max_depth == 0:
        verdict = "NO-COMPOUND"
    elif not decreasing:
        verdict = "SATURATED"
    else:
        verdict = "COMPOUNDING"
    report["verdict"] = verdict
    report["runtime_s"] = round(time.time() - t0, 1)

    with open(os.path.join(args.out, "loop.jsonl"), "w") as fh:
        for r in records:
            fh.write(json.dumps(r, default=str) + "\n")
    with open(os.path.join(args.out, "quarantine.jsonl"), "w") as fh:
        for k in sorted(real["quarantine"]):
            fh.write(json.dumps({"key": k}) + "\n")
    with open(os.path.join(args.out, "summary.json"), "w") as fh:
        json.dump(report, fh, indent=1, default=str)
    print(f"\n== VERDICT: {verdict} == ({report['runtime_s']}s)")


if __name__ == "__main__":
    main()
