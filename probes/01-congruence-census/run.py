#!/usr/bin/env python3
"""Probe 01 — mod-l congruence census over Cremona elliptic curves, N < 10000.

Mechanical implementation of PREREG.md (committed before this file was run).
Cells: H1 harness, H2 null calibration, A1 census, A2 conductor structure,
A3 depth. Verdict is read off the preregistered ladder with no discretion.

Usage:
    python3 run.py --ecdata /path/to/ecdata --out output/

Requires: numpy. Data: ecdata @ 25cec5ecfec8b9f016eb1631ac633194c2bed39f,
files aplist/aplist.00000-09999 and allcurves/allcurves.00000-09999.
"""

import argparse
import json
import math
import os
import subprocess
import sys
import time
from collections import defaultdict

import numpy as np

PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
          67, 71, 73, 79, 83, 89, 97]  # the 25 aplist columns
MODULI = [5, 7, 11, 13]
SEED = 20260728
MIN_SHARED = 10          # minimum shared good positions for an eligible pair
EXT_LIMIT = 1000         # A3 extension: primes 97 < p < EXT_LIMIT
PAIR_ENUM_CAP = 5_000_000
PAIR_SAMPLE = 1_000_000
PINNED_SHA = "25cec5ecfec8b9f016eb1631ac633194c2bed39f"


def primes_upto(n):
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    return [int(p) for p in np.nonzero(sieve)[0]]


EXT_PRIMES = [p for p in primes_upto(EXT_LIMIT - 1) if p > 97]  # 143 primes


class DSU:
    def __init__(self, n):
        self.parent = list(range(n))

    def find(self, x):
        p = self.parent
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


# ---------------------------------------------------------------- parsing

def parse_aplist(path):
    """Return labels[(N,cls)], N array, vals int8 (bad->0), good bool mask."""
    labels, Ns, rows, goods = [], [], [], []
    with open(path) as fh:
        for line in fh:
            tok = line.split()
            if not tok:
                continue
            if len(tok) != 2 + len(PRIMES):
                raise ValueError(f"bad aplist line (got {len(tok)} tokens): {line!r}")
            N, cls = int(tok[0]), tok[1]
            vals, good = [], []
            for t in tok[2:]:
                if t in ("+", "-", "?"):
                    vals.append(0)
                    good.append(False)
                else:
                    vals.append(int(t))
                    good.append(True)
            labels.append((N, cls))
            Ns.append(N)
            rows.append(vals)
            goods.append(good)
    return (labels, np.array(Ns, dtype=np.int64),
            np.array(rows, dtype=np.int8), np.array(goods, dtype=bool))


def parse_allcurves(path):
    """Return dict (N, cls) -> list of (number, ainvs)."""
    out = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            tok = line.split()
            if not tok:
                continue
            N, cls, num = int(tok[0]), tok[1], int(tok[2])
            ainvs = json.loads(tok[3])
            out[(N, cls)].append((num, ainvs))
    return out


# ------------------------------------------------------- point counting

_QR_CACHE = {}


def _qr_table(p):
    tab = _QR_CACHE.get(p)
    if tab is None:
        x = np.arange(p, dtype=np.int64)
        tab = np.full(p, -1, dtype=np.int8)
        tab[(x * x) % p] = 1
        tab[0] = 0
        _QR_CACHE[p] = tab
    return tab


def ap_point_count(ainvs, p):
    """a_p of E: y^2 + a1 xy + a3 y = x^3 + a2 x^2 + a4 x + a6 over F_p.

    Only valid at good primes (p not dividing the conductor). For odd p the
    substitution Y = 2y + a1 x + a3 gives Y^2 = 4x^3 + b2 x^2 + 2 b4 x + b6,
    so a_p = -sum_x chi(f(x)).
    """
    a1, a2, a3, a4, a6 = ainvs
    if p == 2:
        naff = 0
        for x in (0, 1):
            for y in (0, 1):
                if (y * y + a1 * x * y + a3 * y
                        - (x ** 3 + a2 * x * x + a4 * x + a6)) % 2 == 0:
                    naff += 1
        return 2 - naff
    a1, a2, a3, a4, a6 = a1 % p, a2 % p, a3 % p, a4 % p, a6 % p
    b2 = (a1 * a1 + 4 * a2) % p
    b4 = (2 * a4 + a1 * a3) % p
    b6 = (a3 * a3 + 4 * a6) % p
    x = np.arange(p, dtype=np.int64)
    fx = (4 * x ** 3 + b2 * x ** 2 + 2 * b4 * x + b6) % p
    return -int(_qr_table(p)[fx].sum())


# ------------------------------------------------------------- the scan

def scan(vals, good, min_shared=MIN_SHARED, tag="", log=print):
    """Exhaustive pairwise scan: edge iff equal at every shared good position
    and >= min_shared shared good positions. Returns (dsu, stats)."""
    n, K = vals.shape
    dsu = DSU(n)
    n_edges = 0
    n_ineligible = 0
    min_sh = K + 1
    t0 = time.time()
    for i in range(n - 1):
        gi = good[i]
        vi = vals[i]
        g = good[i + 1:]
        v = vals[i + 1:]
        both = g & gi
        neq = v != vi
        neq &= both
        cand = ~neq.any(axis=1)
        if cand.any():
            js = np.nonzero(cand)[0]
            sh = both[js].sum(axis=1)
            for j, s in zip(js.tolist(), sh.tolist()):
                if s >= min_shared:
                    dsu.union(i, i + 1 + j)
                    n_edges += 1
                    if s < min_sh:
                        min_sh = s
                else:
                    n_ineligible += 1
        if i % 10000 == 0 and i:
            log(f"    [{tag}] anchor {i}/{n} edges={n_edges} "
                f"({time.time() - t0:.0f}s)")
    comps = defaultdict(list)
    for i in range(n):
        comps[dsu.find(i)].append(i)
    nontrivial = sorted((m for m in comps.values() if len(m) > 1),
                        key=lambda m: (-len(m), m[0]))
    stats = {"edges": n_edges, "ineligible_pairs": n_ineligible,
             "min_shared_on_edges": (None if min_sh > K else min_sh),
             "seconds": round(time.time() - t0, 1)}
    return nontrivial, stats


# ------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ecdata", required=True)
    ap.add_argument("--out", default="output")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    report = {"prereg_seed": SEED, "moduli": MODULI, "min_shared": MIN_SHARED,
              "primes_parsed": len(PRIMES), "primes_used_per_modulus": len(PRIMES) - 1,
              "ext_primes": len(EXT_PRIMES), "ext_limit": EXT_LIMIT}

    try:
        sha = subprocess.run(["git", "-C", args.ecdata, "rev-parse", "HEAD"],
                             capture_output=True, text=True).stdout.strip()
    except OSError:
        sha = "unknown"
    report["ecdata_sha"] = sha
    if sha != PINNED_SHA:
        print(f"WARNING: ecdata HEAD {sha} != pinned {PINNED_SHA}")

    print("== load ==")
    labels, Ns, vals, good = parse_aplist(
        os.path.join(args.ecdata, "aplist", "aplist.00000-09999"))
    curves = parse_allcurves(
        os.path.join(args.ecdata, "allcurves", "allcurves.00000-09999"))
    n = len(labels)
    lab = {t: i for i, t in enumerate(labels)}
    ncurves = sum(len(v) for v in curves.values())
    all_good = int((good.sum(axis=1) == len(PRIMES)).sum())
    report["n_classes"] = n
    report["n_curves"] = ncurves
    report["classes_with_no_bad_prime_below_97"] = all_good
    print(f"classes={n} curves={ncurves} fully-good classes={all_good}")

    # ---------------- H1 ----------------
    print("== H1 harness ==")
    h1 = {}
    h1["join_curves_without_aplist_row"] = sorted(
        f"{N}{c}" for (N, c) in curves if (N, c) not in lab)
    h1["join_aplist_rows_without_curve"] = sorted(
        f"{N}{c}" for (N, c) in labels if (N, c) not in curves)
    h1["dup_aplist_rows"] = n - len(lab)
    h1["curve1_missing"] = sorted(
        f"{N}{c}" for (N, c), lst in curves.items()
        if 1 not in {num for num, _ in lst})
    h1_join_ok = (not h1["join_curves_without_aplist_row"]
                  and not h1["join_aplist_rows_without_curve"]
                  and h1["dup_aplist_rows"] == 0 and not h1["curve1_missing"])

    exact_comps, exact_stats = scan(vals, good, tag="H1-exact")
    comp_of = {}
    for k, members in enumerate(exact_comps):
        for m in members:
            comp_of[m] = k
    # every isogeny class inside a single group: walk the curve table, map
    # each curve through the join to its fingerprint row and component, and
    # require one component per class.
    class_comp = defaultdict(set)
    for (N, c), lst in curves.items():
        i = lab.get((N, c))
        if i is None:
            continue
        for _num, _ainvs in lst:
            class_comp[(N, c)].add(comp_of.get(i, -1 - i))  # singleton -> unique id
    split_classes = sorted(f"{N}{c}" for (N, c), s in class_comp.items()
                           if len(s) != 1)
    h1["split_isogeny_classes"] = split_classes
    h1["cross_class_exact_agreement_pairs"] = exact_stats["edges"]
    h1["cross_class_exact_agreement_components"] = [
        [f"{labels[m][0]}{labels[m][1]}" for m in members]
        for members in exact_comps]

    rng = np.random.default_rng(SEED)
    sample = rng.choice(n, size=200, replace=False)
    mismatches = []
    for i in sample.tolist():
        N, c = labels[i]
        ainvs = dict(curves[(N, c)])[1]
        for k, p in enumerate(PRIMES):
            if not good[i, k]:
                continue
            computed = ap_point_count(ainvs, p)
            if computed != int(vals[i, k]):
                mismatches.append((f"{N}{c}", p, computed, int(vals[i, k])))
    h1["pointcount_sample_size"] = 200
    h1["pointcount_mismatches"] = mismatches
    h1_pass = h1_join_ok and not split_classes and not mismatches
    h1["pass"] = h1_pass
    report["H1"] = h1
    print(f"H1 pass={h1_pass} join_ok={h1_join_ok} split={len(split_classes)} "
          f"pointcount_mismatches={len(mismatches)} "
          f"exact_agreement_pairs={exact_stats['edges']}")

    # ---------------- A1 + H2 ----------------
    a1 = {}
    h2 = {}
    components = {}
    for ell in MODULI:
        kcol = PRIMES.index(ell)
        cols = [k for k in range(len(PRIMES)) if k != kcol]
        vm = (vals.astype(np.int16) % ell).astype(np.int8)[:, cols]
        gm = good[:, cols]
        print(f"== A1 census mod {ell} ==")
        comps, stats = scan(vm, gm, tag=f"A1-l{ell}")
        components[ell] = comps
        a1[ell] = {"nontrivial_components": len(comps),
                   "classes_involved": sum(len(m) for m in comps),
                   "sizes": sorted((len(m) for m in comps), reverse=True),
                   **stats}
        print(f"mod {ell}: {len(comps)} non-trivial components, "
              f"{a1[ell]['classes_involved']} classes, "
              f"largest={a1[ell]['sizes'][:5]}")

        print(f"== H2 null mod {ell} ==")
        rng = np.random.default_rng(SEED + ell)
        vsh = (vals.astype(np.int16) % ell).astype(np.int8).copy()
        gsh = good.copy()
        for k in range(len(PRIMES)):  # permute all 25 (value, flag) columns
            perm = rng.permutation(n)
            vsh[:, k] = vsh[perm, k]
            gsh[:, k] = gsh[perm, k]
        comps_sh, stats_sh = scan(vsh[:, cols], gsh[:, cols], tag=f"H2-l{ell}")
        h2[ell] = {"nontrivial_components": len(comps_sh),
                   "sizes": sorted((len(m) for m in comps_sh), reverse=True),
                   **stats_sh}
        print(f"H2 mod {ell}: {len(comps_sh)} non-trivial components")
    report["A1"] = {str(k): v for k, v in a1.items()}
    report["H2"] = {str(k): v for k, v in h2.items()}

    # ---------------- A3 depth machinery ----------------
    print("== A3 extension: point counts for members of surviving components ==")
    needed = sorted({m for ell in MODULI for comp in components[ell] for m in comp})
    ext = {}
    t0 = time.time()
    for cnt, i in enumerate(needed):
        N, c = labels[i]
        ainvs = dict(curves[(N, c)])[1]
        ext[i] = {p: ap_point_count(ainvs, p) for p in EXT_PRIMES if N % p != 0}
        if cnt % 500 == 0 and cnt:
            print(f"    ext {cnt}/{len(needed)} ({time.time() - t0:.0f}s)")
    report["A3_ext_classes"] = len(needed)

    # per-member test vector over the combined prime list, per modulus
    def member_vec(ell, i, cols_primes):
        v = np.full(len(cols_primes) + len(EXT_PRIMES), -1, dtype=np.int16)
        for k, p in enumerate(cols_primes):
            kk = PRIMES.index(p)
            if good[i, kk]:
                v[k] = int(vals[i, kk]) % ell
        for k, p in enumerate(EXT_PRIMES):
            a = ext[i].get(p)
            if a is not None:
                v[len(cols_primes) + k] = a % ell
        return v

    # ---------------- A2 + A3 pair evaluation ----------------
    print("== A2/A3 pair evaluation ==")
    pooled = {"persistent": 0, "consistent": 0,
              "by_class": defaultdict(int), "breaking": 0}
    groups_rows = []
    anomalies = []
    breaking_rows = []

    for ell in MODULI:
        cols_primes = [p for p in PRIMES if p != ell]
        all_primes = cols_primes + EXT_PRIMES
        pplus1 = np.array([(p + 1) % ell for p in all_primes], dtype=np.int16)
        # explicit pair budget check (prereg: cap 5e6, else sample 1e6)
        explicit_total = 0
        vecs, eis = {}, {}
        for comp in components[ell]:
            for m in comp:
                v = member_vec(ell, m, cols_primes)
                vecs[m] = v
                have = v >= 0
                eis[m] = bool((v[have] == pplus1[have]).all())
            ne = sum(1 for m in comp if not eis[m])
            k = len(comp)
            explicit_total += k * (k - 1) // 2 - (k - ne) * (k - ne - 1) // 2
        sampling = explicit_total > PAIR_ENUM_CAP
        keep_prob = (PAIR_SAMPLE / explicit_total) if sampling else 1.0
        if sampling:
            print(f"    mod {ell}: explicit pairs {explicit_total} > cap; "
                  f"Bernoulli-sampling ~{PAIR_SAMPLE} (preregistered fallback)")
        srng = np.random.default_rng(SEED)

        for comp in components[ell]:
            k = len(comp)
            eis_members = [m for m in comp if eis[m]]
            n_eis_pairs = len(eis_members) * (len(eis_members) - 1) // 2
            cls_count = defaultdict(int)
            cls_count["EIS"] = n_eis_pairs           # persistent+consistent by construction
            n_persist = n_eis_pairs
            n_consist = n_eis_pairs
            comp_break = []

            def explicit_pairs(comp=comp):
                for ai, a in enumerate(comp):
                    for b in comp[ai + 1:]:
                        if not (eis[a] and eis[b]):
                            yield a, b

            for a, b in explicit_pairs():
                if sampling and srng.random() > keep_prob:
                    continue
                va, vb = vecs[a], vecs[b]
                diff = (va != vb) & (va >= 0) & (vb >= 0)
                if diff.any():
                    bp = all_primes[int(np.argmax(diff))]
                    comp_break.append(bp)
                    breaking_rows.append((ell,
                                          f"{labels[a][0]}{labels[a][1]}",
                                          f"{labels[b][0]}{labels[b][1]}", bp))
                    continue
                n_persist += 1
                N1, N2 = labels[a][0], labels[b][0]
                g = math.gcd(N1, N2)
                if N1 == N2:
                    cl = "P1"
                elif max(N1, N2) % min(N1, N2) == 0:
                    cl = "P2"
                elif g > 1:
                    cl = "P3"
                else:
                    cl = "P4"
                cls_count[cl] += 1
                if cl != "P4":
                    n_consist += 1
                else:
                    anomalies.append((ell,
                                      f"{labels[a][0]}{labels[a][1]}", N1,
                                      f"{labels[b][0]}{labels[b][1]}", N2,
                                      g, N1 // g, N2 // g))
            pooled["persistent"] += n_persist
            pooled["consistent"] += n_consist
            pooled["breaking"] += len(comp_break)
            for key, val in cls_count.items():
                pooled["by_class"][key] += val
            members_lbl = [f"{labels[m][0]}{labels[m][1]}" for m in comp]
            conds = sorted({labels[m][0] for m in comp})
            groups_rows.append({
                "l": ell, "size": k, "n_eis_members": len(eis_members),
                "members": members_lbl[:20] + (
                    [f"...+{k - 20} more"] if k > 20 else []),
                "conductors": conds[:20] + (
                    [f"...+{len(conds) - 20} more"] if len(conds) > 20 else []),
                "pair_classes": dict(cls_count),
                "persistent_pairs": n_persist,
                "breaking_pairs": len(comp_break),
                "min_break_prime": min(comp_break) if comp_break else None})
        report.setdefault("A2A3_sampling", {})[str(ell)] = sampling

    pooled["by_class"] = dict(pooled["by_class"])
    consistency = (pooled["consistent"] / pooled["persistent"]
                   if pooled["persistent"] else None)
    pooled["consistency_fraction"] = consistency
    report["pooled"] = pooled

    # ---------------- verdict (prereg §7, in order) ----------------
    void_reasons = []
    if not h1_pass:
        void_reasons.append("H1 failed")
    for ell in MODULI:
        h2l, a1l = h2[ell]["nontrivial_components"], a1[ell]["nontrivial_components"]
        if h2l >= 3:
            void_reasons.append(f"H2 mod {ell} = {h2l} >= 3")
        elif h2l >= 1 and 4 * h2l >= a1l:
            void_reasons.append(f"H2 mod {ell} = {h2l} comparable to A1 = {a1l}")
    total_a1 = sum(a1[ell]["nontrivial_components"] for ell in MODULI)
    if void_reasons:
        verdict = "VOID"
    elif total_a1 == 0:
        verdict = "KNOWN-ONLY"
    elif pooled["persistent"] == 0:
        verdict = "CLOSED-INCONCLUSIVE"
    elif consistency >= 0.90:
        verdict = "VALIDATED"
    else:
        verdict = "OPEN"
    report["void_reasons"] = void_reasons
    report["verdict"] = verdict

    # ---------------- outputs ----------------
    def jsonable(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        raise TypeError(o)

    with open(os.path.join(args.out, "summary.json"), "w") as fh:
        json.dump(report, fh, indent=1, default=jsonable)
    with open(os.path.join(args.out, "groups.csv"), "w") as fh:
        fh.write("l,size,n_eis_members,persistent_pairs,breaking_pairs,"
                 "min_break_prime,pair_classes,conductors,members\n")
        for r in groups_rows:
            fh.write(f"{r['l']},{r['size']},{r['n_eis_members']},"
                     f"{r['persistent_pairs']},{r['breaking_pairs']},"
                     f"{r['min_break_prime']},"
                     f"\"{r['pair_classes']}\",\"{r['conductors']}\","
                     f"\"{r['members']}\"\n")
    with open(os.path.join(args.out, "anomalies.csv"), "w") as fh:
        fh.write("l,class1,N1,class2,N2,gcd,N1_over_gcd,N2_over_gcd\n")
        for row in anomalies:
            fh.write(",".join(map(str, row)) + "\n")
    with open(os.path.join(args.out, "breaking_pairs.csv"), "w") as fh:
        fh.write("l,class1,class2,break_prime\n")
        for row in breaking_rows:
            fh.write(",".join(map(str, row)) + "\n")

    print("\n== VERDICT ==")
    print(f"verdict: {verdict}")
    if void_reasons:
        print("void reasons:", void_reasons)
    print(f"A1 per modulus: "
          + ", ".join(f"l={e}: {a1[e]['nontrivial_components']}" for e in MODULI))
    print(f"H2 per modulus: "
          + ", ".join(f"l={e}: {h2[e]['nontrivial_components']}" for e in MODULI))
    print(f"persistent pairs={pooled['persistent']} "
          f"consistent={pooled['consistent']} "
          f"fraction={consistency} breaking={pooled['breaking']}")
    print(f"pair classes: {pooled['by_class']}")


if __name__ == "__main__":
    main()
