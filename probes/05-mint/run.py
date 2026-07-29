#!/usr/bin/env python3
"""Probe 05 — ANCHOR minting on S-A (EC/Q) and S-B (Artin-joined).

Usage:
    python3 run.py --data data --p03data ../03-transport/data \
        --ecdata /path/to/ecdata --out output

Implements PREREG.md verbatim. Known disclosed defect, implemented as
written per the binding constraints: identifiability rule (a) is satisfied
by an anchor's own members (each member's good columns reduce to r_h by
construction), so K3 is structurally empty under this prereg's definition;
RESULTS records this as a registry defect for any successor prereg.
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
EPS = 6
MIN_SHARED = 10
MINT_NORM = 43                       # alphabet bound 13 >= l-1 for all moduli
MINT_CELL_BITS = math.log2(2 * 13 + 1 + 3)   # uniform code, 30 symbols
REPLICATES = 20


def load_module(name, relpath):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, relpath))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


p01 = load_module("p01r", "../01-congruence-census/run.py")
p02 = load_module("p02r", "../02-encoder-energy/run.py")
p03 = load_module("p03r", "../03-transport/run.py")


# ---------------------------------------------------------------- helpers

class Fam:
    """Family of rows with values (int|None) and a p02 Table."""

    def __init__(self, name, norms, rows, labels):
        self.name = name
        self.norms = list(norms)
        self.rows = rows
        self.labels = labels
        self.N = len(rows)
        with p03.p02_norms(self.norms):
            self.table = p02.Table(name, [[("?" if v is None else str(int(v)))
                                           for v in r] for r in rows])

    def val_good(self):
        V = np.array([[0 if v is None else v for v in r] for r in self.rows],
                     dtype=np.int64)
        G = np.array([[v is not None for v in r] for r in self.rows])
        return V, G


def uniformize(table):
    """MINT-block coding: uniform per column, model cost 0 (prereg §1)."""
    for c in range(C):
        n = table.nsym[c]
        table.q[c] = np.full(n, 1.0 / n)
        table.logq_sym[c] = np.full(n, math.log2(n))
        table.logq_cell[:, c] = math.log2(n)
    table.model_bits_col = [0.0] * C
    table.L0_model = 0.0
    table.L0_data = float(table.N * sum(math.log2(table.nsym[c])
                                        for c in range(C)))
    table.L0 = table.L0_data
    for l in MODULI:
        Ql = np.zeros((C, l))
        for c in range(C):
            B = (table.nsym[c] - 3 - 1) // 2
            for i, v in enumerate(range(-B, B + 1)):
                Ql[c, v % l] += table.q[c][i]
        table.Q[l] = Ql
        with np.errstate(divide="ignore"):
            nl = -np.log2(Ql)
        table.negLogQ[l] = np.where(np.isfinite(nl), nl, 0.0)
    return table


def clusters(fams, l):
    """Probe-01 components over the concatenation of families, modulus l."""
    kcol = PRIMES.index(l)
    cols = [k for k in range(C) if k != kcol]
    Vs, Gs, owners = [], [], []
    for f in fams:
        V, G = f.val_good()
        Vs.append(V)
        Gs.append(G)
        owners += [(f.name, i) for i in range(f.N)]
    V = np.concatenate(Vs) if Vs else np.zeros((0, C), dtype=np.int64)
    G = np.concatenate(Gs) if Gs else np.zeros((0, C), dtype=bool)
    vm = (V % l).astype(np.int8)[:, cols]
    gm = G[:, cols]
    comps, _ = p01.scan(vm, gm, tag=f"l{l}", log=lambda *a, **k: None)
    return [(l, [owners[i] for i in comp]) for comp in comps]


class AnchorEngine:
    """Sequential ΔE acceptance with cross-modulus ownership + reuse pass."""

    def __init__(self, fams, log2N0):
        self.fams = {f.name: f for f in fams}
        self.log2N0 = log2N0

    def residue_vec(self, members, l):
        rh = np.full(C, 0, dtype=np.int16)
        kcol = PRIMES.index(l)
        for c in range(C):
            if c == kcol:
                rh[c] = -1
                continue
            for (fn, i) in members:
                v = self.fams[fn].rows[i][c]
                if v is not None:
                    rh[c] = v % l
                    break
        return rh

    def anchor_delta(self, members, l, rh, owned, charge_row=True,
                     claim=False):
        fb = p02.flag_bits(EPS)
        k = len(members)
        cost = (p02.H_BITS + p02.gamma_bits(k) + k * self.log2N0
                + 24 * math.log2(l))
        if charge_row:
            cost += 25 * MINT_CELL_BITS
        delta = 0.0
        newly = []
        for (fn, i) in members:
            t = self.fams[fn].table
            rows_f = self.fams[fn].rows
            for c in range(C):
                if rh[c] < 0 or (fn, i, c) in owned:
                    continue
                v = rows_f[i][c]
                if v is not None and v % l == rh[c]:
                    delta += fb - t.negLogQ[l][c, rh[c]]
                else:
                    delta += EPS
                newly.append((fn, i, c))
        if claim:
            owned.update(newly)
        return cost + delta

    def run(self, cands, order_rng=None):
        """cands: list of (l, members). Returns accepted list + mint rows."""
        order = sorted(range(len(cands)),
                       key=lambda t: (cands[t][0], -len(cands[t][1]),
                                      cands[t][1][0]))
        if order_rng is not None:
            order = list(order_rng.permutation(len(cands)))
        owned = set()
        accepted = []
        for idx in order:
            l, members = cands[idx]
            if len(members) < 2:
                continue
            rh = self.residue_vec(members, l)
            dE = self.anchor_delta(members, l, rh, owned, claim=False)
            if dE < 0:
                self.anchor_delta(members, l, rh, owned, claim=True)
                accepted.append({"l": l, "k": len(members), "dE": dE,
                                 "rh": rh, "members": members,
                                 "id": hashlib.sha256(
                                     (str(l) + str(sorted(members)))
                                     .encode()).hexdigest()[:16]})
        return accepted


def reuse_pass(accepted, log2N0):
    """Menu: DUP among mint rows; PATTERN over MINT; chained ANCHOR."""
    n = len(accepted)
    reuse = [0] * n
    chain_depth = [0] * n
    if n == 0:
        return reuse, chain_depth, {"dups": 0, "patterns": 0, "chained": 0}
    rows = [[int(a["rh"][c]) if a["rh"][c] >= 0 else None for c in range(C)]
            for a in accepted]
    mint = Fam("MINT", [MINT_NORM] * C, rows, [a["id"] for a in accepted])
    uniformize(mint.table)
    fb = p02.flag_bits(EPS)
    stats = {"dups": 0, "patterns": 0, "chained": 0}
    # DUP among mint rows
    for j in range(n):
        for i in range(j):
            if rows[i] == rows[j]:
                dE = (p02.H_BITS + 2 * log2N0
                      + sum(fb - MINT_CELL_BITS if rows[j][c] is not None
                            else fb - MINT_CELL_BITS for c in range(C)))
                if dE < 0:
                    reuse[i] += 1
                    reuse[j] += 1
                    stats["dups"] += 1
                break
    # PATTERN over MINT block
    for pi in p02.PI_ORDER:
        fn_ = p02.PI_FUN[pi]
        for l in MODULI:
            members = []
            for j in range(n):
                ok = any(rows[j][c] is not None for c in range(C))
                for c in range(C):
                    v = rows[j][c]
                    if v is not None and v % l != fn_(PRIMES[c]) % l:
                        ok = False
                        break
                if ok:
                    members.append(j)
            if len(members) < 2:
                continue
            k = len(members)
            savings = 0.0
            for j in members:
                for c in range(C):
                    if rows[j][c] is not None:
                        res = fn_(PRIMES[c]) % l
                        savings += mint.table.negLogQ[l][c, res] - fb
            dE = (p02.H_BITS + 2 + p02.gamma_bits(k) + k * log2N0) - savings
            if dE < 0:
                stats["patterns"] += 1
                for j in members:
                    reuse[j] += 1
    # chained ANCHOR over MINT clusters
    eng = AnchorEngine([mint], log2N0)
    cands = []
    for l in MODULI:
        cands += clusters([mint], l)
    chained = eng.run([(l, m) for (l, m) in cands if len(m) >= 2])
    stats["chained"] = len(chained)
    for a in chained:
        for (_fn, j) in a["members"]:
            reuse[j] += 1
            chain_depth[j] = max(chain_depth[j], 1)
    return reuse, chain_depth, stats


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--p03data", default=os.path.join(HERE, "..",
                                                      "03-transport", "data"))
    ap.add_argument("--ecdata", required=True)
    ap.add_argument("--out", default="output")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    report = {"seed": SEED,
              "disclosed_defect": ("identifiability rule (a) is satisfied by "
                                   "an anchor's own members; K3 is "
                                   "structurally empty under this prereg — "
                                   "registry defect for any successor")}

    # ---------- S-A ----------
    print("== load S-A ==")
    apl, _ = p02.parse_aplist(os.path.join(args.ecdata, "aplist",
                                           "aplist.00000-09999"))
    sa_rows = [[None if s in ("+", "-", "?") else int(s) for s in syms]
               for (_N, _c, syms, _t) in apl]
    sa_labels = [f"{N}{c}" for (N, c, _s, _t) in apl]
    SA = Fam("SA", PRIMES, sa_rows, sa_labels)

    # ---------- S-B ----------
    print("== load S-B ==")
    art_rows_csv = p03.read_csv(os.path.join(args.data, "artin.csv"))
    art_rows, art_labels, n_noL, n_sealed = [], [], 0, 0
    import re as _re
    for r in art_rows_csv:
        if int(r["sealed"]):
            n_sealed += 1
            continue
        ef = json.loads(r["euler_factors"])
        if ef is None:
            n_noL += 1
            continue
        bad = set(json.loads(r["BadPrimes"]) or []) \
            | set(json.loads(r["HardPrimes"]) or [])
        row = []
        for i, p_ in enumerate(PRIMES):
            if p_ in bad or i >= len(ef) or not ef[i] or len(ef[i]) < 2:
                row.append(None)
            else:
                row.append(-ef[i][1])
        art_rows.append(row)
        art_labels.append(r["Baselabel"])
    ART = Fam("ART", [1] * C, art_rows, art_labels)

    f1_rows_csv = p03.read_csv(os.path.join(args.p03data, "f1_ec_q.csv"))
    f1_classes = {}
    for r in f1_rows_csv:
        m = _re.match(r"^(\d+)\.([a-z]+)(\d+)$", r["lmfdb_label"])
        if m and int(m.group(3)) == 1:
            f1_classes[f"{m.group(1)}.{m.group(2)}"] = r
    f1_rows, f1_labels = [], []
    for cls in sorted(f1_classes, key=lambda s: (int(s.split(".")[0]),
                                                 s.split(".")[1])):
        r = f1_classes[cls]
        N = int(r["conductor"])
        ain = json.loads(r["ainvs"])
        f1_rows.append([None if N % p_ == 0 else p01.ap_point_count(ain, p_)
                        for p_ in PRIMES])
        f1_labels.append(cls)
    F1 = Fam("F1", PRIMES, f1_rows, f1_labels)

    f2_rows_csv = p03.read_csv(os.path.join(args.p03data, "f2_newforms.csv"))
    f2_rows, f2_labels = [], []
    for r in sorted(f2_rows_csv, key=lambda r: (int(r["level"]), r["label"])):
        if int(r["sealed"]):
            continue
        if int(r["dim"]) == 1:
            f2_rows.append([int(r[f"a{p_}"]) for p_ in PRIMES])
        else:
            f2_rows.append([None] * C)
        f2_labels.append(r["label"])
    F2 = Fam("F2", PRIMES, f2_rows, f2_labels)

    SB = [ART, F1, F2]
    report["rows"] = {"SA": SA.N, "ART": ART.N, "F1": F1.N, "F2": F2.N,
                      "artin_no_L": n_noL, "artin_sealed": n_sealed}
    print(report["rows"], f"({time.time()-t0:.0f}s)")

    # ---------- H1 ----------
    print("== H1 ==")
    rng = np.random.default_rng(SEED)
    h1 = {"sa_mismatches": 0, "artin_bound_violations": 0,
          "artin_factor_mismatches": 0, "artin_rows_checked": 0}
    # S-A: 200 rows re-counted from Cremona curve-1 models
    curves = p01.parse_allcurves(os.path.join(args.ecdata, "allcurves",
                                              "allcurves.00000-09999"))
    for i in rng.choice(SA.N, 200, replace=False).tolist():
        N, ccls = apl[i][0], apl[i][1]
        ain = dict(curves[(N, ccls)])[1]
        for k, p_ in enumerate(PRIMES):
            v = sa_rows[i][k]
            if v is not None and p01.ap_point_count(ain, p_) != v:
                h1["sa_mismatches"] += 1
    # Artin: full-table dim bound; 200 rows factor-membership check
    for row in art_rows:
        for v in row:
            if v is not None and abs(v) > 2:
                h1["artin_bound_violations"] += 1
    lf_by_label = {r["Baselabel"]: json.loads(r["LocalFactors"])
                   for r in art_rows_csv}
    art_idx = rng.choice(ART.N, min(200, ART.N), replace=False).tolist()
    for i in art_idx:
        lf = lf_by_label.get(art_labels[i])
        if not lf:
            continue
        allowed = set()
        for fac in lf:
            if len(fac) >= 2:
                c1 = fac[1][0] if isinstance(fac[1], list) else fac[1]
                allowed.add(-c1)
        h1["artin_rows_checked"] += 1
        for v in ART.rows[i]:
            if v is not None and v not in allowed:
                h1["artin_factor_mismatches"] += 1
    h1["pass"] = (h1["sa_mismatches"] == 0
                  and h1["artin_bound_violations"] == 0
                  and h1["artin_factor_mismatches"] == 0)
    report["H1"] = h1
    print(h1, f"({time.time()-t0:.0f}s)")

    # ---------- candidates ----------
    print("== clusters ==")
    cands_A = []
    for l in MODULI:
        cands_A += [(l, m) for (l, m) in clusters([SA], l) if len(m) >= 2]
    cands_B = []
    for l in MODULI:
        cands_B += [(l, m) for (l, m) in clusters(SB, l) if len(m) >= 2]
    report["K5"] = {"candidates_SA": len(cands_A),
                    "candidates_SB": len(cands_B),
                    "case": "SAMPLE-SIZE" if len(cands_B) < 500 else "ADEQUATE"}
    print(report["K5"], f"({time.time()-t0:.0f}s)")

    # ---------- accept pass ----------
    print("== anchors ==")
    engA = AnchorEngine([SA], math.log2(SA.N))
    accA = engA.run(cands_A)
    reuseA, depthA, statsA = reuse_pass(accA, math.log2(SA.N))
    N0B = sum(f.N for f in SB)
    engB = AnchorEngine(SB, math.log2(N0B))
    accB = engB.run(cands_B)
    reuseB, depthB, statsB = reuse_pass(accB, math.log2(N0B))
    print(f"S-A accepted {len(accA)} (reuse stats {statsA}); "
          f"S-B accepted {len(accB)} ({statsB})")

    # ---------- H2 (decode incl. minted rows) ----------
    print("== H2 ==")
    h2 = {"streams": [], "pass": True}
    mint_rows_B = [[int(a["rh"][c]) if a["rh"][c] >= 0 else None
                    for c in range(C)] for a in accB]
    test_fams = [SA, ART, F1, F2]
    if mint_rows_B:
        MB = Fam("MINTB", [MINT_NORM] * C, mint_rows_B,
                 [a["id"] for a in accB])
        uniformize(MB.table)
        test_fams.append(MB)
    for f in test_fams:
        with p03.p02_norms(f.norms):
            codec = p02.StreamCodec(f.table, [], EPS)
            blob = codec.encode()
            ok = bool((codec.decode(blob) == f.table.SYM).all())
        rel = abs(8 * len(blob) - f.table.L0_data) / f.table.L0_data
        h2["streams"].append({"fam": f.name, "identical": ok,
                              "rel": round(rel, 7)})
        h2["pass"] &= ok and rel <= 0.005
    report["H2"] = h2
    print("H2 pass:", h2["pass"], f"({time.time()-t0:.0f}s)")

    # ---------- H3 exact reduction to probe-02 HUB ----------
    print("== H3 ==")
    k4csv = os.path.join(HERE, "..", "02-encoder-energy", "output",
                         "k4_groups.csv")
    import csv as _csv
    p02rows = list(_csv.DictReader(open(k4csv)))
    comps02 = {l: [sorted(i for (_f, i) in m)
                   for (ll, m) in cands_A if ll == l and len(m) >= 3]
               for l in MODULI}
    mine = []
    for l in MODULI:
        for ms in comps02[l]:
            hub = p02.make_hub(ms, l, SA.table)
            evh = p02.eval_mod_entry(SA.table, hub, EPS)
            mst = [p02.Entry("BRIDGE", l=l, i=ms[0], j=m) for m in ms[1:]]
            evm = [p02.eval_mod_entry(SA.table, e, EPS) for e in mst]
            A = ((hub.cost(SA.table) - evh["S"])
                 - (sum(e.cost(SA.table) for e in mst)
                    - sum(x["S"] for x in evm)))
            Ap = (p02.gamma_bits(len(ms)) + (2 - len(ms)) * SA.table.log2N
                  - 5 * (len(ms) - 1))
            mine.append((l, len(ms), A, Ap))
    inband = sum(1 for (_l, k, A, Ap) in mine if abs(A - Ap) / abs(Ap) <= 0.20)
    frac = inband / len(mine) if mine else 0.0
    csv_pairs = sorted((int(r["l"]), int(r["k"]), round(float(r["A_hub_meas"]), 2))
                       for r in p02rows)
    my_pairs = sorted((l, k, round(A, 2)) for (l, k, A, _ap) in mine)
    h3 = {"groups": len(mine), "in_band_frac": frac,
          "probe02_in_band_frac": 251 / 714,
          "exact_match_to_probe02_csv": csv_pairs == my_pairs,
          "pass": csv_pairs == my_pairs and abs(frac - 251 / 714) < 1e-9}
    report["H3"] = h3
    print({k: v for k, v in h3.items()})

    H_pass = h1["pass"] and h2["pass"] and h3["pass"]

    # ---------- N1 shuffled (FATAL if any anchor) ----------
    print("== N1 ==")
    def shuffled(fam, fid):
        rngs = np.random.default_rng(SEED + fid)
        V, G = fam.val_good()
        rows = [list(r) for r in fam.rows]
        arr = np.empty((fam.N, C), dtype=object)
        for i, r in enumerate(rows):
            arr[i] = r
        for c in range(C):
            arr[:, c] = arr[rngs.permutation(fam.N), c]
        return Fam(fam.name + "sh", fam.norms,
                   [list(arr[i]) for i in range(fam.N)], fam.labels)

    SAsh = shuffled(SA, 1)
    cands_sh = []
    for l in MODULI:
        cands_sh += [(l, m) for (l, m) in clusters([SAsh], l) if len(m) >= 2]
    acc_sh = AnchorEngine([SAsh], math.log2(SA.N)).run(cands_sh)
    SBsh = [shuffled(ART, 2), shuffled(F1, 3), shuffled(F2, 4)]
    cands_shB = []
    for l in MODULI:
        cands_shB += [(l, m) for (l, m) in clusters(SBsh, l) if len(m) >= 2]
    acc_shB = AnchorEngine(SBsh, math.log2(N0B)).run(cands_shB)
    n1 = {"SA_shuffled_accepted": len(acc_sh),
          "SB_shuffled_accepted": len(acc_shB),
          "pass": len(acc_sh) == 0 and len(acc_shB) == 0}
    report["N1"] = n1
    print(n1, f"({time.time()-t0:.0f}s)")

    # ---------- N2 wobble (genuine) ----------
    print("== N2 ==")
    sets, bits = [], []
    for rrep in range(REPLICATES):
        rngr = np.random.default_rng(SEED + rrep)
        acc = AnchorEngine([SA], math.log2(SA.N)).run(cands_A, order_rng=rngr)
        sets.append({a["id"] for a in acc})
        bits.append(sum(a["dE"] for a in acc))
    jac = []
    for i in range(REPLICATES):
        for j in range(i + 1, REPLICATES):
            u = sets[i] | sets[j]
            jac.append(len(sets[i] & sets[j]) / len(u) if u else 1.0)
    spread = [abs(bits[i] - bits[j]) for i in range(REPLICATES)
              for j in range(i + 1, REPLICATES)]
    n2 = {"jaccard_min": min(jac), "jaccard_mean": float(np.mean(jac)),
          "bits_values": [round(b, 1) for b in bits[:5]],
          "bits_spread_q95": float(np.quantile(spread, 0.95)),
          "floor_bits": float(np.quantile(spread, 0.95))}
    report["N2"] = n2
    print(n2, f"({time.time()-t0:.0f}s)")

    # ---------- N3 foils (full path) ----------
    print("== N3 ==")
    sizes_A = [len(m) for (_l, m) in cands_A] or [3]
    sizes_B = [len(m) for (_l, m) in cands_B] or [3]
    rngf = np.random.default_rng(SEED)
    n3_acc = 0
    for fams, eng, sizes, npool in ((["SA"], engA, sizes_A, SA.N),
                                    (["SB"], engB, sizes_B, N0B)):
        owners = ([("SA", i) for i in range(SA.N)] if fams == ["SA"] else
                  [(f.name, i) for f in SB for i in range(f.N)])
        for _ in range(1000):
            l = int(rngf.choice(MODULI))
            k = int(sizes[int(rngf.integers(0, len(sizes)))])
            members = [owners[int(x)] for x in
                       rngf.choice(len(owners), min(k, len(owners)),
                                   replace=False)]
            rh = np.full(C, 0, dtype=np.int16)
            kcol = PRIMES.index(l)
            for c in range(C):
                rh[c] = -1 if c == kcol else int(rngf.integers(0, l))
            dE = eng.anchor_delta(members, l, rh, set(), claim=False)
            if dE < 0:
                n3_acc += 1        # reuse pass irrelevant: dE<0 is necessary
    n3 = {"foils": 2000, "accepted": n3_acc, "pass": n3_acc == 0}
    report["N3"] = n3
    print(n3)

    # ---------- K cells ----------
    eis = max((a for a in accA if a["l"] == 5), key=lambda a: a["k"],
              default=None)
    K1 = {"eis_minted": bool(eis and eis["k"] == 143),
          "k": eis["k"] if eis else None,
          "dE": eis["dE"] if eis else None,
          "reuse": (reuseA[accA.index(eis)] if eis else None),
          "pass": bool(eis and eis["k"] == 143 and eis["dE"] < 0)}
    report["K1"] = K1

    def identifiable(a, fams):
        for f in fams:
            for i in range(f.N):
                shared = 0
                ok = True
                for c in range(C):
                    v = f.rows[i][c]
                    if v is None or a["rh"][c] < 0:
                        continue
                    shared += 1
                    if v % a["l"] != a["rh"][c]:
                        ok = False
                        break
                if ok and shared >= MIN_SHARED:
                    return "row:" + f.labels[i]
        for pi in p02.PI_ORDER:
            fn_ = p02.PI_FUN[pi]
            if all(a["rh"][c] < 0 or a["rh"][c] == fn_(PRIMES[c]) % a["l"]
                   for c in range(C)):
                return "pattern:" + pi
        return None

    K2 = []
    for ai, a in enumerate(accB):
        ident = identifiable(a, SB)
        K2.append({"id": a["id"], "l": a["l"], "k": a["k"],
                   "dE": round(a["dE"], 1), "reuse": reuseB[ai],
                   "identifiable": ident})
    K3 = [e for e in K2 if e["reuse"] >= 1 and e["identifiable"] is None]
    report["K2"] = {"count": len(K2), "anchors": K2[:50]}
    report["K3"] = {"count": len(K3), "list": K3}
    report["K4"] = {"SA_depths": sorted(set(depthA)) if depthA else [],
                    "SB_depths": sorted(set(depthB)) if depthB else [],
                    "SA_depth_hist": {str(d): depthA.count(d)
                                      for d in set(depthA)},
                    "SB_depth_hist": {str(d): depthB.count(d)
                                      for d in set(depthB)}}
    print("K1:", K1)
    print("K2 count:", len(K2), "K3 count:", len(K3), "K4:", report["K4"])

    # ---------- ladder ----------
    if not H_pass:
        verdict = "VOID"
    elif not n1["pass"]:
        verdict = "FATAL"
    elif not n3["pass"]:
        verdict = "BLIND"
    elif not K1["pass"]:
        verdict = "MOVE-DEAD"
    elif len(K3) == 0 and len(cands_B) < 500:
        verdict = "UNDERPOWERED"
    elif len(K3) == 0:
        verdict = "NO-MINT"
    else:
        verdict = "MINTED"
    report["verdict"] = verdict
    report["runtime_s"] = round(time.time() - t0, 1)

    with open(os.path.join(args.out, "anchors_SB.json"), "w") as fh:
        json.dump(K2, fh, indent=1)
    with open(os.path.join(args.out, "summary.json"), "w") as fh:
        json.dump(report, fh, indent=1, default=str)
    print(f"\n== VERDICT: {verdict} == ({report['runtime_s']}s)")


if __name__ == "__main__":
    main()
