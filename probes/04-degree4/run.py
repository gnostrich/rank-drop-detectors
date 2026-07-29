#!/usr/bin/env python3
"""Probe 04 — degree 4, incomplete correspondence. Implements PREREG.md.

Usage:
    python3 run.py --data data --p03data ../03-transport/data --out output

Inherits probe 03's pipeline by module import (Block, column-norm
substitution, scans, pattern entries) and probe 02's energy through it.
New here: genus-2 point counting over F_p and F_p^2 (H1) and the H4
product-decomposability sieve. No minting move exists in this file.
"""

import argparse
import csv
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
SEED = 20260728
EPS = 6
MIN_SHARED = 10
DEG4_NORMS = [4 * p for p in PRIMES]      # alphabet bound floor(4*sqrt(p))


def load_module(name, relpath):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, relpath))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


p03 = load_module("p03run", "../03-transport/run.py")
p02 = p03.p02
p01 = p03.p01


# ---------------------------------------------------- genus-2 point counts

def poly_g(f, h, p):
    """g = 4f + h^2 mod p, coefficient list (low to high), trimmed."""
    g = [0] * 13
    for i, cf in enumerate(f):
        g[i] = (g[i] + 4 * cf) % p
    for i, a in enumerate(h):
        for j, b in enumerate(h):
            g[i + j] = (g[i + j] + a * b) % p
    while len(g) > 1 and g[-1] % p == 0:
        g.pop()
    return [c % p for c in g]


def count_Fp(f, h, p):
    """#C(F_p) for y^2 + h y = f, p odd, good curve reduction assumed."""
    g = poly_g(f, h, p)
    qr = np.zeros(p, dtype=np.int8)
    z = np.arange(1, p, dtype=np.int64)
    qr[(z * z) % p] = 1
    x = np.arange(p, dtype=np.int64)
    v = np.zeros(p, dtype=np.int64)
    for cf in reversed(g):
        v = (v * x + cf) % p
    chi = np.where(v == 0, 0, np.where(qr[v] == 1, 1, -1))
    naff = p + int(chi.sum())
    deg = len(g) - 1
    lc = g[-1]
    if deg == 6:
        ninf = 1 + (0 if lc == 0 else (1 if qr[lc] else -1))
    elif deg == 5:
        ninf = 1
    else:                      # degenerate degree at this p: cannot count
        return None
    return naff + ninf


def count_Fp2(f, h, p):
    """#C(F_{p^2}); field = F_p[s]/(s^2 - r), r least non-residue."""
    qr = np.zeros(p, dtype=np.int8)
    z = np.arange(1, p, dtype=np.int64)
    qr[(z * z) % p] = 1
    r = next(t for t in range(2, p) if not qr[t])
    g = poly_g(f, h, p)
    # all q = p^2 elements as (u, v)
    u = np.repeat(np.arange(p, dtype=np.int64), p)
    v = np.tile(np.arange(p, dtype=np.int64), p)

    def mul(a0, a1, b0, b1):
        return (a0 * b0 + r * a1 * b1) % p, (a0 * b1 + a1 * b0) % p

    # squares table of F_q indexed by u*p+v
    s0, s1 = mul(u, v, u, v)
    sq = np.zeros(p * p, dtype=np.int8)
    sq[s0 * p + s1] = 1
    nz = ~((u == 0) & (v == 0))
    sq0 = sq.copy()
    # evaluate g at every x = (u, v)
    e0 = np.zeros(p * p, dtype=np.int64)
    e1 = np.zeros(p * p, dtype=np.int64)
    for cf in reversed(g):
        e0, e1 = mul(e0, e1, u, v)
        e0 = (e0 + cf) % p
    idx = e0 * p + e1
    zero = (e0 == 0) & (e1 == 0)
    chi = np.where(zero, 0, np.where(sq0[idx] == 1, 1, -1))
    q = p * p
    naff = q + int(chi.sum())
    deg = len(g) - 1
    lc = g[-1] % p             # leading coeff is in F_p
    if deg == 6:
        ninf = 1 + (0 if lc == 0 else (1 if sq0[lc * p] else -1))
    elif deg == 5:
        ninf = 1
    else:
        return None
    return naff + ninf


def euler_c1_c2(f, h, p):
    n1 = count_Fp(f, h, p)
    n2 = count_Fp2(f, h, p)
    if n1 is None or n2 is None:
        return None
    t1 = p + 1 - n1
    t2 = p * p + 1 - n2
    return -t1, (t1 * t1 - t2) // 2


# ---------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--p03data", default=os.path.join(HERE, "..",
                                                      "03-transport", "data"))
    ap.add_argument("--out", default="output")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    report = {"seed": SEED}

    # ---------- load ----------
    print("== load ==")
    f5_rows = p03.read_csv(os.path.join(args.data, "f5_g2c.csv"))
    f6_rows = p03.read_csv(os.path.join(args.data, "f6_class_euler.csv"))
    f7_rows = p03.read_csv(os.path.join(args.data, "f7_smf.csv"))
    class_euler = {}
    for r in f6_rows:
        ef = json.loads(r["euler_factors"])
        class_euler[(int(r["cond"]), r["class"])] = ef

    def cond_markers(N):
        return [p_ for p_ in PRIMES if N % p_ == 0]

    def trace_row(N, ef):
        row = []
        for i, p_ in enumerate(PRIMES):
            if N % p_ == 0 or ef is None or ef[i] is None or len(ef[i]) < 5:
                row.append(None)
            else:
                row.append(-ef[i][1])
        return row

    f5_vals, f5_meta = [], []
    for r in sorted(f5_rows, key=lambda r: (int(r["cond"]), r["label"])):
        N = int(r["cond"])
        ef = class_euler.get((N, r["class"]))
        f5_vals.append(trace_row(N, ef))
        f5_meta.append({"label": r["label"], "class": f"{N}.{r['class']}",
                        "cond": N, "abs_disc": int(r["abs_disc"]),
                        "eqn": r["eqn"], "gl2": r["is_gl2_type"],
                        "end_alg": r["end_alg"],
                        "geom_end_alg": r["geom_end_alg"],
                        "sealed": int(r["sealed"])})
    F5 = p03.Block("F5", DEG4_NORMS, f5_vals, f5_meta)

    f6_vals, f6_meta = [], []
    for (N, cls), ef in sorted(class_euler.items()):
        f6_vals.append(trace_row(N, ef))
        f6_meta.append({"class": f"{N}.{cls}", "cond": N,
                        "sealed": int(N > 9000)})
    F6 = p03.Block("F6", DEG4_NORMS, f6_vals, f6_meta)

    f7_vals, f7_meta = [], []
    for r in f7_rows:
        lam = json.loads(r["trace_lambda_p"]) if r["trace_lambda_p"] else []
        row = []
        for i, p_ in enumerate(PRIMES):
            v = lam[i] if i < len(lam) else None
            row.append(int(v) if isinstance(v, int) else None)
        f7_vals.append(row)
        f7_meta.append({"label": r["label"], "level": int(r["level"]),
                        "family": r["family"], "dim": int(r["dim"]),
                        "type": r["aut_rep_type"],
                        "related": r["related_objects"], "sealed": 0})
    F7 = p03.Block("F7", DEG4_NORMS, f7_vals, f7_meta)

    # degree-2 reference blocks from probe 03's committed snapshots
    f1_rows = p03.read_csv(os.path.join(args.p03data, "f1_ec_q.csv"))
    import re as _re
    f1_classes = {}
    for r in f1_rows:
        m = _re.match(r"^(\d+)\.([a-z]+)(\d+)$", r["lmfdb_label"])
        if m and int(m.group(3)) == 1:
            f1_classes[f"{m.group(1)}.{m.group(2)}"] = r
    f1_vals, f1_meta = [], []
    for cls in sorted(f1_classes, key=lambda s: (int(s.split(".")[0]),
                                                 s.split(".")[1])):
        r = f1_classes[cls]
        N = int(r["conductor"])
        ain = json.loads(r["ainvs"])
        f1_vals.append([None if N % p_ == 0 else p01.ap_point_count(ain, p_)
                        for p_ in PRIMES])
        f1_meta.append({"cls": cls, "N": N, "sealed": 0})
    F1 = p03.Block("F1", PRIMES, f1_vals, f1_meta)
    f2_rows = p03.read_csv(os.path.join(args.p03data, "f2_newforms.csv"))
    n_f2 = len(f2_rows)

    def split_sealed(block):
        kv, km, ns = [], [], 0
        for v, m in zip(block.rows, block.meta):
            if m.get("sealed"):
                ns += 1
            else:
                kv.append(v)
                km.append(m)
        return p03.Block(block.name, block.norms, kv, km), ns

    blocks, sealed = {}, {}
    for name, b in (("F5", F5), ("F6", F6), ("F7", F7), ("F1", F1)):
        blocks[name], sealed[name] = split_sealed(b)
        blocks[name].build_table()
    N_total = sum(b.N for b in blocks.values())
    log2N = math.log2(N_total)
    report["rows_open"] = {k: b.N for k, b in blocks.items()}
    report["rows_sealed"] = sealed
    report["F2_reference_rows"] = n_f2
    report["F7_weight20_rows_found"] = F7.N
    print("open:", report["rows_open"], "sealed:", sealed,
          f"({time.time()-t0:.0f}s)")

    # ---------- H1 recomputation ----------
    print("== H1 ==")
    rng = np.random.default_rng(SEED)
    idx = rng.choice(blocks["F5"].N, size=200, replace=False)
    mism = checked = fe_bad = 0
    for i in idx.tolist():
        m = blocks["F5"].meta[i]
        eqn = json.loads(m["eqn"])
        if isinstance(eqn, str):
            eqn = json.loads(eqn)
        f, h = eqn
        ef = class_euler[(m["cond"], m["class"].split(".")[1])]
        for k, p_ in enumerate(PRIMES):
            if p_ == 2 or m["cond"] % p_ == 0 or m["abs_disc"] % p_ == 0:
                continue
            got = euler_c1_c2(f, h, p_)
            if got is None:
                continue
            c1, c2 = got
            st = ef[k]
            checked += 1
            if int(st[1]) != c1 or int(st[2]) != c2:
                mism += 1
            if int(st[3]) != p_ * int(st[1]) or int(st[4]) != p_ * p_:
                fe_bad += 1
    h1 = {"rows": 200, "columns_checked": checked, "mismatches": mism,
          "functional_eq_violations": fe_bad,
          "pass": mism == 0 and fe_bad == 0}
    report["H1"] = h1
    print(h1, f"({time.time()-t0:.0f}s)")

    # ---------- H4 sieve (before any residual count) ----------
    print("== H4 sieve ==")
    f1_by_cond = defaultdict(list)
    for i, m in enumerate(blocks["F1"].meta):
        f1_by_cond[m["N"]].append(i)

    def sieve_row(vals, N, depth=C):
        """Return matching (E1, E2) class pair or None; depth = #columns."""
        for N1 in range(11, int(math.isqrt(N)) + 1):
            if N % N1 or N1 not in f1_by_cond:
                continue
            N2 = N // N1
            if N2 not in f1_by_cond:
                continue
            for i1 in f1_by_cond[N1]:
                a1 = blocks["F1"].rows[i1]
                for i2 in f1_by_cond[N2]:
                    a2 = blocks["F1"].rows[i2]
                    shared = ok = 0
                    for c in range(depth):
                        if vals[c] is None or a1[c] is None or a2[c] is None:
                            continue
                        shared += 1
                        if vals[c] == a1[c] + a2[c]:
                            ok += 1
                    if shared >= MIN_SHARED and ok == shared:
                        return (blocks["F1"].meta[i1]["cls"],
                                blocks["F1"].meta[i2]["cls"])
        return None

    sieve_hits = {}
    for j, m in enumerate(blocks["F6"].meta):
        hit = sieve_row(blocks["F6"].rows[j], m["cond"])
        if hit:
            sieve_hits[m["class"]] = hit
    gl2_classes = {m["class"] for m in blocks["F5"].meta
                   if str(m["gl2"]).lower() in ("true", "t", "1")}
    product_end = {m["class"] for m in blocks["F5"].meta
                   if "x" in str(m["end_alg"]) or "M_2" in str(m["end_alg"])}
    hits_set = set(sieve_hits)
    h4 = {"decomposable": len(sieve_hits),
          "gl2_type_classes": len(gl2_classes),
          "product_endalg_classes": len(product_end),
          "hits_with_product_endalg": len(hits_set & product_end),
          "product_endalg_not_hit": len(product_end - hits_set),
          "implemented": True}
    report["H4"] = h4
    print(h4, f"({time.time()-t0:.0f}s)")

    # ---------- codebook + H2 ----------
    print("== codebook + H2 ==")
    patterns = {n: (p03.pattern_entries(b, log2N) if b.N else [])
                for n, b in blocks.items() if n != "F1"}
    cross = {}
    for a, bnm in (("F6", "F5"),):      # class -> curve rows (bookkeeping)
        ms = p03.cross_dup_scan(blocks[a], blocks[bnm])
        kept = [(i, j, s, p03.dup_delta(blocks[a], blocks[bnm], i, j, log2N))
                for (i, j, s) in ms]
        cross[(a, bnm)] = [t for t in kept if t[3] < 0]
    for a, bnm in (("F5", "F7"), ("F6", "F7")):
        if blocks["F7"].N:
            ms = p03.cross_dup_scan(blocks[a], blocks[bnm])
            kept = [(i, j, s,
                     p03.dup_delta(blocks[a], blocks[bnm], i, j, log2N))
                    for (i, j, s) in ms]
            cross[(a, bnm)] = [t for t in kept if t[3] < 0]
        else:
            cross[(a, bnm)] = []
    report["codebook"] = {
        "cross_dups": {f"{a}x{b}": len(v) for (a, b), v in cross.items()},
        "patterns": {k: len(v) for k, v in patterns.items()}}
    print(report["codebook"])

    h2 = {"streams": [], "pass": True}
    for name in ("F5", "F6", "F7", "F1"):
        b = blocks[name]
        if b.N == 0:
            h2["streams"].append({"block": name, "book": "empty",
                                  "identical": True, "rel": 0.0,
                                  "note": "0 rows"})
            continue
        books = [("empty", [])]
        pats = [p02.Entry("PATTERN", l=e["l"], members=e["members"],
                          pi=e["pi"]) for e in patterns.get(name, [])]
        if pats:
            books.append((f"patterns({len(pats)})", pats))
        for bname, book in books:
            with p03.p02_norms(b.norms):
                codec = p02.StreamCodec(b.table, book, EPS)
                blob = codec.encode()
                ok = bool((codec.decode(blob) == b.table.SYM).all())
                if book:
                    owned = p02.book_delta_E_owned(b.table, book, EPS)
                    analytic = b.table.L0_data + owned["cells_delta_bits"]
                else:
                    analytic = b.table.L0_data
            rel = abs(8 * len(blob) - analytic) / analytic
            h2["streams"].append({"block": name, "book": bname,
                                  "identical": ok, "rel": round(rel, 7)})
            h2["pass"] &= ok and rel <= 0.005
    report["H2"] = h2
    print("H2 pass:", h2["pass"], f"({time.time()-t0:.0f}s)")

    # ---------- H3 known-link recovery ----------
    curated = set()
    for j, m in enumerate(blocks["F7"].meta):
        rel = json.loads(m["related"]) if m["related"] else []
        for u in rel:
            if "Genus2Curve/Q/" in u:
                parts = u.split("Genus2Curve/Q/")[1].split("/")
                curated.add((f"{parts[0]}.{parts[1]}", m["label"]))
    blind = set()
    for (a, bnm) in (("F5", "F7"), ("F6", "F7")):
        for (i, j, s, dE) in cross[(a, bnm)]:
            cls = (blocks[a].meta[i]["class"] if a in ("F5", "F6")
                   else blocks[a].meta[i]["cls"])
            blind.add((cls, blocks["F7"].meta[j]["label"]))
    rec = (len(blind & curated) / len(curated)) if curated else None
    h3 = {"curated_subset_size": len(curated),
          "blind_recovered": len(blind & curated), "recall": rec,
          "smoke_test": len(curated) < 50,
          "pass": True if not curated else (rec >= 0.95)}
    report["H3"] = h3
    print("H3:", h3)

    H_pass = h1["pass"] and h2["pass"] and h3["pass"] and h4["implemented"]

    # ---------- N1 shuffled null ----------
    n1_total = 0
    fam_ids = {"F5": 7, "F6": 8, "F7": 9}
    shuf = {}
    for name, fid in fam_ids.items():
        b = blocks[name]
        rngs = np.random.default_rng(SEED + fid)
        if b.N == 0:
            shuf[name] = (np.zeros((0, C), dtype=np.int64),
                          np.zeros((0, C), dtype=bool))
            continue
        V, G = p03.val_arr(b), p03.good_mask(b)
        for c in range(C):
            perm = rngs.permutation(b.N)
            V[:, c], G[:, c] = V[perm, c], G[perm, c]
        shuf[name] = (V, G)
    for (a, bnm) in (("F6", "F5"), ("F5", "F7"), ("F6", "F7")):
        VP, GP = shuf[a]
        VC, GC = shuf[bnm]
        for j in range(blocks[bnm].N):
            both = GP & GC[j]
            neq = (VP != VC[j]) & both
            ok = (~neq.any(axis=1)) & (both.sum(axis=1) >= MIN_SHARED)
            n1_total += int(ok.sum())
    report["N1"] = {"cross_matches_on_shuffled": n1_total,
                    "pass": n1_total == 0}
    print("N1:", report["N1"])

    # ---------- N2 foil power ----------
    rngf = np.random.default_rng(SEED)
    detected = total = 0
    src = blocks["F6"]
    B = [math.isqrt(16 * p_) for p_ in PRIMES]

    def dup_match(rowA, rowB):
        shared = 0
        for c in range(C):
            if rowA[c] is None or rowB[c] is None:
                continue
            shared += 1
            if rowA[c] != rowB[c]:
                return False
        return shared >= MIN_SHARED

    for ftype in ("perturb", "wrong-level", "wrong-weight"):
        for _ in range(1000):
            i = int(rngf.integers(0, src.N))
            row = list(src.rows[i])
            if ftype == "perturb":
                goods = [c for c in range(C) if row[c] is not None]
                if not goods:
                    continue
                c = goods[int(rngf.integers(0, len(goods)))]
                v = row[c] + 1
                row[c] = -B[c] if v > B[c] else v
            elif ftype == "wrong-level":
                j = (i + 1) % src.N
                while src.meta[j]["cond"] == src.meta[i]["cond"]:
                    j = (j + 1) % src.N
                row = list(src.rows[j])
            else:
                row = [None if v is None else
                       (min(max(2 * v, -B[c]), B[c])) for c, v in
                       enumerate(row)]
            total += 1
            if not dup_match(row, src.rows[i]):
                detected += 1
    n2 = {"foils": total, "detected": detected,
          "power": detected / total if total else 0.0,
          "pass": total > 0 and detected / total >= 0.95}
    report["N2"] = n2
    print("N2:", n2)

    # ---------- N3 depth curve ----------
    print("== N3 depth curve ==")
    depth_curve = {}
    for d in (10, 25):
        VP, GP = p03.val_arr(blocks["F6"]), p03.good_mask(blocks["F6"])
        VC, GC = p03.val_arr(blocks["F5"]), p03.good_mask(blocks["F5"])
        m56 = 0
        for j in range(blocks["F5"].N):
            both = GP[:, :d] & GC[j, :d]
            neq = (VP[:, :d] != VC[j, :d]) & both
            ok = (~neq.any(axis=1)) & (both.sum(axis=1) >= min(MIN_SHARED, d))
            m56 += int(ok.sum())
        m7 = 0
        if blocks["F7"].N:
            V7, G7 = p03.val_arr(blocks["F7"]), p03.good_mask(blocks["F7"])
            for j in range(blocks["F7"].N):
                both = GP[:, :d] & G7[j, :d]
                neq = (VP[:, :d] != V7[j, :d]) & both
                ok = (~neq.any(axis=1)) & (both.sum(axis=1)
                                           >= min(MIN_SHARED, d))
                m7 += int(ok.sum())
        sieve_d = 0
        for j, m in enumerate(blocks["F6"].meta):
            if sieve_row(blocks["F6"].rows[j], m["cond"], depth=d):
                sieve_d += 1
        depth_curve[str(d)] = {"F6xF5": m56, "F6xF7_and_F5xF7": m7,
                               "sieve_decomposable": sieve_d}
    depth_curve["50"] = depth_curve["all"] = dict(
        depth_curve["25"], note="only 25 primes available; 50/all = 25")
    report["N3_depth_curve"] = depth_curve
    print(depth_curve)

    # ---------- K1 residual census ----------
    label_links = set()
    for (i, j, s, dE) in cross[("F6", "F5")]:
        if blocks["F5"].meta[j]["class"] == blocks["F6"].meta[i]["class"]:
            label_links.add((blocks["F6"].meta[i]["class"],
                             blocks["F5"].meta[j]["label"]))
    all_cross = set()
    for (a, bnm), v in cross.items():
        for (i, j, s, dE) in v:
            la = (blocks[a].meta[i].get("class")
                  or blocks[a].meta[i].get("label"))
            lb = (blocks[bnm].meta[j].get("label")
                  or blocks[bnm].meta[j].get("class"))
            all_cross.add((a, bnm, la, lb))
    residuals = []
    for (a, bnm, la, lb) in sorted(all_cross):
        if a == "F6" and bnm == "F5":
            same_class = lb.startswith(la.split(".")[0] + ".") and \
                la.split(".")[1] == lb.split(".")[1]
            if same_class:
                continue                      # curated by label (bookkeeping)
            if la in sieve_hits:
                continue                      # decomposable
            residuals.append((a, bnm, la, lb))
        else:
            if (la, lb) in curated:
                continue
            if la in sieve_hits:
                continue
            residuals.append((a, bnm, la, lb))
    K1 = {"cross_matches": len(all_cross),
          "label_bookkeeping": len(label_links),
          "residual_count": len(residuals),
          "residuals": residuals[:100]}
    report["K1"] = K1
    print("K1:", {k: v for k, v in K1.items() if k != "residuals"})

    # ---------- K2 ----------
    comparable = (blocks["F5"].N + blocks["F6"].N) * blocks["F7"].N
    K2 = {"rows": report["rows_open"], "primes_per_row": C,
          "comparable_pairs_deg4_x_F7": comparable,
          "curated_links_for_H3": len(curated),
          "case": ("SAMPLE-SIZE" if comparable < 10000 else "ADEQUATE")}
    report["K2"] = K2
    print("K2:", K2)

    # ---------- K3 ----------
    k3 = {}
    for name in ("F5", "F6", "F7"):
        b = blocks[name]
        if b.N == 0:
            k3[name] = {"L0": 0, "bits_saved": 0, "fraction": None,
                        "billing": "empty family"}
            continue
        saved = -sum(e["dE"] for e in patterns.get(name, []))
        billing = "structure only"
        if name == "F5":
            dup_saved = -sum(dE for (_i, _j, _s, dE) in cross[("F6", "F5")])
            saved += dup_saved
            billing = ("duplication artifact: F5 rows share class "
                       "L-functions with F6; this fraction is a copy, "
                       "NOT comparable to probe 02's 1.67%")
        k3[name] = {"L0": b.table.L0, "bits_saved": saved,
                    "fraction": saved / b.table.L0, "billing": billing}
    report["K3"] = k3
    print("K3:", {k: (round(v["fraction"], 4) if v["fraction"] else v["fraction"])
                  for k, v in k3.items()})

    # ---------- ladder ----------
    if not H_pass:
        verdict = "VOID"
    elif not n2["pass"]:
        verdict = "BLIND"
    elif K1["residual_count"] == 0 and comparable < 10000:
        verdict = "UNDERPOWERED"
    elif K1["residual_count"] == 0:
        verdict = "NULL"
    else:
        verdict = "RESIDUALS"
    report["verdict"] = verdict
    report["runtime_s"] = round(time.time() - t0, 1)

    with open(os.path.join(args.out, "sieve.csv"), "w") as fh:
        fh.write("class,E1,E2\n")
        for cls, (e1, e2) in sorted(sieve_hits.items()):
            fh.write(f"{cls},{e1},{e2}\n")
    with open(os.path.join(args.out, "summary.json"), "w") as fh:
        json.dump(report, fh, indent=1, default=str)
    print(f"\n== VERDICT: {verdict} == ({report['runtime_s']}s)")


if __name__ == "__main__":
    main()
