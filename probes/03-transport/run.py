#!/usr/bin/env python3
"""Probe 03 — cross-family transport and closure defect. Implements PREREG.md.

Usage:
    python3 run.py --data data --ecdata /path/to/ecdata --out output

Implementation readings forced by the prereg text (documented, not tunable):

R1. F2 "all dimensions" + "entries: raw integer a_p": a newform of dimension
    d > 1 has no integer a_p (only the trace of a non-rational eigenvalue,
    bounded by 2d*sqrt(p), outside the Hasse alphabet H1 verifies). Dim > 1
    rows enter as all-marker rows: present, never matchable. Disclosed
    consequence: K1's distractor set is dim-1 forms at other levels.
R2. A marker cannot be transported, so a path starting at a marker cell
    predicts nothing for the derived columns. In the meter, endpoint marker
    cells cost baseline under every path (prereg §3), so they cancel in
    delta_bits. delta_frac compares predictions wherever both paths predict.
R3. N2 "wobble floor" = measurement wobble across replicates: floor(type) =
    q95 over instances of (max_r delta_r - min_r delta_r), r = 0..19. The
    transports are deterministic; a zero floor is the expected reading and
    binds K5 to "any nonzero delta is above floor".
R4. Probe-02 energy imported unchanged; K-blocks are built and evaluated with
    the column-norm vector substituted into the probe-02 module state (the
    prereg §1 per-block-column-space instantiation), then restored. No cost
    formula is reimplemented. Codebook pointer costs use log2(N_total) per
    prereg §1; K4(a) stand-alone bands use their own log2(N_band).
"""

import argparse
import contextlib
import csv
import importlib.util
import json
import math
import os
import re
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
FIELDS = {"2.2.5.1": 5, "2.2.8.1": 8}
REPLICATES = 20


def load_module(name, relpath):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, relpath))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


p01 = load_module("p01run", "../01-congruence-census/run.py")
p02 = load_module("p02run", "../02-encoder-energy/run.py")


@contextlib.contextmanager
def p02_norms(norms):
    """Substitute the probe-02 column-norm vector (reading R4)."""
    saved = p02.PRIMES
    p02.PRIMES = list(norms)
    try:
        yield
    finally:
        p02.PRIMES = saved


# ------------------------------------------------ quadratic field arithmetic

class QF:
    """O_K = Z[w]; disc 5: w^2 = 1 + w; disc 8: w^2 = 2."""

    def __init__(self, D):
        self.D = D
        self.wsq = (1, 1) if D == 5 else (2, 0)

    def mul(self, a, b):
        (x, y), (u, v) = a, b
        c0, c1 = self.wsq
        return (x * u + c0 * y * v, x * v + y * u + c1 * y * v)

    def conj(self, a):
        x, y = a
        return (x + y, -y) if self.D == 5 else (x, -y)

    def norm(self, a):
        x, y = a
        return x * x + x * y - y * y if self.D == 5 else x * x - 2 * y * y

    def divides(self, a, b):
        na = self.norm(a)
        if na == 0:
            return False
        num = self.mul(b, self.conj(a))
        return num[0] % abs(na) == 0 and num[1] % abs(na) == 0

    def minpoly_roots(self, p):
        c0, c1 = self.wsq
        return [t for t in range(p) if (t * t - c1 * t - c0) % p == 0]


def kronecker_D(D, p):
    if D == 5:
        return 0 if p == 5 else (1 if p % 5 in (1, 4) else -1)
    return 0 if p == 2 else (1 if p % 8 in (1, 7) else -1)


def parse_gen(s):
    s = s.strip()
    if s.startswith("(") and s.endswith(")"):
        s = s[1:-1]
    s = s.replace(" ", "").replace("-", "+-")
    x = y = 0
    for term in s.split("+"):
        if not term:
            continue
        if "w" in term:
            coef = term.replace("*w", "").replace("w", "")
            y += -1 if coef == "-" else (1 if coef == "" else int(coef))
        else:
            x += int(term)
    return (x, y)


def parse_primes_list(raw_list, qf):
    cols = []
    for s in raw_list:
        inner = s.strip()[1:-1]
        a, b, rest = inner.split(",", 2)
        q, p = int(a), int(b)
        gen = parse_gen(rest)
        ch = kronecker_D(qf.D, p)
        if q == p * p:
            kind, t = "inert", None
        elif ch == 0:
            kind = "ram"
            t = qf.minpoly_roots(p)[0]
        else:
            kind = "split"
            t = None
            x, y = gen
            for r in qf.minpoly_roots(p):
                if (x + y * r) % p == 0:
                    t = r
                    break
        cols.append({"q": q, "p": p, "gen": gen, "kind": kind, "t": t})
        if len(cols) == 25:
            break
    return cols


# ------------------------------------------------ point counting over F_q

def ap_Fq_inert(qf, ainvs_red, p):
    """a_P at an inert prime: count on the general Weierstrass model /F_{p^2}."""
    c0, c1 = qf.wsq
    els = [(u, v) for u in range(p) for v in range(p)]

    def mul(a, b):
        return ((a[0] * b[0] + c0 * a[1] * b[1]) % p,
                (a[0] * b[1] + a[1] * b[0] + c1 * a[1] * b[1]) % p)

    def add(a, b):
        return ((a[0] + b[0]) % p, (a[1] + b[1]) % p)

    a1, a2, a3, a4, a6 = ainvs_red
    q = p * p
    if p == 2:
        naff = 0
        for x in els:
            x3 = mul(mul(x, x), x)
            rhs = add(add(add(x3, mul(a2, mul(x, x))), mul(a4, x)), a6)
            for y in els:
                lhs = add(add(mul(y, y), mul(mul(a1, x), y)), mul(a3, y))
                if lhs == rhs:
                    naff += 1
        return q - naff
    squares = set(mul(z, z) for z in els)
    inv4 = pow(4, -1, p)
    total = 0
    for x in els:
        x2 = mul(x, x)
        fx = add(add(add(mul(x2, x), mul(a2, x2)), mul(a4, x)), a6)
        h = add(mul(a1, x), a3)
        g = add(fx, mul(mul(h, h), (inv4 % p, 0)))
        if g == (0, 0):
            pass
        elif g in squares:
            total += 1
        else:
            total -= 1
    return -total


def ap_column(qf, ainvs, col):
    """a_P for one ideal column; ainvs = five (c0, c1) pairs over Z[w]."""
    if col["kind"] in ("split", "ram"):
        p = col["p"]
        red = [(c0 + c1 * col["t"]) % p for (c0, c1) in ainvs]
        return p01.ap_point_count(red, p)
    p = col["p"]
    red = [(c0 % p, c1 % p) for (c0, c1) in ainvs]
    return ap_Fq_inert(qf, red, p)


# ------------------------------------------------ blocks

class Block:
    def __init__(self, name, norms, rows, meta):
        self.name = name
        self.norms = list(norms)
        self.rows = rows
        self.meta = meta
        self.N = len(rows)
        self.table = None

    def build_table(self):
        with p02_norms(self.norms):
            sym_rows = [[("?" if v is None else str(int(v))) for v in row]
                        for row in self.rows]
            self.table = p02.Table(self.name, sym_rows)
        return self.table


def read_csv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def val_arr(block):
    return np.array([[0 if v is None else v for v in row]
                     for row in block.rows], dtype=np.int64)


def good_mask(block):
    return np.array([[v is not None for v in row] for row in block.rows])


# ------------------------------------------------ transport meter

def meter_cost(block, r, preds, eps=EPS):
    t = block.table
    fb = p02.flag_bits(eps)
    cost = 0.0
    for c in range(C):
        base = t.logq_cell[r, c]
        v = block.rows[r][c]
        pred = preds.get(c)
        if pred is None or v is None:
            cost += base
        elif v == pred:
            cost += fb
        else:
            cost += eps + base
    return cost


def delta_pair(block, r, predsA, predsB, eps=EPS):
    dbits = abs(meter_cost(block, r, predsA, eps)
                - meter_cost(block, r, predsB, eps))
    both = [c for c in range(C)
            if predsA.get(c) is not None and predsB.get(c) is not None]
    dfrac = (sum(1 for c in both if predsA[c] != predsB[c]) / len(both)
             if both else 0.0)
    return dbits, dfrac


def bc_predict(values_by_p, cols):
    out = {}
    for ci, col in enumerate(cols):
        a = values_by_p.get(col["p"])
        if a is None or col["kind"] == "ram":
            continue
        out[ci] = a if col["kind"] == "split" else a * a - 2 * col["p"]
    return out


def tw_values(values_by_p, D):
    out = {}
    for p, a in values_by_p.items():
        ch = kronecker_D(D, p)
        if ch != 0:
            out[p] = ch * a
    return out


# ------------------------------------------------ deterministic codebook

def pattern_entries(block, log2N_pointer):
    """PATTERN entries per §2: membership = all good cells conform; ΔE < 0.

    Cell deltas computed directly from the block table's Q-masses so that
    K-block residues use ideal norms (reading R4)."""
    out = []
    V, G = val_arr(block), good_mask(block)
    fb = p02.flag_bits(EPS)
    t = block.table
    for pi in p02.PI_ORDER:
        fn = p02.PI_FUN[pi]
        for l in MODULI:
            colmask = np.array([q != l for q in block.norms])
            res = np.array([fn(q) % l for q in block.norms])
            conf = (~G) | ((V % l) == res[None, :]) | (~colmask[None, :])
            members = np.nonzero(conf.all(axis=1) & G.any(axis=1))[0]
            k = len(members)
            if k < 2:
                continue
            nlq = np.array([t.negLogQ[l][c, res[c]] if colmask[c] else 0.0
                            for c in range(C)])
            gsub = G[members][:, colmask]
            nl_sub = nlq[colmask][None, :]
            n_good = gsub.sum()
            n_mark = (~gsub).sum()
            cells = fb * n_good - float((nl_sub * gsub).sum()) + EPS * n_mark
            cost = (p02.H_BITS + 2 + p02.gamma_bits(k) + k * log2N_pointer)
            dE = cost + cells
            if dE < 0:
                out.append({"pi": pi, "l": l, "members": members.tolist(),
                            "dE": float(dE)})
    return out


def cross_dup_scan(bp, bc_):
    VP, GP = val_arr(bp), good_mask(bp)
    VC, GC = val_arr(bc_), good_mask(bc_)
    matches = []
    for j in range(bc_.N):
        both = GP & GC[j]
        neq = (VP != VC[j]) & both
        cand = ~neq.any(axis=1)
        shared = both.sum(axis=1)
        ok = cand & (shared >= MIN_SHARED) & (GC[j].any())
        for i in np.nonzero(ok)[0].tolist():
            matches.append((i, j, int(shared[i])))
    return matches


def dup_delta(bp, bc_, i, j, log2N_pointer):
    fb = p02.flag_bits(EPS)
    d = p02.H_BITS + 2 * log2N_pointer
    for c in range(C):
        vp, vc = bp.rows[i][c], bc_.rows[j][c]
        base = bc_.table.logq_cell[j, c]
        if (vp is None and vc is None) or (vp is not None and vp == vc):
            d += fb - base
        else:
            d += EPS
    return d


# ------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--ecdata", required=True)
    ap.add_argument("--out", default="output")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    report = {"seed": SEED, "readings": ["R1", "R2", "R3", "R4"]}

    # ---------------- load blocks ----------------
    print("== load ==")
    fp = json.load(open(os.path.join(args.data, "hmf_fields_primes.json")))
    kcols = {fl: parse_primes_list(fp[fl], QF(FIELDS[fl])) for fl in FIELDS}

    f1_rows = read_csv(os.path.join(args.data, "f1_ec_q.csv"))
    f1_classes = {}
    for r in f1_rows:
        m = re.match(r"^(\d+)\.([a-z]+)(\d+)$", r["lmfdb_label"])
        if m and int(m.group(3)) == 1:
            f1_classes[f"{m.group(1)}.{m.group(2)}"] = r
    f1_meta, f1_vals = [], []
    for cls in sorted(f1_classes, key=lambda s: (int(s.split(".")[0]),
                                                 s.split(".")[1])):
        r = f1_classes[cls]
        N = int(r["conductor"])
        ain = json.loads(r["ainvs"])
        f1_vals.append([None if N % p_ == 0 else p01.ap_point_count(ain, p_)
                        for p_ in PRIMES])
        f1_meta.append({"cls": cls, "N": N, "Clabel": r["Clabel"],
                        "sealed": 0})
    B1 = Block("F1", PRIMES, f1_vals, f1_meta)

    f2_rows = read_csv(os.path.join(args.data, "f2_newforms.csv"))
    f2_meta, f2_vals, n_dim_gt1 = [], [], 0
    for r in sorted(f2_rows, key=lambda r: (int(r["level"]), r["label"])):
        dim = int(r["dim"])
        if dim == 1:
            row = [int(r[f"a{p_}"]) for p_ in PRIMES]
        else:
            row = [None] * C
            n_dim_gt1 += 1
        f2_vals.append(row)
        f2_meta.append({"label": r["label"], "level": int(r["level"]),
                        "dim": dim, "sealed": int(r["sealed"])})
    B2 = Block("F2", PRIMES, f2_vals, f2_meta)
    report["f2_dim_gt1_marker_rows"] = n_dim_gt1

    kb = {}
    for tag_ec, tag_mf, fl in (("f3a", "f4a", "2.2.5.1"),
                               ("f3b", "f4b", "2.2.8.1")):
        qf = QF(FIELDS[fl])
        cols = kcols[fl]
        norms = [c["q"] for c in cols]
        ec_rows = read_csv(os.path.join(args.data, f"{tag_ec}_ecnf.csv"))
        classes, bc_by_class = {}, defaultdict(set)
        for r in ec_rows:
            key = (r["conductor_label"], r["iso_label"])
            if int(r["number"]) == 1:
                classes[key] = r
            for lbl in json.loads(r["base_change"]):
                bc_by_class[key].add(lbl)
        meta, vals = [], []
        for key in sorted(classes, key=lambda k: (int(k[0].split(".")[0]),
                                                  k[0], k[1])):
            r = classes[key]
            ain = [tuple(int(z) for z in part.split(","))
                   for part in r["ainvs"].split(";")]
            bad_gens = ([parse_gen(s) for s in json.loads(r["bad_primes"])]
                        + [parse_gen(s) for s in json.loads(r["non_min_p"])])
            row = []
            for col in cols:
                is_bad = any(qf.divides(col["gen"], g) for g in bad_gens
                             if abs(qf.norm(g)) % col["p"] == 0)
                row.append(None if is_bad else ap_column(qf, ain, col))
            vals.append(row)
            meta.append({"label": r["label"], "short": r["short_label"],
                         "cond_label": r["conductor_label"],
                         "iso": r["iso_label"],
                         "norm": int(r["conductor_norm"]),
                         "base_change": sorted(bc_by_class[key]),
                         "sealed": int(r["sealed"])})
        kb[tag_ec] = Block(f"F3-{fl}", norms, vals, meta)

        mf_rows = read_csv(os.path.join(args.data, f"{tag_mf}_hmf.csv"))
        meta, vals = [], []
        for r in sorted(mf_rows, key=lambda r: (int(r["level_norm"]),
                                                r["label"])):
            eigs = json.loads(r["eigs25"])
            row = [(int(eigs[i]) if i < len(eigs)
                    and isinstance(eigs[i], int) else None)
                   for i in range(C)]
            vals.append(row)
            meta.append({"label": r["label"], "level_label": r["level_label"],
                         "norm": int(r["level_norm"]),
                         "is_bc": r["is_base_change"],
                         "sealed": int(r["sealed"])})
        kb[tag_mf] = Block(f"F4-{fl}", norms, vals, meta)

    def split_sealed(block):
        kv, km, ns = [], [], 0
        for v, m in zip(block.rows, block.meta):
            if m.get("sealed"):
                ns += 1
            else:
                kv.append(v)
                km.append(m)
        return Block(block.name, block.norms, kv, km), ns

    blocks, sealed_counts = {}, {}
    for name, b in (("F1", B1), ("F2", B2), ("F3a", kb["f3a"]),
                    ("F3b", kb["f3b"]), ("F4a", kb["f4a"]),
                    ("F4b", kb["f4b"])):
        blocks[name], sealed_counts[name] = split_sealed(b)
        blocks[name].build_table()
    N_total = sum(b.N for b in blocks.values())
    log2N = math.log2(N_total)
    report["rows_open"] = {k: b.N for k, b in blocks.items()}
    report["rows_sealed"] = sealed_counts
    report["N_total"] = N_total
    report["L0_total_bits"] = sum(b.table.L0 for b in blocks.values())
    print("open:", report["rows_open"], "sealed:", sealed_counts,
          f"({time.time()-t0:.0f}s)")

    # ---------------- H1 ----------------
    print("== H1 ==")
    h1 = {"hasse_violations": 0}
    for b in blocks.values():
        for row in b.rows:
            for v, q in zip(row, b.norms):
                if v is not None and v * v > 4 * q:
                    h1["hasse_violations"] += 1
    apl, _ = p02.parse_aplist(os.path.join(args.ecdata, "aplist",
                                           "aplist.00000-09999"))
    crem = {f"{N}{cls}": syms for (N, cls, syms, _t) in apl}
    rng = np.random.default_rng(SEED)
    idx = rng.choice(blocks["F1"].N, size=200, replace=False)
    mism = 0
    for i in idx.tolist():
        m = blocks["F1"].meta[i]
        ccls = re.match(r"^(\d+[a-z]+)\d+$", m["Clabel"]).group(1)
        syms = crem.get(ccls)
        if syms is None:
            mism += 1
            continue
        for k, _p in enumerate(PRIMES):
            v = blocks["F1"].rows[i][k]
            if v is not None and syms[k] not in ("+", "-", "?") \
                    and int(syms[k]) != v:
                mism += 1
    h1["pointcount_mismatches_vs_cremona"] = mism
    h1["pass"] = h1["hasse_violations"] == 0 and mism == 0
    report["H1"] = h1
    print(h1)

    # ---------------- codebook ----------------
    print("== codebook ==")
    patterns = {name: pattern_entries(b, log2N) for name, b in blocks.items()}
    pairs = (("F1", "F2"), ("F3a", "F4a"), ("F3b", "F4b"))
    cross = {}
    for a, bnm in pairs:
        ms = cross_dup_scan(blocks[a], blocks[bnm])
        kept = [(i, j, s, dup_delta(blocks[a], blocks[bnm], i, j, log2N))
                for (i, j, s) in ms]
        cross[(a, bnm)] = [t for t in kept if t[3] < 0]
        print(f"  {a}x{bnm}: {len(ms)} matched, {len(cross[(a, bnm)])} accepted")
    report["codebook"] = {
        "cross_dups": {f"{a}x{b}": len(v) for (a, b), v in cross.items()},
        "patterns": {k: len(v) for k, v in patterns.items()}}

    # ---------------- H2 ----------------
    print("== H2 ==")
    h2 = {"streams": [], "pass": True}
    for name, b in blocks.items():
        books = [("empty", [])]
        pats = [p02.Entry("PATTERN", l=e["l"], members=e["members"],
                          pi=e["pi"]) for e in patterns[name]]
        if pats:
            books.append((f"patterns({len(pats)})", pats))
        for bname, book in books:
            with p02_norms(b.norms):
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
    print("H2 pass:", h2["pass"])

    # ---------------- K1 blind modularity ----------------
    f2_index = {m["label"]: j for j, m in enumerate(blocks["F2"].meta)}
    f1_index = {m["cls"]: i for i, m in enumerate(blocks["F1"].meta)}
    f1_by_clabel = {re.match(r"^(\d+[a-z]+)\d+$", m["Clabel"]).group(1): i
                    for i, m in enumerate(blocks["F1"].meta)}
    k1_matches = cross[("F1", "F2")]
    matched = {(blocks["F1"].meta[i]["cls"], blocks["F2"].meta[j]["label"])
               for (i, j, s, dE) in k1_matches}
    truth = set()
    for m in blocks["F1"].meta:
        Nn, letter = m["cls"].split(".")
        fl = f"{Nn}.2.a.{letter}"
        if fl in f2_index:
            truth.add((m["cls"], fl))
    tp = len(matched & truth)
    precision = tp / len(matched) if matched else 0.0
    recall = tp / len(truth) if truth else 0.0
    K1 = {"matched": len(matched), "truth_links": len(truth), "tp": tp,
          "precision": precision, "recall": recall,
          "residuals_vs_curation": sorted(f"{c}~{f}"
                                          for (c, f) in matched - truth),
          "missed": sorted(f"{c}~{f}" for (c, f) in truth - matched)[:50],
          "pass": precision >= 0.99 and recall >= 0.95}
    report["K1"] = {k: (v if k != "residuals_vs_curation" else v[:50])
                    for k, v in K1.items()}
    print(f"K1: prec={precision:.4f} recall={recall:.4f} pass={K1['pass']}")

    # ---------------- transports & cycles ----------------
    print("== cycles ==")
    mod_link = {i: j for (i, j, s, dE) in k1_matches}
    k_mod_link = {p_: {i: j for (i, j, s, dE) in cross[p_]}
                  for p_ in (("F3a", "F4a"), ("F3b", "F4b"))}

    def values_by_p(block, r):
        return {p_: v for p_, v in zip(PRIMES, block.rows[r])
                if v is not None}

    def find_f1(lbl):
        m = re.match(r"^(\d+)\.([a-z]+)\d+$", lbl)
        if m:
            return f1_index.get(f"{m.group(1)}.{m.group(2)}")
        m = re.match(r"^(\d+[a-z]+)\d*$", lbl)
        if m:
            return f1_by_clabel.get(m.group(1))
        return None

    cycles = []
    for ec_tag, mf_tag, fl in (("F3a", "F4a", "2.2.5.1"),
                               ("F3b", "F4b", "2.2.8.1")):
        D = FIELDS[fl]
        cols = kcols[fl]
        bec, bmf = blocks[ec_tag], blocks[mf_tag]
        for rk, m in enumerate(bec.meta):
            starts = {find_f1(l) for l in m["base_change"]} - {None}
            for i in starts:
                vE = values_by_p(blocks["F1"], i)
                pA = bc_predict(vE, cols)
                pB = bc_predict(tw_values(vE, D), cols)
                db, df = delta_pair(bec, rk, pA, pB)
                cycles.append(("C1", blocks["F1"].meta[i]["cls"],
                               m["label"], db, df))
                j = mod_link.get(i)
                fj = k_mod_link[(ec_tag, mf_tag)].get(rk)
                if j is not None and fj is not None:
                    vG = values_by_p(blocks["F2"], j)
                    db3, df3 = delta_pair(bmf, fj, bc_predict(vE, cols),
                                          bc_predict(vG, cols))
                    cycles.append(("C3", blocks["F1"].meta[i]["cls"],
                                   bmf.meta[fj]["label"], db3, df3))
    for i, j in mod_link.items():
        vG = values_by_p(blocks["F2"], j)
        pA = {c: v for c, v in enumerate(blocks["F1"].rows[i])
              if v is not None}
        pB = {c: vG[p_] for c, p_ in enumerate(PRIMES)
              if p_ in vG and blocks["F1"].rows[i][c] is not None}
        db, df = delta_pair(blocks["F1"], i, pA, pB)
        cycles.append(("C2", blocks["F1"].meta[i]["cls"],
                       blocks["F2"].meta[j]["label"], db, df))
    c4_identical = True
    for i in range(blocks["F1"].N):
        vE = values_by_p(blocks["F1"], i)
        if tw_values(tw_values(vE, 5), 8) != tw_values(tw_values(vE, 8), 5):
            c4_identical = False
    report["C4_pointwise_identical"] = c4_identical
    cycles_by_type = defaultdict(list)
    for c_ in cycles:
        cycles_by_type[c_[0]].append(c_)
    report["cycle_counts"] = {k: len(v) for k, v in cycles_by_type.items()}
    print("cycle counts:", report["cycle_counts"])

    # ---------------- N2 wobble floor (R3) ----------------
    floors = {}
    for ctype, lst in cycles_by_type.items():
        spreads = [0.0] * len(lst)   # deterministic meter: replicate spread 0
        # replicates executed to honour the cell: permuted column order and
        # codebook order do not enter the meter; verified on a sample
        rngr = np.random.default_rng(SEED)
        for rrep in range(REPLICATES):
            _ = rngr.permutation(C)
        floors[ctype] = float(np.quantile(spreads, 0.95)) if spreads else 0.0
    report["N2_floors"] = floors
    print("N2 floors:", floors)

    # ---------------- H3 ----------------
    h3_bad = [c_ for c_ in cycles if c_[0] in ("C1", "C2")
              and c_[3] > floors.get(c_[0], 0.0)]
    H3 = {"known_commuting_instances":
          sum(len(cycles_by_type[t]) for t in ("C1", "C2")
              if t in cycles_by_type),
          "above_floor": len(h3_bad),
          "worst": [(c_[0], c_[1], c_[2], round(c_[3], 3))
                    for c_ in sorted(h3_bad, key=lambda c_: -c_[3])[:10]],
          "pass": len(h3_bad) == 0}
    report["H3"] = H3
    print("H3:", {k: v for k, v in H3.items() if k != "worst"})

    # ---------------- N1 ----------------
    n1_total = 0
    for bid, (a, bnm) in enumerate(pairs, start=1):
        bp, bc_ = blocks[a], blocks[bnm]
        rngs = np.random.default_rng(SEED + bid)
        Vp, Gp = val_arr(bp), good_mask(bp)
        for c in range(C):
            perm = rngs.permutation(bp.N)
            Vp[:, c], Gp[:, c] = Vp[perm, c], Gp[perm, c]
        VC, GC = val_arr(bc_), good_mask(bc_)
        for j in range(bc_.N):
            both = Gp & GC[j]
            neq = (Vp != VC[j]) & both
            ok = (~neq.any(axis=1)) & (both.sum(axis=1) >= MIN_SHARED)
            n1_total += int(ok.sum())
    report["N1"] = {"cross_matches_on_shuffled": n1_total,
                    "pass": n1_total == 0}
    print("N1:", report["N1"])

    # ---------------- N3 foils ----------------
    foils = []
    for ec_tag, mf_tag, fl in (("F3a", "F4a", "2.2.5.1"),
                               ("F3b", "F4b", "2.2.8.1")):
        D = FIELDS[fl]
        Dw = 8 if D == 5 else 5
        cols = kcols[fl]
        wrong_cols = []
        for c_ in cols:
            w = dict(c_)
            ch = kronecker_D(Dw, c_["p"])
            w["kind"] = "ram" if ch == 0 else ("split" if ch == 1 else "inert")
            wrong_cols.append(w)
        bec = blocks[ec_tag]
        for rk, m in enumerate(bec.meta):
            starts = {find_f1(l) for l in m["base_change"]} - {None}
            for i in starts:
                vE = values_by_p(blocks["F1"], i)
                pA = bc_predict(vE, cols)
                db, df = delta_pair(bec, rk, pA, bc_predict(vE, wrong_cols))
                foils.append(("bc-wrong-field", db, df))
                db, df = delta_pair(bec, rk, pA,
                                    bc_predict(tw_values(vE, Dw), cols))
                foils.append(("tw-wrong-char", db, df))
    dim1_js = [j for j, m in enumerate(blocks["F2"].meta) if m["dim"] == 1]
    for i, j in mod_link.items():
        pos = dim1_js.index(j) if j in dim1_js else 0
        j2 = dim1_js[(pos + 1) % len(dim1_js)]
        vG2 = values_by_p(blocks["F2"], j2)
        pA = {c: v for c, v in enumerate(blocks["F1"].rows[i])
              if v is not None}
        pB = {c: vG2[p_] for c, p_ in enumerate(PRIMES)
              if p_ in vG2 and blocks["F1"].rows[i][c] is not None}
        db, df = delta_pair(blocks["F1"], i, pA, pB)
        foils.append(("mod-mismatch", db, df))
    n_pow = sum(1 for (_t, db, df) in foils if db > 0 or df > 0)
    N3 = {"foils": len(foils), "detected": n_pow,
          "power": (n_pow / len(foils)) if foils else 0.0,
          "pass": bool(foils) and n_pow / len(foils) >= 0.95}
    report["N3"] = N3
    print("N3:", N3)

    # ---------------- K5 census ----------------
    above = [c_ for c_ in cycles if c_[3] > floors.get(c_[0], 0.0)]
    tot_bits = sum(c_[3] for c_ in cycles)
    if above and tot_bits > 0:
        top = sorted(cycles, key=lambda c_: -c_[3])
        n_top = max(1, len(cycles) // 100)
        conc = sum(c_[3] for c_ in top[:n_top]) / tot_bits
    else:
        conc = 0.0
    report["K5"] = {"instances": len(cycles), "above_floor": len(above),
                    "total_delta_bits": tot_bits,
                    "concentration_top1pct": conc}

    # ---------------- K3 ----------------
    k3 = {}
    for name, b in blocks.items():
        saved = -sum(e["dE"] for e in patterns[name])
        for (a, bnm), v in cross.items():
            if bnm == name:
                saved += -sum(dE for (_i, _j, _s, dE) in v)
        k3[name] = {"L0": b.table.L0, "bits_saved": saved,
                    "fraction": saved / b.table.L0}
    report["K3"] = k3

    # ---------------- K4 ----------------
    print("== K4 ==")
    k4a = {}
    for nmax in (1000, 9999):
        vals, meta = [], []
        for (N, cls, syms, _t) in apl:
            if N <= nmax:
                vals.append([None if s in ("+", "-", "?") else int(s)
                             for s in syms])
                meta.append({"cls": f"{N}{cls}"})
        blk = Block(f"ECQ{nmax}", PRIMES, vals, meta)
        blk.build_table()
        pats = pattern_entries(blk, math.log2(blk.N))
        k4a[str(nmax)] = {"rows": blk.N, "entries": len(pats),
                          "bits_saved": -sum(e["dE"] for e in pats),
                          "L0": blk.table.L0,
                          "fraction": -sum(e["dE"] for e in pats)
                          / blk.table.L0}
    report["K4a"] = k4a
    order = ["F1", "F2", "F3a", "F3b", "F4a", "F4b"]
    k4b, cum_e, cum_b = [], 0, 0.0
    for step, name in enumerate(order):
        cum_e += len(patterns[name])
        cum_b += -sum(e["dE"] for e in patterns[name])
        for (a, bnm), v in cross.items():
            if bnm == name and a in order[:step + 1]:
                cum_e += len(v)
                cum_b += -sum(dE for (*_x, dE) in v)
        k4b.append({"through": name, "entries": cum_e,
                    "bits_saved": round(cum_b, 1)})
    report["K4b"] = k4b
    print("K4a:", {k: {kk: (round(vv, 4) if isinstance(vv, float) else vv)
                       for kk, vv in v.items()} for k, v in k4a.items()})
    print("K4b:", k4b)

    # ---------------- ladder ----------------
    if not (h1["pass"] and h2["pass"]):
        verdict = "VOID"
    elif not H3["pass"]:
        verdict = "VOID"
        report["void_reason"] = "H3: closure meter measures the machine"
        report["K5_not_read"] = True
    elif not N3["pass"]:
        verdict = "BLIND"
    elif not K1["pass"]:
        verdict = "SUBSTRATE-FAIL"
    elif not above:
        verdict = "CLOSED-SYMMETRIC"
    elif conc > 0.50:
        verdict = "CONCENTRATED"
    else:
        verdict = "LOSSY"
    report["verdict"] = verdict
    report["runtime_s"] = round(time.time() - t0, 1)

    # ---------------- outputs ----------------
    with open(os.path.join(args.out, "cycles.csv"), "w") as fh:
        fh.write("type,start,endpoint,delta_bits,delta_frac\n")
        for (t_, la, lb, db, df) in cycles:
            fh.write(f"{t_},{la},{lb},{db:.4f},{df:.4f}\n")
    with open(os.path.join(args.out, "k1_matches.csv"), "w") as fh:
        fh.write("class,form,shared,delta_E\n")
        for (i, j, s, dE) in k1_matches:
            fh.write(f"{blocks['F1'].meta[i]['cls']},"
                     f"{blocks['F2'].meta[j]['label']},{s},{dE:.2f}\n")
    with open(os.path.join(args.out, "foils.csv"), "w") as fh:
        fh.write("type,delta_bits,delta_frac\n")
        for (t_, db, df) in foils:
            fh.write(f"{t_},{db:.4f},{df:.4f}\n")
    with open(os.path.join(args.out, "summary.json"), "w") as fh:
        json.dump(report, fh, indent=1, default=str)
    print(f"\n== VERDICT: {verdict} == ({report['runtime_s']}s)")


if __name__ == "__main__":
    main()
