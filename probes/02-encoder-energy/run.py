#!/usr/bin/env python3
"""Probe 02 — the encoder and its energy. Mechanical implementation of PREREG.md.

Usage:
    python3 run.py --ecdata /path/to/ecdata --out output/

Python 3 + numpy only. Data pinned at ecdata
25cec5ecfec8b9f016eb1631ac633194c2bed39f, conductors 1-9999 only (SEAL.md).

Implementation readings forced by the prereg text (documented, not tunable):

R1. The cell code is the two-branch flag scheme of PREREG §4.3 at every claimed
    cell. A claimed cell holding a bad marker (+ - ?) cannot conform (no residue
    class contains a marker: Q_c sums integers only), so it takes the escape
    branch: flag + literal at the full baseline distribution. This is the only
    reading a *real decoder* can implement — the decoder cannot know a cell
    holds a marker before decoding it — and it is the prereg's own mechanism
    ("a false claim is not an error; it is 6 bits"). "Fall through to baseline"
    = the fall-through (escape) branch.
R2. Claim sets depend only on decoder-available data: DUP claims all 25 cells
    of its target row; BRIDGE claims cells where the (already decoded) parent
    cell is an integer, c != column(l); PATTERN/HUB/TEST claim all cells of
    each member with c != column(l). Ownership: first entry in codebook order,
    per §4.3.
R3. Acceptance/ΔE gates (K1, K2's ΔE<0, N1, N2, K7) use the full
    coder-consistent energy of §4.4, escapes included. The quantities named
    A_PAT / A_HUB / per-member saving in K2-K4 are computed by their defining
    formulas in §5.1 (S-based: savings over conforming good cells), evaluated
    on measured marginals. Both accountings are recorded.
R4. H2's (T, K-book) matrix runs on the well-typed combinations: each book on
    its defining substrate plus the empty book on both substrates. The literal
    cross product would index curve rows into the class table (row spaces
    differ) and is undefined; disclosed, not repaired.
R5. N3's byte-stream comparison runs on the order-independent books (baseline,
    PATTERN, HUB); BRIDGE/DUP reference previously decoded rows, so an
    arbitrary row permutation breaks decode order. Analytic-total equality is
    verified for every book. "Exact" equality is exact up to float summation
    order (< 1e-6 bits), reported.
"""

import argparse
import hashlib
import importlib.util
import json
import math
import os
import subprocess
import sys
import time
from collections import defaultdict

import numpy as np

# ------------------------------------------------------------ constants (§7)
PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
          67, 71, 73, 79, 83, 89, 97]
C = 25
MODULI = [5, 7, 11, 13]
SEED = 20260728
EPS_BITS_MAIN = 6                      # eps0 = 2^-6
EPS_BITS_SWEEP = [4, 6, 8]             # K6
H_BITS = 5.0                           # 3 type + 2 modulus
PHI_SIZE = 14
FOIL_M = 1000
PINNED_SHA = "25cec5ecfec8b9f016eb1631ac633194c2bed39f"
TOT = 1 << 18                          # coder frequency resolution (implementation)
MARKERS = ["+", "-", "?"]
PI_ORDER = ["EIS", "0", "1", "ID"]     # 2-bit pattern ids, |Pi| = 4
PI_FUN = {"EIS": lambda p: p + 1, "0": lambda p: 0, "1": lambda p: 1,
          "ID": lambda p: p}
PROBE01_PUBLISHED = {5: 3521, 7: 896, 11: 32, 13: 5}


def gamma_bits(k):
    if k <= 0:
        return 1.0
    return 2 * math.floor(math.log2(k)) + 1


def flag_bits(eps_bits):
    return -math.log2(1.0 - 2.0 ** (-eps_bits))


def seal_guard(path):
    """SEAL.md: only conductor-range 00000-09999 files may be read."""
    if not os.path.basename(path).endswith("00000-09999"):
        raise RuntimeError(f"SEAL VIOLATION attempted: {path}")
    return path


# ---------------------------------------------------------------- parsing

def parse_aplist(path):
    rows, raw_lines = [], []
    with open(seal_guard(path)) as fh:
        for line in fh:
            tok = line.split()
            if not tok:
                continue
            raw_lines.append(tok)
            N, cls = int(tok[0]), tok[1]
            syms = tok[2:2 + C]
            trailing = tok[2 + C:]
            for t in trailing:
                if not (t[0] in "+-?" and t[1] == "(" and t[-1] == ")"):
                    raise ValueError(f"bad trailing token {t!r}")
                p = int(t[2:-1])
                if p <= 97 or N % p != 0:
                    raise ValueError(f"inconsistent trailing marker {t!r}")
            for s in syms:
                if s not in MARKERS:
                    int(s)
            rows.append((N, cls, syms, trailing))
    return rows, raw_lines


def serialise_aplist_row(row):
    N, cls, syms, trailing = row
    out = [str(N), cls]
    for s in syms:
        out.append(s if s in MARKERS else str(int(s)))
    out.extend(trailing)
    return out


def parse_allcurves(path):
    rows, raw_lines = [], []
    with open(seal_guard(path)) as fh:
        for line in fh:
            tok = line.split()
            if not tok:
                continue
            raw_lines.append(tok)
            rows.append((int(tok[0]), tok[1], int(tok[2]), json.loads(tok[3]),
                         [int(t) for t in tok[4:]]))
    return rows, raw_lines


def serialise_allcurves_row(row):
    N, cls, num, ainvs, rest = row
    return [str(N), cls, str(num),
            "[" + ",".join(str(a) for a in ainvs) + "]"] + [str(x) for x in rest]


def parse_allisog_5isog(path):
    out = {}
    with open(seal_guard(path)) as fh:
        for line in fh:
            tok = line.split()
            if not tok:
                continue
            mat = json.loads(tok[5])
            out[(int(tok[0]), tok[1])] = any(
                d % 5 == 0 for rowm in mat for d in rowm)
    return out


# ---------------------------------------------------------------- table

class Table:
    """Canonical table over per-column alphabets with KT baseline (§4.1-4.2)."""

    def __init__(self, name, sym_rows):
        self.name = name
        self.N = len(sym_rows)
        self.log2N = math.log2(self.N)
        self.alpha, self.sym_index = [], []
        for c in range(C):
            B = math.isqrt(4 * PRIMES[c])            # floor(2*sqrt(p))
            symbols = [str(v) for v in range(-B, B + 1)] + MARKERS
            self.alpha.append(symbols)
            self.sym_index.append({s: i for i, s in enumerate(symbols)})
        self.nsym = [len(a) for a in self.alpha]
        self.SYM = np.zeros((self.N, C), dtype=np.int16)
        self.ISINT = np.zeros((self.N, C), dtype=bool)
        self.VAL = np.zeros((self.N, C), dtype=np.int16)
        for r, syms in enumerate(sym_rows):
            for c, s in enumerate(syms):
                idx = self.sym_index[c].get(s)
                if idx is None:
                    raise ValueError(f"{name}: symbol {s!r} outside alphabet c={c}")
                self.SYM[r, c] = idx
                if s not in MARKERS:
                    self.ISINT[r, c] = True
                    self.VAL[r, c] = int(s)
        self.counts = [np.bincount(self.SYM[:, c], minlength=self.nsym[c])
                       for c in range(C)]
        self.q, self.logq_sym = [], []
        for c in range(C):
            qc = (self.counts[c] + 0.5) / (self.N + 0.5 * self.nsym[c])
            self.q.append(qc)
            self.logq_sym.append(-np.log2(qc))
        self.logq_cell = np.zeros((self.N, C))
        for c in range(C):
            self.logq_cell[:, c] = self.logq_sym[c][self.SYM[:, c]]
        self.model_bits_col = [(self.nsym[c] - 1) / 2.0 * self.log2N
                               for c in range(C)]
        self.L0_model = float(sum(self.model_bits_col))
        # order-independent data term: per-column counts . logq
        self.L0_data = float(sum(float(np.dot(self.counts[c], self.logq_sym[c]))
                                 for c in range(C)))
        self.L0 = self.L0_model + self.L0_data
        self.Q, self.negLogQ = {}, {}
        for l in MODULI:
            Ql = np.zeros((C, l))
            for c in range(C):
                B = (self.nsym[c] - 3 - 1) // 2
                for i, v in enumerate(range(-B, B + 1)):
                    Ql[c, v % l] += self.q[c][i]
            self.Q[l] = Ql
            with np.errstate(divide="ignore"):
                nl = -np.log2(Ql)
            self.negLogQ[l] = np.where(np.isfinite(nl), nl, 0.0)  # empty class:
            # -log2 Q only ever multiplies a conform indicator, and conforming
            # in an empty class is impossible, so 0 is a safe placeholder.
            self.Q0 = None


# ---------------------------------------------------------------- entries

class Entry:
    __slots__ = ("etype", "l", "i", "j", "members", "pi", "rh", "nX", "label")

    def __init__(self, etype, l=None, i=None, j=None, members=None, pi=None,
                 rh=None, nX=0, label=""):
        self.etype, self.l, self.i, self.j = etype, l, i, j
        self.members, self.pi, self.rh, self.nX = members, pi, rh, nX
        self.label = label

    def cost(self, table):
        L = table.log2N
        if self.etype in ("DUP", "BRIDGE"):
            return H_BITS + 2 * L
        k = len(self.members)
        if self.etype == "PATTERN":
            return H_BITS + 2 + gamma_bits(k) + k * L
        if self.etype == "HUB":
            return H_BITS + gamma_bits(k) + k * L + 24 * math.log2(self.l)
        if self.etype == "TEST":
            return (H_BITS + 2 + math.log2(PHI_SIZE)
                    + gamma_bits(self.nX) + self.nX * L)
        raise ValueError(self.etype)

    def residue_row(self, table):
        """Length-C residue vector, -1 = column not claimed (reading R2)."""
        col_l = PRIMES.index(self.l) if self.l in PRIMES else -1
        rh = np.full(C, -1, dtype=np.int16)
        if self.etype == "BRIDGE":
            for c in range(C):
                if c != col_l and table.ISINT[self.i, c]:
                    rh[c] = int(table.VAL[self.i, c]) % self.l
        elif self.etype in ("PATTERN", "TEST"):
            fn = PI_FUN[self.pi]
            for c in range(C):
                if c != col_l:
                    rh[c] = fn(PRIMES[c]) % self.l
        elif self.etype == "HUB":
            rh = self.rh
        return rh

    def target_rows(self):
        if self.etype in ("DUP", "BRIDGE"):
            return [self.j]
        return list(self.members)


def make_hub(members, l, table, label=""):
    """r_h[c] = residue of the lowest-indexed member with an integer at c."""
    col_l = PRIMES.index(l)
    rh = np.full(C, -1, dtype=np.int16)
    for c in range(C):
        if c == col_l:
            continue
        rh[c] = 0
        for m in members:
            if table.ISINT[m, c]:
                rh[c] = int(table.VAL[m, c]) % l
                break
    return Entry("HUB", l=l, members=list(members), rh=rh, label=label)


def eval_mod_entry(table, e, eps_bits):
    """Standalone ΔE (true, escapes incl.) + S-based savings for a mod-l entry.

    Returns dict(delta_E, S, conform, escapes, per_row_S, per_row_deltaE)."""
    fb = flag_bits(eps_bits)
    rh = e.residue_row(table)
    cols = np.nonzero(rh >= 0)[0]
    rows = np.asarray(e.target_rows(), dtype=np.int64)
    if len(cols) == 0 or len(rows) == 0:
        z = np.zeros(len(rows))
        return {"delta_E": e.cost(table), "S": 0.0, "conform": 0, "escapes": 0,
                "per_row_S": z, "per_row_cells": z}
    V = table.VAL[np.ix_(rows, cols)]
    I = table.ISINT[np.ix_(rows, cols)]
    conform = I & ((np.mod(V, e.l)) == rh[cols][None, :])
    nlq = table.negLogQ[e.l][cols, rh[cols]][None, :]   # -log2 Q >= 0
    conf_rows = conform.sum(axis=1)
    esc_rows = (~conform).sum(axis=1)
    S_rows = (nlq * conform).sum(axis=1) - fb * conf_rows
    cells_rows = fb * conf_rows - (nlq * conform).sum(axis=1) + eps_bits * esc_rows
    return {"delta_E": float(e.cost(table) + cells_rows.sum()),
            "S": float(S_rows.sum()),
            "conform": int(conf_rows.sum()), "escapes": int(esc_rows.sum()),
            "per_row_S": S_rows, "per_row_cells": cells_rows}


def eval_dup_book(table, book, eps_bits):
    """Vectorized DUP-book evaluation (disjoint target rows).

    Per claimed cell: conform costs fb (delta fb - baseline), escape costs
    eps + baseline literal (delta +eps)."""
    fb = flag_bits(eps_bits)
    parents = np.array([e.i for e in book], dtype=np.int64)
    children = np.array([e.j for e in book], dtype=np.int64)
    conform = table.SYM[children] == table.SYM[parents]
    conf_rows = conform.sum(axis=1)
    esc_rows = C - conf_rows
    logq_children = table.logq_cell[children]
    conf_logq = (logq_children * conform).sum(axis=1)
    esc_logq = (logq_children * ~conform).sum(axis=1)
    delta_rows = fb * conf_rows - conf_logq + eps_bits * esc_rows
    abs_cost_rows = fb * conf_rows + eps_bits * esc_rows + esc_logq
    cost = np.full(len(book), H_BITS + 2 * table.log2N)
    return {"delta_E": float((cost + delta_rows).sum()),
            "per_entry_delta_E": cost + delta_rows,
            "cells_bits": float(abs_cost_rows.sum()),
            "conform": int(conf_rows.sum()), "escapes": int(esc_rows.sum())}


# ------------------------------------------------- ownership book evaluation

def build_claims(table, book):
    owner = np.full((table.N, C), -1, dtype=np.int32)
    rho = np.full((table.N, C), -1, dtype=np.int16)
    dup_parent = np.full(table.N, -1, dtype=np.int64)
    for idx, e in enumerate(book):
        if e.etype == "DUP":
            r = e.j
            free = owner[r] == -1
            owner[r][free] = idx
            dup_parent[r] = e.i
            continue
        rh = e.residue_row(table)
        cols = np.nonzero(rh >= 0)[0]
        for r in e.target_rows():
            free = owner[r][cols] == -1
            cc = cols[free]
            owner[r][cc] = idx
            rho[r][cc] = rh[cc]
    return owner, rho, dup_parent


def book_cells_delta(table, book, owner, rho, dup_parent, eps_bits):
    """Sum over owned claimed cells of (cost - baseline); python per cell."""
    fb = flag_bits(eps_bits)
    total = 0.0
    esc = conf = 0
    rs, cs = np.nonzero(owner >= 0)
    ll = np.array([getattr(e, "l") or 0 for e in book])
    types = [e.etype for e in book]
    for r, c in zip(rs.tolist(), cs.tolist()):
        idx = owner[r, c]
        if types[idx] == "DUP":
            if table.SYM[r, c] == table.SYM[dup_parent[r], c]:
                total += fb - table.logq_cell[r, c]
                conf += 1
            else:
                total += eps_bits
                esc += 1
        else:
            l = int(ll[idx])
            p = rho[r, c]
            if table.ISINT[r, c] and int(table.VAL[r, c]) % l == p:
                total += fb - table.negLogQ[l][c, p]
                conf += 1
            else:
                total += eps_bits
                esc += 1
    return total, conf, esc


def book_delta_E_owned(table, book, eps_bits):
    owner, rho, dup_parent = build_claims(table, book)
    cells, conf, esc = book_cells_delta(table, book, owner, rho, dup_parent,
                                        eps_bits)
    cost = float(sum(e.cost(table) for e in book))
    return {"cost_bits": cost, "cells_delta_bits": cells,
            "delta_E_bits": cost + cells, "conform": conf, "escapes": esc,
            "owner": owner, "rho": rho, "dup_parent": dup_parent}


# ---------------------------------------------------------------- coder

MASK64 = (1 << 64) - 1
TOPV = 1 << 56


class RangeEncoder:
    """64-bit carry-handling range coder; truncation loss ~2^-38/symbol."""

    def __init__(self):
        self.low = 0
        self.range = MASK64
        self.out = bytearray()

    def _carry(self):
        i = len(self.out) - 1
        while i >= 0:
            if self.out[i] == 0xFF:
                self.out[i] = 0
                i -= 1
            else:
                self.out[i] += 1
                return
        self.out.insert(0, 1)   # cannot occur in practice; correctness guard

    def encode(self, cum, freq):
        r = self.range // TOT
        new_low = self.low + r * cum
        if new_low > MASK64:
            self._carry()
            new_low &= MASK64
        self.low = new_low
        self.range = r * freq
        while self.range < TOPV:
            self.out.append((self.low >> 56) & 0xFF)
            self.low = (self.low << 8) & MASK64
            self.range <<= 8

    def finish(self):
        for _ in range(8):
            self.out.append((self.low >> 56) & 0xFF)
            self.low = (self.low << 8) & MASK64
        return bytes(self.out)


class RangeDecoder:
    def __init__(self, data):
        self.data = data
        self.pos = 8
        self.range = MASK64
        self.code = int.from_bytes(data[:8], "big")

    def _byte(self):
        b = self.data[self.pos] if self.pos < len(self.data) else 0
        self.pos += 1
        return b

    def decode_target(self):
        self.r = self.range // TOT
        t = self.code // self.r
        return TOT - 1 if t >= TOT else t

    def consume(self, cum, freq):
        self.code -= self.r * cum
        self.range = self.r * freq
        while self.range < TOPV:
            self.code = ((self.code << 8) | self._byte())
            self.range <<= 8


class FreqTable:
    def __init__(self, probs):
        p = np.asarray(probs, dtype=float)
        f = np.maximum(1, np.round(p * TOT).astype(np.int64))
        diff = int(TOT - f.sum())
        order = np.argsort(-f)
        i = 0
        while diff != 0:
            j = int(order[i % len(f)])
            step = 1 if diff > 0 else -1
            if f[j] + step >= 1:
                f[j] += step
                diff -= step
            i += 1
        self.freq = f
        self.cum = np.concatenate([[0], np.cumsum(f)])

    def find(self, target):
        return int(np.searchsorted(self.cum, target, side="right") - 1)


def coder_self_test(rng):
    for _ in range(5):
        nsym = int(rng.integers(2, 50))
        ft = FreqTable(rng.dirichlet(np.ones(nsym) * 0.3))
        syms = rng.integers(0, nsym, size=3000)
        enc = RangeEncoder()
        for s in syms:
            enc.encode(int(ft.cum[s]), int(ft.freq[s]))
        blob = enc.finish()
        dec = RangeDecoder(blob)
        for s in syms:
            g = ft.find(dec.decode_target())
            assert g == s, "coder self-test failed"
            dec.consume(int(ft.cum[g]), int(ft.freq[g]))


class StreamCodec:
    def __init__(self, table, book, eps_bits):
        self.t = table
        self.book = book
        self.eps_bits = eps_bits
        self.col_ft = [FreqTable(table.q[c]) for c in range(C)]
        e0 = 2.0 ** (-eps_bits)
        self.flag_ft = FreqTable([1.0 - e0, e0])
        self._restricted = {}
        self.owner, self.rho, self.dup_parent = build_claims(table, book)

    def restricted_ft(self, c, l, p):
        key = (c, l, p)
        hit = self._restricted.get(key)
        if hit is not None:
            return hit
        B = (self.t.nsym[c] - 3 - 1) // 2
        idxs = [i for i, v in enumerate(range(-B, B + 1)) if v % l == p]
        probs = np.array([self.t.q[c][i] for i in idxs])
        ft = FreqTable(probs / probs.sum())
        self._restricted[key] = (ft, idxs)
        return ft, idxs

    def encode(self):
        t, owner, rho, dup_parent = self.t, self.owner, self.rho, self.dup_parent
        enc = RangeEncoder()
        fF, fC = self.flag_ft, self.col_ft
        for r in range(t.N):
            SYMr = t.SYM[r]
            OWNr = owner[r]
            for c in range(C):
                idx = OWNr[c]
                s = int(SYMr[c])
                if idx < 0:
                    ft = fC[c]
                    enc.encode(int(ft.cum[s]), int(ft.freq[s]))
                    continue
                e = self.book[idx]
                if e.etype == "DUP":
                    if s == int(t.SYM[dup_parent[r], c]):
                        enc.encode(int(fF.cum[0]), int(fF.freq[0]))
                    else:
                        enc.encode(int(fF.cum[1]), int(fF.freq[1]))
                        ft = fC[c]
                        enc.encode(int(ft.cum[s]), int(ft.freq[s]))
                    continue
                l, p = e.l, int(rho[r, c])
                if bool(t.ISINT[r, c]) and int(t.VAL[r, c]) % l == p:
                    enc.encode(int(fF.cum[0]), int(fF.freq[0]))
                    ft, idxs = self.restricted_ft(c, l, p)
                    j = idxs.index(s)
                    enc.encode(int(ft.cum[j]), int(ft.freq[j]))
                else:
                    enc.encode(int(fF.cum[1]), int(fF.freq[1]))
                    ft = fC[c]
                    enc.encode(int(ft.cum[s]), int(ft.freq[s]))
        return enc.finish()

    def decode(self, blob):
        t, owner, rho, dup_parent = self.t, self.owner, self.rho, self.dup_parent
        SYM = np.zeros((t.N, C), dtype=np.int16)
        dec = RangeDecoder(blob)
        fF, fC = self.flag_ft, self.col_ft

        def read(ft):
            g = ft.find(dec.decode_target())
            dec.consume(int(ft.cum[g]), int(ft.freq[g]))
            return g

        for r in range(t.N):
            for c in range(C):
                idx = owner[r, c]
                if idx < 0:
                    s = read(fC[c])
                else:
                    e = self.book[idx]
                    flag = read(fF)
                    if e.etype == "DUP":
                        s = int(SYM[dup_parent[r], c]) if flag == 0 else read(fC[c])
                    elif flag == 0:
                        ft, idxs = self.restricted_ft(c, e.l, int(rho[r, c]))
                        s = idxs[read(ft)]
                    else:
                        s = read(fC[c])
                SYM[r, c] = s
        return SYM


# ------------------------------------------------------------ probe-01 groups

def load_probe01_module():
    spec = importlib.util.spec_from_file_location(
        "p01run", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "..", "01-congruence-census", "run.py"))
    p01 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(p01)
    return p01


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ecdata", required=True)
    ap.add_argument("--out", default="output")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    report = {"seed": SEED, "eps_bits_main": EPS_BITS_MAIN, "TOT": TOT,
              "readings": ["R1", "R2", "R3", "R4", "R5"]}
    try:
        sha = subprocess.run(["git", "-C", args.ecdata, "rev-parse", "HEAD"],
                             capture_output=True, text=True).stdout.strip()
    except OSError:
        sha = "unknown"
    report["ecdata_sha"] = sha
    if sha != PINNED_SHA:
        print(f"WARNING: ecdata HEAD {sha} != pinned {PINNED_SHA}")

    coder_self_test(np.random.default_rng(12345))
    print("coder self-test OK")

    # ---------------- load + H1 ----------------
    print("== load ==")
    ap_rows, ap_raw = parse_aplist(
        os.path.join(args.ecdata, "aplist", "aplist.00000-09999"))
    ac_rows, ac_raw = parse_allcurves(
        os.path.join(args.ecdata, "allcurves", "allcurves.00000-09999"))
    phi5 = parse_allisog_5isog(
        os.path.join(args.ecdata, "allisog", "allisog.00000-09999"))

    h1 = {}
    h1["aplist_roundtrip_failures"] = sum(
        1 for row, raw in zip(ap_rows, ap_raw)
        if serialise_aplist_row(row) != raw)
    h1["allcurves_roundtrip_failures"] = sum(
        1 for row, raw in zip(ac_rows, ac_raw)
        if serialise_allcurves_row(row) != raw)
    lab = {(N, cls): i for i, (N, cls, _s, _t) in enumerate(ap_rows)}
    labels = [(N, cls) for (N, cls, _s, _t) in ap_rows]
    h1["dup_aplist_rows"] = len(ap_rows) - len(lab)
    h1["curves_without_class_row"] = sum(
        1 for (N, cls, *_r) in ac_rows if (N, cls) not in lab)
    classes_with_curves = {(N, cls) for (N, cls, *_r) in ac_rows}
    h1["class_rows_without_curve"] = sum(
        1 for key in lab if key not in classes_with_curves)

    tbl_A = Table("T_cls", [syms for (_N, _c, syms, _t) in ap_rows])
    crv_syms, crv_class_row, crv_is_rep, crv_labels = [], [], [], []
    for (N, cls, num, _ai, _rest) in ac_rows:
        i = lab[(N, cls)]
        crv_syms.append(ap_rows[i][2])
        crv_class_row.append(i)
        crv_is_rep.append(num == 1)
        crv_labels.append(f"{N}{cls}{num}")
    tbl_B = Table("T_crv", crv_syms)
    print(f"tables: A {tbl_A.N}x25 L0={tbl_A.L0:.0f}b | "
          f"B {tbl_B.N}x25 L0={tbl_B.L0:.0f}b")

    print("== probe-01 group reconstruction ==")
    p01 = load_probe01_module()
    vals01 = np.where(tbl_A.ISINT, tbl_A.VAL, 0).astype(np.int8)
    good01 = tbl_A.ISINT.copy()
    comps = {}
    for l in MODULI:
        kcol = PRIMES.index(l)
        cols = [k for k in range(C) if k != kcol]
        vm = (vals01.astype(np.int16) % l).astype(np.int8)[:, cols]
        comps[l], _stats = p01.scan(vm, good01[:, cols], tag=f"p01-l{l}",
                                    log=lambda *a, **k: None)
    recon = {l: len(comps[l]) for l in MODULI}
    h1["probe01_reconstruction"] = {str(l): recon[l] for l in MODULI}
    h1["probe01_reconstruction_ok"] = recon == PROBE01_PUBLISHED
    g143 = sorted(max(comps[5], key=len))
    h1["g143_size"] = len(g143)
    h1["g143_contains_11a"] = lab[(11, "a")] in g143
    h1_pass = (h1["aplist_roundtrip_failures"] == 0
               and h1["allcurves_roundtrip_failures"] == 0
               and h1["dup_aplist_rows"] == 0
               and h1["curves_without_class_row"] == 0
               and h1["class_rows_without_curve"] == 0
               and h1["probe01_reconstruction_ok"]
               and h1["g143_size"] == 143 and h1["g143_contains_11a"])
    h1["pass"] = h1_pass
    report["H1"] = h1
    print(f"H1 pass={h1_pass} {h1}")

    # ---------------- books ----------------
    def lbl(r):
        return f"{labels[r][0]}{labels[r][1]}"

    rep_row_of_class = {}
    for rr, (i, is_rep) in enumerate(zip(crv_class_row, crv_is_rep)):
        if is_rep:
            rep_row_of_class[i] = rr
    book_K1 = []
    for rr, (i, is_rep) in enumerate(zip(crv_class_row, crv_is_rep)):
        if not is_rep:
            par = rep_row_of_class[i]
            assert par < rr, "DUP parent must precede child"
            book_K1.append(Entry("DUP", i=par, j=rr, label="K1"))

    def star_book(members, l):
        ms = sorted(members)
        return [Entry("BRIDGE", l=l, i=ms[0], j=m, label="MST") for m in ms[1:]]

    book_K2 = star_book(g143, 5)
    book_K3 = [Entry("PATTERN", l=5, members=g143, pi="EIS", label="K3")]
    hubs_k3plus, hubs_k2 = [], []
    for l in MODULI:
        for comp in comps[l]:
            (hubs_k3plus if len(comp) >= 3 else hubs_k2).append((l, sorted(comp)))
    book_K4 = [make_hub(ms, l, tbl_A, f"K4-l{l}") for (l, ms) in hubs_k3plus]
    phi_rows = {lab[key] for key, v in phi5.items() if v and key in lab}
    Xset = phi_rows.symmetric_difference(set(g143))
    book_K5 = [Entry("TEST", l=5, members=g143, pi="EIS", nX=len(Xset),
                     label="K5")]
    report["K5_phi_size"] = len(phi_rows)
    report["K5_X_size"] = len(Xset)
    print(f"books: K1 {len(book_K1)} DUPs | K2 {len(book_K2)} bridges | "
          f"K4 {len(book_K4)} hubs (k>=3) + {len(hubs_k2)} k=2 groups | "
          f"K5 |phi|={len(phi_rows)} |X|={len(Xset)}")

    # ---------------- H2 + H3 ----------------
    print("== H2/H3 decoder + accounting fidelity ==")
    combos = [("T_cls", tbl_A, "empty", []),
              ("T_crv", tbl_B, "empty", []),
              ("T_crv", tbl_B, "K1-book", book_K1),
              ("T_cls", tbl_A, "K2-book", book_K2),
              ("T_cls", tbl_A, "K3-book", book_K3),
              ("T_cls", tbl_A, "K4-book", book_K4),
              ("T_cls", tbl_A, "K5-book", book_K5)]
    decoder_results, h2_pass, h3_pass = [], True, True
    for tname, tbl, bname, book in combos:
        codec = StreamCodec(tbl, book, EPS_BITS_MAIN)
        blob = codec.encode()
        identical = bool((codec.decode(blob) == tbl.SYM).all())
        emitted_bits = 8.0 * len(blob)
        if book:
            owned = book_delta_E_owned(tbl, book, EPS_BITS_MAIN)
            analytic = tbl.L0_data + owned["cells_delta_bits"]
        else:
            analytic = tbl.L0_data
        rel = abs(emitted_bits - analytic) / analytic
        decoder_results.append({"table": tname, "book": bname,
                                "byte_identical": identical,
                                "emitted_bytes": len(blob),
                                "analytic_bits": round(analytic, 1),
                                "rel_err": rel})
        h2_pass &= identical
        h3_pass &= rel <= 0.005
        print(f"  {tname} x {bname}: identical={identical} "
              f"emitted={emitted_bits:.0f}b analytic={analytic:.0f}b "
              f"rel={rel:.6f} ({time.time()-t0:.0f}s)")
    report["H2"] = {"pass": h2_pass, "note": "well-typed combos (reading R4)"}
    report["H3"] = {"pass": h3_pass}

    h4 = {}
    for tbl in (tbl_A, tbl_B):
        per_cell = (tbl.L0 - tbl.L0_model) / (tbl.N * C)
        xent = float(np.mean([np.dot(tbl.counts[c] / tbl.N, tbl.logq_sym[c])
                              for c in range(C)]))
        h4[tbl.name] = {"per_cell_bits": per_cell, "mean_plugin_xent": xent,
                        "diff": abs(per_cell - xent)}
    h4["pass"] = all(v["diff"] <= 0.01 for k, v in h4.items() if k != "pass")
    report["H4"] = h4
    print(f"H4 pass={h4['pass']}")

    H_pass = h1_pass and h2_pass and h3_pass and h4["pass"]
    if not H_pass:
        report["verdict"] = "VOID"
        report["void_reason"] = [k for k, ok in
                                 [("H1", h1_pass), ("H2", h2_pass),
                                  ("H3", h3_pass), ("H4", h4["pass"])] if not ok]
        with open(os.path.join(args.out, "decoder.json"), "w") as fh:
            json.dump(decoder_results, fh, indent=1)
        with open(os.path.join(args.out, "summary.json"), "w") as fh:
            json.dump(report, fh, indent=1, default=str)
        print("VERDICT: VOID —", report["void_reason"])
        return

    # ---------------- N1 ----------------
    print("== N1 foil gate ==")
    rng = np.random.default_rng(SEED)
    all_sizes = [len(comp) for l in MODULI for comp in comps[l]]
    foil_rows, n1_bad = [], defaultdict(int)
    for etype in ("BRIDGE", "PATTERN", "HUB"):
        for _ in range(FOIL_M):
            l = int(rng.choice(MODULI))
            if etype == "BRIDGE":
                i, j = sorted(int(x) for x in rng.choice(tbl_A.N, 2, replace=False))
                e, k = Entry("BRIDGE", l=l, i=i, j=j), 2
            else:
                k = int(all_sizes[int(rng.integers(0, len(all_sizes)))])
                members = sorted(int(x) for x in
                                 rng.choice(tbl_A.N, k, replace=False))
                if etype == "PATTERN":
                    e = Entry("PATTERN", l=l, members=members,
                              pi=PI_ORDER[int(rng.integers(0, 4))])
                else:
                    rh = np.full(C, -1, dtype=np.int16)
                    col_l = PRIMES.index(l)
                    for c in range(C):
                        if c != col_l:
                            rh[c] = int(rng.integers(0, l))
                    e = Entry("HUB", l=l, members=members, rh=rh)
            dE = eval_mod_entry(tbl_A, e, EPS_BITS_MAIN)["delta_E"]
            if dE < 0:
                n1_bad[etype] += 1
            foil_rows.append((etype, l, k, dE))
    with open(os.path.join(args.out, "foils.csv"), "w") as fh:
        fh.write("type,l,k,delta_E_bits\n")
        for row in foil_rows:
            fh.write(f"{row[0]},{row[1]},{row[2]},{row[3]:.3f}\n")
    n1_pass = all(n1_bad[t] == 0 for t in ("BRIDGE", "PATTERN", "HUB"))
    report["N1"] = {"pass": n1_pass, "compressing_foils": dict(n1_bad)}
    print(f"N1 pass={n1_pass} {dict(n1_bad)}")

    # ---------------- N2 ----------------
    print("== N2 shuffled-table gate ==")
    n2 = {"replayed": 0, "compressing": 0, "detail": {}}
    for l in MODULI:
        rngl = np.random.default_rng(SEED + l)
        arr = np.empty((tbl_A.N, C), dtype=object)
        for r, (_N, _c, syms, _t) in enumerate(ap_rows):
            arr[r] = syms
        for c in range(C):
            arr[:, c] = arr[rngl.permutation(tbl_A.N), c]
        tbl_sh = Table(f"T_sh{l}", [list(arr[r]) for r in range(tbl_A.N)])
        entries = [e for e in book_K4 if e.l == l]
        entries += [make_hub(ms, ll, tbl_A, "k2hub")
                    for (ll, ms) in hubs_k2 if ll == l]
        if l == 5:
            entries = book_K2 + book_K3 + book_K5 + entries
        bad = sum(1 for e in entries
                  if eval_mod_entry(tbl_sh, e, EPS_BITS_MAIN)["delta_E"] < 0)
        n2["replayed"] += len(entries)
        n2["compressing"] += bad
        n2["detail"][str(l)] = {"entries": len(entries), "compressing": bad}
    rngd = np.random.default_rng(SEED)
    arrB = np.empty((tbl_B.N, C), dtype=object)
    for r, syms in enumerate(crv_syms):
        arrB[r] = syms
    for c in range(C):
        arrB[:, c] = arrB[rngd.permutation(tbl_B.N), c]
    tbl_B_sh = Table("T_crv_sh", [list(arrB[r]) for r in range(tbl_B.N)])
    dsh = eval_dup_book(tbl_B_sh, book_K1, EPS_BITS_MAIN)
    dup_bad = int((dsh["per_entry_delta_E"] < 0).sum())
    n2["replayed"] += len(book_K1)
    n2["compressing"] += dup_bad
    n2["detail"]["DUP"] = {"entries": len(book_K1), "compressing": dup_bad}
    n2["pass"] = n2["compressing"] == 0
    report["N2"] = n2
    print(f"N2 pass={n2['pass']} {n2['detail']}")

    # ---------------- N3 ----------------
    print("== N3 relabelling invariance ==")
    rngp = np.random.default_rng(SEED)
    perm = rngp.permutation(tbl_A.N)
    inv = np.empty_like(perm)
    inv[perm] = np.arange(tbl_A.N)
    tbl_Ap = Table("T_cls_perm", [ap_rows[perm[r]][2] for r in range(tbl_A.N)])

    def relabel(book):
        out = []
        for e in book:
            if e.etype in ("DUP", "BRIDGE"):
                out.append(Entry(e.etype, l=e.l, i=int(inv[e.i]), j=int(inv[e.j])))
            else:
                out.append(Entry(e.etype, l=e.l,
                                 members=sorted(int(inv[m]) for m in e.members),
                                 pi=e.pi, rh=e.rh, nX=e.nX))
        return out

    n3 = {"L0_diff": abs(tbl_A.L0 - tbl_Ap.L0), "analytic": {}, "streams": {}}
    for name, book in (("K2", book_K2), ("K3", book_K3),
                       ("K4", book_K4), ("K5", book_K5)):
        a = book_delta_E_owned(tbl_A, book, EPS_BITS_MAIN)["delta_E_bits"]
        b = book_delta_E_owned(tbl_Ap, relabel(book), EPS_BITS_MAIN)["delta_E_bits"]
        n3["analytic"][name] = abs(a - b)
    for name, book in (("empty-A", []), ("K3", book_K3), ("K4", book_K4)):
        b1 = StreamCodec(tbl_A, book, EPS_BITS_MAIN).encode()
        b2 = StreamCodec(tbl_Ap, relabel(book), EPS_BITS_MAIN).encode()
        n3["streams"][name] = abs(len(b1) - len(b2))
    n3["pass"] = (n3["L0_diff"] == 0.0
                  and all(v <= 1e-6 for v in n3["analytic"].values())
                  and all(v <= 1 for v in n3["streams"].values()))
    report["N3"] = n3
    print(f"N3 pass={n3['pass']} {n3}")

    N_pass = n1_pass and n2["pass"] and n3["pass"]
    if not N_pass:
        report["verdict"] = ("GAMEABLE" if not n1_pass else
                             "NULL-LEAK" if not n2["pass"] else "INVARIANCE-FAIL")
        with open(os.path.join(args.out, "summary.json"), "w") as fh:
            json.dump(report, fh, indent=1, default=str)
        print("VERDICT:", report["verdict"])
        return

    # ---------------- K cells ----------------
    print("== K1 DUP ==")
    d1 = eval_dup_book(tbl_B, book_K1, EPS_BITS_MAIN)
    dup_rows = [e.j for e in book_K1]
    base_bits = float(tbl_B.logq_cell[dup_rows].sum())
    ratio = d1["cells_bits"] / base_bits
    K = {"K1": {"delta_E": d1["delta_E"], "dup_cell_ratio": ratio,
                "n_dups": len(book_K1), "escapes": d1["escapes"],
                "pass": d1["delta_E"] < 0 and ratio <= 0.05}}
    print(K["K1"])

    def k_eval(eps_bits):
        fb = flag_bits(eps_bits)
        out = {}
        # K2 — MST of bridges over the 143-family
        evs = [eval_mod_entry(tbl_A, e, eps_bits) for e in book_K2]
        dE2 = sum(ev["delta_E"] for ev in evs)
        preds, meas = [], []
        for e, ev in zip(book_K2, evs):
            rh = e.residue_row(tbl_A)
            cols = np.nonzero(rh >= 0)[0]
            preds.append(float(tbl_A.negLogQ[5][cols, rh[cols]].sum()) - fb * 24)
            meas.append(ev["S"])
        mean_pred, mean_meas = float(np.mean(preds)), float(np.mean(meas))
        rel2 = abs(mean_meas - mean_pred) / abs(mean_pred)
        S2 = sum(meas)
        cost2 = sum(e.cost(tbl_A) for e in book_K2)
        out["K2"] = {"delta_E": dE2, "mean_saving_pred": mean_pred,
                     "mean_saving_meas": mean_meas, "rel_err": rel2,
                     "escapes": sum(ev["escapes"] for ev in evs),
                     "pass": dE2 < 0 and rel2 <= 0.15}
        # K3 — PATTERN vs MST, S-based per P-b (prereg formula verbatim)
        ev3 = eval_mod_entry(tbl_A, book_K3[0], eps_bits)
        k = len(g143)
        A_meas = (book_K3[0].cost(tbl_A) - ev3["S"]) - (cost2 - S2)
        A_pred = 2 + gamma_bits(k) + (2 - k) * tbl_A.log2N - 5 * (k - 1)
        rel3 = abs(A_meas - A_pred) / abs(A_pred)
        out["K3"] = {"A_pat_meas": A_meas, "A_pat_pred": A_pred, "rel_err": rel3,
                     "delta_E_true": ev3["delta_E"],
                     "pass": rel3 <= 0.20 and A_meas < 0}
        # K4 — HUB vs MST for every k>=3 group, S-based per P-c
        rows4, n_neg, n_band = [], 0, 0
        for (l, ms) in hubs_k3plus:
            k4 = len(ms)
            hub = make_hub(ms, l, tbl_A)
            evh = eval_mod_entry(tbl_A, hub, eps_bits)
            mst = star_book(ms, l)
            evm = [eval_mod_entry(tbl_A, e, eps_bits) for e in mst]
            A = ((hub.cost(tbl_A) - evh["S"])
                 - (sum(e.cost(tbl_A) for e in mst) - sum(x["S"] for x in evm)))
            Ap = gamma_bits(k4) + (2 - k4) * tbl_A.log2N - 5 * (k4 - 1)
            neg = A < 0
            band = abs(A - Ap) / abs(Ap) <= 0.20
            n_neg += neg
            n_band += band
            rows4.append((l, k4, A, Ap, int(neg), int(band)))
        ngr = len(hubs_k3plus)
        A2list = []
        for (l, ms) in hubs_k2:
            hub = make_hub(ms, l, tbl_A)
            evh = eval_mod_entry(tbl_A, hub, eps_bits)
            mst = star_book(ms, l)
            evm = [eval_mod_entry(tbl_A, e, eps_bits) for e in mst]
            A2list.append((hub.cost(tbl_A) - evh["S"])
                          - (sum(e.cost(tbl_A) for e in mst)
                             - sum(x["S"] for x in evm)))
        out["K4"] = {"groups": ngr, "frac_negative": n_neg / ngr,
                     "frac_in_band": n_band / ngr,
                     "k2_groups": len(A2list),
                     "k2_mean_A": float(np.mean(A2list)),
                     "k2_frac_negative": float(np.mean([a < 0 for a in A2list])),
                     "pass": n_neg / ngr >= 0.99 and n_band / ngr >= 0.95,
                     "rows": rows4}
        # K5 — TEST vs PATTERN (identical claims; description-cost difference)
        dE_test = book_K5[0].cost(tbl_A) + (ev3["delta_E"] - book_K3[0].cost(tbl_A))
        lhs = dE_test < ev3["delta_E"]
        rhs = len(Xset) < len(g143) - math.log2(PHI_SIZE) / tbl_A.log2N
        out["K5"] = {"X": len(Xset), "delta_E_test": dE_test,
                     "delta_E_pattern": ev3["delta_E"],
                     "test_beats_pattern": lhs, "iff_rhs": rhs,
                     "pass": lhs == rhs}
        return out

    print("== K2-K5 (eps0 = 2^-6) ==")
    Kmain = k_eval(EPS_BITS_MAIN)
    for name in ("K2", "K3", "K4", "K5"):
        K[name] = Kmain[name]
        print(f"{name}: " + str({k: v for k, v in Kmain[name].items()
                                 if k != "rows"}))

    print("== K6 eps sweep ==")
    k6 = {"flips": []}
    for eb in EPS_BITS_SWEEP:
        if eb == EPS_BITS_MAIN:
            continue
        Ke = k_eval(eb)
        for name in ("K2", "K3", "K4", "K5"):
            if Ke[name]["pass"] != Kmain[name]["pass"]:
                k6["flips"].append({"cell": name, "eps_bits": eb})
    k6["pass"] = len(k6["flips"]) == 0
    K["K6"] = k6
    print(f"K6 pass={k6['pass']} flips={k6['flips']}")

    print("== K7 headroom (report-only) ==")
    full_book = []
    for l in MODULI:
        for comp in comps[l]:
            ms = sorted(comp)
            mst = star_book(ms, l)
            d_mst = sum(eval_mod_entry(tbl_A, e, EPS_BITS_MAIN)["delta_E"]
                        for e in mst)
            hub = make_hub(ms, l, tbl_A)
            d_hub = eval_mod_entry(tbl_A, hub, EPS_BITS_MAIN)["delta_E"]
            best = ("MST", d_mst, mst)
            if d_hub < best[1]:
                best = ("HUB", d_hub, [hub])
            for pi in PI_ORDER:
                pe = Entry("PATTERN", l=l, members=ms, pi=pi)
                d = eval_mod_entry(tbl_A, pe, EPS_BITS_MAIN)["delta_E"]
                if d < best[1]:
                    best = (f"PATTERN-{pi}", d, [pe])
            full_book.extend(best[2])

    def sort_key(e):
        if e.members is not None:
            return (e.l, -len(e.members), e.members[0])
        return (e.l, -2, e.j)

    full_book.sort(key=sort_key)
    res7 = book_delta_E_owned(tbl_A, full_book, EPS_BITS_MAIN)
    K["K7"] = {"entries": len(full_book), "delta_E": res7["delta_E_bits"],
               "fraction_of_L0": -res7["delta_E_bits"] / tbl_A.L0}
    print(K["K7"])

    # ---------------- verdict (§9) ----------------
    if not K["K1"]["pass"]:
        verdict = "UNDERPOWERED"
    elif not all(K[n]["pass"] for n in ("K2", "K3", "K4", "K5")):
        verdict = "PARTIAL"
    elif not k6["pass"]:
        verdict = "SENSITIVE"
    else:
        verdict = "VALIDATED"
    report["K"] = {n: {kk: vv for kk, vv in K[n].items() if kk != "rows"}
                   for n in K}
    report["verdict"] = verdict
    report["failed_K_cells"] = [n for n in ("K2", "K3", "K4", "K5")
                                if not K[n]["pass"]]
    report["runtime_s"] = round(time.time() - t0, 1)

    # ---------------- outputs (§8) ----------------
    def members_hash(labels_list):
        return "sha256:" + hashlib.sha256(
            ",".join(sorted(labels_list)).encode()).hexdigest()[:16]

    with open(os.path.join(args.out, "k4_groups.csv"), "w") as fh:
        fh.write("l,k,A_hub_meas,A_hub_pred,negative,in_band\n")
        for (l, k4, A, Ap, s, b) in Kmain["K4"]["rows"]:
            fh.write(f"{l},{k4},{A:.2f},{Ap:.2f},{s},{b}\n")
    with open(os.path.join(args.out, "decoder.json"), "w") as fh:
        json.dump(decoder_results, fh, indent=1)
    with open(os.path.join(args.out, "baseline.json"), "w") as fh:
        json.dump({"A": {"L0": tbl_A.L0, "L0_model": tbl_A.L0_model,
                         "L0_data": tbl_A.L0_data, "nsym": tbl_A.nsym},
                   "B": {"L0": tbl_B.L0, "L0_model": tbl_B.L0_model,
                         "L0_data": tbl_B.L0_data}}, fh, indent=1)
    with open(os.path.join(args.out, "headroom.json"), "w") as fh:
        json.dump(K["K7"], fh, indent=1)
    with open(os.path.join(args.out, "energy.jsonl"), "w") as fh:
        recs = [
            {"cell": "K1", "eps0": -EPS_BITS_MAIN, "entry_type": "DUP",
             "k": len(book_K1),
             "members_hash": members_hash([crv_labels[e.j] for e in book_K1]),
             "delta_E_bits": K["K1"]["delta_E"],
             "escapes": K["K1"]["escapes"],
             "gate": "PASS" if K["K1"]["pass"] else "FAIL"},
            {"cell": "K2", "eps0": -EPS_BITS_MAIN, "entry_type": "BRIDGE",
             "l": 5, "k": len(g143), "members_hash": members_hash(
                 [lbl(m) for m in g143]),
             "delta_E_bits": K["K2"]["delta_E"],
             "saving_pred_mean": K["K2"]["mean_saving_pred"],
             "saving_meas_mean": K["K2"]["mean_saving_meas"],
             "rel_err": K["K2"]["rel_err"], "escapes": K["K2"]["escapes"],
             "gate": "PASS" if K["K2"]["pass"] else "FAIL"},
            {"cell": "K3", "eps0": -EPS_BITS_MAIN, "entry_type": "PATTERN",
             "l": 5, "k": len(g143),
             "members_hash": members_hash([lbl(m) for m in g143]),
             "delta_E_bits": K["K3"]["delta_E_true"],
             "A_pat_meas": K["K3"]["A_pat_meas"],
             "A_pat_pred": K["K3"]["A_pat_pred"],
             "rel_err": K["K3"]["rel_err"],
             "gate": "PASS" if K["K3"]["pass"] else "FAIL"},
            {"cell": "K4", "eps0": -EPS_BITS_MAIN, "entry_type": "HUB",
             "groups": K["K4"]["groups"],
             "frac_negative": K["K4"]["frac_negative"],
             "frac_in_band": K["K4"]["frac_in_band"],
             "gate": "PASS" if K["K4"]["pass"] else "FAIL"},
            {"cell": "K5", "eps0": -EPS_BITS_MAIN, "entry_type": "TEST",
             "l": 5, "X": K["K5"]["X"],
             "delta_E_test": K["K5"]["delta_E_test"],
             "delta_E_pattern": K["K5"]["delta_E_pattern"],
             "gate": "PASS" if K["K5"]["pass"] else "FAIL"},
            {"cell": "K7", "eps0": -EPS_BITS_MAIN,
             "entries": K["K7"]["entries"], "delta_E_bits": K["K7"]["delta_E"],
             "fraction_of_L0": K["K7"]["fraction_of_L0"], "gate": "REPORT"},
        ]
        for rec in recs:
            fh.write(json.dumps(rec, default=float) + "\n")
    with open(os.path.join(args.out, "summary.json"), "w") as fh:
        json.dump(report, fh, indent=1, default=str)

    print(f"\n== VERDICT: {verdict} ==  ({report['runtime_s']}s)")


if __name__ == "__main__":
    main()
