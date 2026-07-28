#!/usr/bin/env python3
"""Probe 03 data acquisition — LMFDB snapshots with SHA256 manifest.

Downloads exactly the families and bands declared in PREREG.md §1 from the
beta.lmfdb.org API (gate cookie human=1; the API is unversioned — the pin is
the committed snapshot plus this manifest). Sealed rows (top decile, 901-1000)
are downloaded, stored, and flagged 'sealed'; nothing here scans or measures.

Usage: python3 download.py --out data/
"""

import argparse
import csv
import hashlib
import json
import os
import time
import urllib.parse
import urllib.request

BASE = "https://beta.lmfdb.org/api"
HDRS = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) probe-03-downloader",
        "Cookie": "human=1"}
PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
          67, 71, 73, 79, 83, 89, 97]
FIELDS = ["2.2.5.1", "2.2.8.1"]
LTE1000 = urllib.parse.quote('py{"$lte":1000}')

manifest = []


def get(url):
    last = None
    for attempt in range(6):
        try:
            req = urllib.request.Request(url, headers=HDRS)
            b = urllib.request.urlopen(req, timeout=120).read()
            d = json.loads(b)
            manifest.append({"url": url,
                             "sha256": hashlib.sha256(b).hexdigest(),
                             "ts": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                                 time.gmtime()),
                             "bytes": len(b)})
            return d
        except Exception as e:  # noqa: BLE001 — retry then surface
            last = e
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"download failed: {url}: {last}")


def paged(table, query, fields, max_pages=400):
    url = f"{BASE}/{table}/?{query}&_format=json&_fields={fields}"
    out = []
    for _ in range(max_pages):
        d = get(url)
        rows = d.get("data", [])
        out.extend(rows)
        nxt = d.get("next", "")
        if len(rows) < 100 or not nxt or "_offset" not in nxt:
            break
        url = "https://beta.lmfdb.org" + nxt if nxt.startswith("/") else nxt
        time.sleep(0.4)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    # F2 — newforms weight 2, trivial character, level <= 1000, all dims
    rows = paged("mf_newforms", f"weight=2&char_order=1&level={LTE1000}",
                 "label,level,dim,traces")
    with open(os.path.join(args.out, "f2_newforms.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["label", "level", "dim", "sealed"] +
                   [f"a{p}" for p in PRIMES])
        for r in sorted(rows, key=lambda r: (r["level"], r["label"])):
            tr = r["traces"]
            w.writerow([r["label"], r["level"], r["dim"],
                        int(r["level"] > 900)] +
                       [tr[p - 1] if tr and len(tr) >= p else "" for p in PRIMES])
    print(f"F2 newforms: {len(rows)}")

    # F1 — EC/Q curves, conductor <= 1000 (class reps chosen in run.py)
    rows = paged("ec_curvedata", f"conductor={LTE1000}",
                 "lmfdb_label,Clabel,conductor,ainvs")
    with open(os.path.join(args.out, "f1_ec_q.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["lmfdb_label", "Clabel", "conductor", "ainvs"])
        for r in sorted(rows, key=lambda r: (r["conductor"], r["lmfdb_label"])):
            w.writerow([r["lmfdb_label"], r["Clabel"], r["conductor"],
                        json.dumps(r["ainvs"])])
    print(f"F1 EC/Q curves: {len(rows)}")

    # F3 — EC over the two real quadratic fields
    for tag, fl in zip(("f3a", "f3b"), FIELDS):
        rows = paged("ec_nfcurves",
                     f"field_label={fl}&conductor_norm={LTE1000}",
                     "label,short_label,conductor_norm,conductor_label,"
                     "iso_label,number,ainvs,bad_primes,non_min_p,base_change")
        with open(os.path.join(args.out, f"{tag}_ecnf.csv"), "w",
                  newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["label", "short_label", "conductor_norm",
                        "conductor_label", "iso_label", "number", "ainvs",
                        "bad_primes", "non_min_p", "base_change", "sealed"])
            for r in sorted(rows, key=lambda r: (r["conductor_norm"],
                                                 r["label"])):
                w.writerow([r["label"], r["short_label"], r["conductor_norm"],
                            r["conductor_label"], r["iso_label"], r["number"],
                            r["ainvs"], json.dumps(r["bad_primes"]),
                            json.dumps(r["non_min_p"]),
                            json.dumps(r["base_change"]),
                            int(r["conductor_norm"] > 900)])
        print(f"{tag} EC/{fl}: {len(rows)}")

    # F4 — HMFs, parallel weight 2, dimension 1; eigenvalues per label
    for tag, fl in zip(("f4a", "f4b"), FIELDS):
        forms = paged("hmf_forms",
                      f"field_label={fl}&level_norm={LTE1000}&dimension=1",
                      "label,level_norm,level_label,dimension,weight,"
                      "is_base_change,is_CM")
        forms = [r for r in forms if str(r.get("weight")) in ("[2, 2]", "[2,2]")]
        out_rows = []
        for i, r in enumerate(forms):
            d = get(f"{BASE}/hmf_hecke/?label=s{urllib.parse.quote(r['label'])}"
                    f"&_format=json&_fields=label,hecke_eigenvalues")
            ev = d["data"][0]["hecke_eigenvalues"] if d.get("data") else []
            out_rows.append((r, ev[:25]))
            if i % 100 == 0:
                print(f"  {tag} eigenvalues {i}/{len(forms)}")
            time.sleep(0.25)
        with open(os.path.join(args.out, f"{tag}_hmf.csv"), "w",
                  newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["label", "level_norm", "level_label", "dimension",
                        "weight", "is_base_change", "is_CM", "sealed",
                        "eigs25"])
            for r, ev in sorted(out_rows, key=lambda t: (t[0]["level_norm"],
                                                         t[0]["label"])):
                w.writerow([r["label"], r["level_norm"], r["level_label"],
                            r["dimension"], r["weight"], r["is_base_change"],
                            r["is_CM"], int(r["level_norm"] > 900),
                            json.dumps(ev)])
        print(f"{tag} HMF/{fl}: {len(forms)}")

    # prime-ideal orderings
    fields_out = {}
    for fl in FIELDS:
        d = get(f"{BASE}/hmf_fields/?label={fl}&_format=json"
                f"&_fields=label,primes")
        fields_out[fl] = d["data"][0]["primes"][:40]
    with open(os.path.join(args.out, "hmf_fields_primes.json"), "w") as fh:
        json.dump(fields_out, fh, indent=1)

    with open(os.path.join(args.out, "manifest.json"), "w") as fh:
        json.dump({"endpoint": BASE, "requests": manifest}, fh, indent=1)
    print(f"manifest: {len(manifest)} requests")


if __name__ == "__main__":
    main()
