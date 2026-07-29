#!/usr/bin/env python3
"""Probe 04 data acquisition — LMFDB snapshots with SHA256 manifest.

Families per PREREG §1: g2c curves cond<=10000 (+ per-class L euler_factors),
smf_newforms weight [2,0] (any family). Sealed rows (cond 9001-10000) are
downloaded, stored, flagged, never scanned. Nothing here scans or measures.

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
HDRS = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) probe-04-downloader",
        "Cookie": "human=1"}
LTE = urllib.parse.quote('py{"$lte":10000}')
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
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"download failed: {url}: {last}")


def paged(table, query, fields, max_pages=800):
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
        time.sleep(0.35)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    # F5 — genus-2 curves
    rows = paged("g2c_curves", f"cond={LTE}",
                 "label,class,cond,abs_disc,eqn,bad_lfactors,is_gl2_type,"
                 "end_alg,geom_end_alg,st_group,Lhash")
    rows.sort(key=lambda r: (r["cond"], r["label"]))
    with open(os.path.join(args.out, "f5_g2c.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["label", "class", "cond", "abs_disc", "eqn",
                    "bad_lfactors", "is_gl2_type", "end_alg", "geom_end_alg",
                    "st_group", "Lhash", "sealed"])
        for r in rows:
            w.writerow([r["label"], r["class"], r["cond"], r["abs_disc"],
                        json.dumps(r["eqn"]), json.dumps(r["bad_lfactors"]),
                        r["is_gl2_type"], r["end_alg"], r["geom_end_alg"],
                        r["st_group"], r["Lhash"], int(r["cond"] > 9000)])
    classes = sorted({(r["cond"], r["class"]) for r in rows})
    print(f"F5 curves: {len(rows)}; classes: {len(classes)}")

    # F6 — class L euler factors, one lfunc request per class
    with open(os.path.join(args.out, "f6_class_euler.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cond", "class", "euler_factors", "sealed"])
        for i, (cond, cls) in enumerate(classes):
            origin = urllib.parse.quote(f"Genus2Curve/Q/{cond}/{cls}")
            d = get(f"{BASE}/lfunc_lfunctions/?origin=s{origin}"
                    f"&_format=json&_fields=origin,euler_factors")
            ef = d["data"][0]["euler_factors"] if d.get("data") else None
            w.writerow([cond, cls, json.dumps(ef), int(cond > 9000)])
            if i % 200 == 0:
                print(f"  euler {i}/{len(classes)}")
            time.sleep(0.25)
    print(f"F6 classes: {len(classes)}")

    # F7 — weight-[2,0] Siegel/paramodular newforms, any family
    w20 = urllib.parse.quote("py[2,0]")
    rows7 = paged("smf_newforms", f"degree=2&weight={w20}",
                  "label,level,weight,family,dim,aut_rep_type,"
                  "trace_lambda_p,traces,related_objects")
    with open(os.path.join(args.out, "f7_smf.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["label", "level", "weight", "family", "dim",
                    "aut_rep_type", "trace_lambda_p", "related_objects"])
        for r in rows7:
            w.writerow([r["label"], r["level"], json.dumps(r["weight"]),
                        r["family"], r["dim"], r["aut_rep_type"],
                        json.dumps(r.get("trace_lambda_p")),
                        json.dumps(r.get("related_objects"))])
    print(f"F7 weight-[2,0] forms: {len(rows7)}")

    with open(os.path.join(args.out, "manifest.json"), "w") as fh:
        json.dump({"endpoint": BASE, "requests": manifest}, fh, indent=1)
    print(f"manifest: {len(manifest)} requests")


if __name__ == "__main__":
    main()
