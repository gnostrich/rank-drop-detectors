#!/usr/bin/env python3
"""Probe 05 data acquisition — Artin representations with SHA256 manifest.

Dim-2, rational-character (CharacterField=1) Artin reps, conductor <= 2000,
a_p via lfunc_instances -> lfunc_lfunctions.euler_factors (route verified on
2.23.3t2.b). Reps without an L-function row are recorded with null factors
and counted. Sealed rows (conductor 1801-2000) stored flagged, never scanned.

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
HDRS = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) probe-05-downloader",
        "Cookie": "human=1"}
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
        time.sleep(0.35)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    lte = urllib.parse.quote('py{"$lte":2000}')
    reps = paged("artin_reps",
                 f"Dim=2&CharacterField=1&Conductor={lte}",
                 "Baselabel,Conductor,BadPrimes,HardPrimes,Is_Even,"
                 "Indicator,GaloisConjugates")
    reps.sort(key=lambda r: (int(r["Conductor"]), r["Baselabel"]))
    print(f"artin reps enumerated: {len(reps)}")

    n_missing = 0
    with open(os.path.join(args.out, "artin.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Baselabel", "Conductor", "BadPrimes", "HardPrimes",
                    "Is_Even", "Indicator", "LocalFactors", "euler_factors",
                    "sealed"])
        for i, r in enumerate(reps):
            url_id = urllib.parse.quote(f"ArtinRepresentation/"
                                        f"{r['Baselabel']}.a")
            d = get(f"{BASE}/lfunc_instances/?url=s{url_id}"
                    f"&_format=json&_fields=label")
            ef = None
            if d.get("data"):
                lbl = d["data"][0]["label"]
                d2 = get(f"{BASE}/lfunc_lfunctions/?label="
                         f"s{urllib.parse.quote(lbl)}"
                         f"&_format=json&_fields=euler_factors")
                if d2.get("data"):
                    ef = d2["data"][0].get("euler_factors")
            if ef is None:
                n_missing += 1
            gc = r.get("GaloisConjugates") or []
            lf = gc[0].get("LocalFactors") if gc else None
            w.writerow([r["Baselabel"], r["Conductor"],
                        json.dumps(r.get("BadPrimes")),
                        json.dumps(r.get("HardPrimes")),
                        r.get("Is_Even"), r.get("Indicator"),
                        json.dumps(lf), json.dumps(ef),
                        int(int(r["Conductor"]) > 1800)])
            if i % 100 == 0:
                print(f"  lfunc {i}/{len(reps)} (missing so far {n_missing})")
            time.sleep(0.2)
    print(f"artin reps: {len(reps)}, without L-function: {n_missing}")

    with open(os.path.join(args.out, "manifest.json"), "w") as fh:
        json.dump({"endpoint": BASE, "requests": manifest}, fh, indent=1)
    print(f"manifest: {len(manifest)} requests")


if __name__ == "__main__":
    main()
