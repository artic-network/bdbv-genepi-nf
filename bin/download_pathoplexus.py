#!/usr/bin/env python3
"""Download aligned FASTA and metadata from Pathoplexus LAPIS API."""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen, Request
from urllib.error import URLError

import pandas as pd

BASE_URL = "https://lapis.pathoplexus.org"

COUNTRY_MAP = {
    "drc": "Democratic Republic of the Congo",
    "kenya": "Kenya",
    "uganda": "Uganda",
}

OUTBREAK_COUNTRIES = ["drc", "kenya", "uganda"]

METADATA_FIELDS = [
    "accessionVersion", "outbreak", "sampleCollectionDate",
    "geoLocCountry", "geoLocAdmin1", "geoLocAdmin2", "geoLocCity",
    "authors", "authorAffiliations", "hostNameScientific",
    "length", "earliestReleaseDate", "submissionId", "groupName",
    "dataUseTerms", "dataUseTermsUrl", "biosampleAccession",
    "completeness", "insdcAccessionFull", "insdcRawReadsAccession",
    "isLabHost", "sampleReceivedDate", "specimenCollectorSampleId",
    "totalAmbiguousNucs", "totalUnknownNucs",
]


def fetch_url(url, timeout=120):
    req = Request(url, headers={"Accept": "*/*"})
    try:
        with urlopen(req, timeout=timeout) as resp:
            return resp.read().decode("utf-8")
    except URLError as e:
        print(f"Error fetching {url}: {e}", file=sys.stderr)
        sys.exit(1)


def build_query(params):
    return urlencode(params, doseq=True)


def download_metadata(organism, outpath, filters=None):
    fields = ",".join(METADATA_FIELDS)
    params = {"dataFormat": "tsv", "versionStatus": "LATEST_VERSION", "fields": fields}
    if filters:
        params.update(filters)
    url = f"{BASE_URL}/{organism}/sample/details?{build_query(params)}"
    print(f"Downloading metadata from {url}")
    tsv = fetch_url(url)
    outpath.write_text(tsv)
    n_rows = tsv.count("\n") - 1
    print(f"  → {n_rows} sequences, saved to {outpath}")
    return n_rows


def download_sequences(organism, outpath, filters=None):
    params = {"versionStatus": "LATEST_VERSION"}
    if filters:
        params.update(filters)
    url = f"{BASE_URL}/{organism}/sample/alignedNucleotideSequences?{build_query(params)}"
    print(f"Downloading aligned sequences from {url}")
    fasta = fetch_url(url, timeout=300)
    outpath.write_text(fasta)
    n_seqs = fasta.count("\n>") + (1 if fasta.startswith(">") else 0)
    print(f"  → {n_seqs} sequences, saved to {outpath}")
    return n_seqs


def parse_scope(scope_str):
    tokens = [t.strip().lower() for t in scope_str.split(",")]
    if "outbreak" in tokens:
        if len(tokens) > 1:
            print("ERROR: --scope 'outbreak' cannot be combined with country names. "
                  "Use either 'outbreak' or a comma-separated list of countries (e.g. 'drc,kenya').",
                  file=sys.stderr)
            sys.exit(1)
        return list(OUTBREAK_COUNTRIES)
    for t in tokens:
        if t not in COUNTRY_MAP:
            print(f"ERROR: unknown country '{t}'. Known: {', '.join(COUNTRY_MAP.keys())}",
                  file=sys.stderr)
            sys.exit(1)
    return tokens


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o", "--outdir", type=Path, default=None,
        help="Output directory (default: data/YYYYMMDD)"
    )
    parser.add_argument(
        "--organism", default="ebola-bdbv",
        help="LAPIS organism slug (default: ebola-bdbv)"
    )
    parser.add_argument(
        "--scope", default="drc",
        help="Comma-separated countries (drc,kenya,uganda) or 'outbreak' for all (default: drc)"
    )
    parser.add_argument(
        "--group", default="Pathogen Genomic Lab",
        help="Filter to this submitting group (only applied for single-country DRC scope)"
    )
    parser.add_argument(
        "--metadata-only", action="store_true",
        help="Download only metadata, skip sequences"
    )
    parser.add_argument(
        "--sequences-only", action="store_true",
        help="Download only sequences, skip metadata"
    )
    args = parser.parse_args()

    countries = parse_scope(args.scope)

    now = datetime.now()
    timestamp = now.strftime("%Y-%m-%dT%H%M")

    if args.outdir is None:
        args.outdir = Path("data") / now.strftime("%Y%m%d")

    args.outdir.mkdir(parents=True, exist_ok=True)

    all_meta = []
    all_fasta = []

    for country_key in countries:
        country_name = COUNTRY_MAP[country_key]
        print(f"\n--- {country_name} ---")
        filters = {"geoLocCountry": country_name}
        if countries == ["drc"] and args.group:
            filters["groupName"] = args.group

        if not args.sequences_only:
            meta_path = args.outdir / f"{country_key}_metadata.tsv"
            download_metadata(args.organism, meta_path, filters)
            all_meta.append(meta_path)

        if not args.metadata_only:
            fasta_path = args.outdir / f"{country_key}_aligned.fasta"
            download_sequences(args.organism, fasta_path, filters)
            all_fasta.append(fasta_path)

    merged_meta = args.outdir / "metadata.tsv"
    merged_fasta = args.outdir / "aligned.fasta"

    if not args.sequences_only:
        if len(all_meta) == 1:
            all_meta[0].rename(merged_meta)
        else:
            dfs = [pd.read_csv(p, sep="\t") for p in all_meta]
            merged = pd.concat(dfs, ignore_index=True)
            merged.to_csv(merged_meta, sep="\t", index=False)
            print(f"\nMerged metadata: {len(merged)} rows → {merged_meta}")

    if not args.metadata_only:
        if len(all_fasta) == 1:
            all_fasta[0].rename(merged_fasta)
        else:
            with open(merged_fasta, "w") as out:
                for fp in all_fasta:
                    out.write(fp.read_text())
            n_total = sum(1 for line in open(merged_fasta) if line.startswith(">"))
            print(f"\nMerged sequences: {n_total} → {merged_fasta}")

    paths = {}
    if not args.sequences_only:
        paths["metadata"] = str(merged_meta)
    if not args.metadata_only:
        paths["fasta"] = str(merged_fasta)
    print(json.dumps(paths))


if __name__ == "__main__":
    main()
