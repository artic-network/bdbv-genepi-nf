#!/usr/bin/env python3
"""Download aligned FASTA and metadata from Pathoplexus LAPIS API."""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from urllib.parse import urlencode, parse_qsl
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError

import pandas as pd

BASE_URL = "https://lapis.pathoplexus.org"

COUNTRY_MAP = {
    "drc": "Democratic Republic of the Congo",
    "kenya": "Kenya",
    "nigeria": "Nigeria",
    "uganda": "Uganda",
}

OUTBREAK_COUNTRIES = ["drc", "kenya", "nigeria", "uganda"]

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


def api_error_detail(e):
    """The message LAPIS returns with an HTTP error (for a bad filter it names the problem)."""
    try:
        detail = json.loads(e.read().decode("utf-8"))["error"]["detail"]
    except Exception:
        return ""
    # for an unknown filter LAPIS lists every valid key: keep the useful first part
    return detail.split(" Valid keys are:")[0]


def fetch_url(url, timeout=120):
    req = Request(url, headers={"Accept": "*/*"})
    try:
        with urlopen(req, timeout=timeout) as resp:
            return resp.read().decode("utf-8")
    except HTTPError as e:
        print(f"Error fetching {url}: {e}", file=sys.stderr)
        detail = api_error_detail(e)
        if detail:
            print(f"  Pathoplexus said: {detail}", file=sys.stderr)
            print("  Check --filters: field names are case-sensitive metadata fields, "
                  "see https://lapis.pathoplexus.org/ebola-bdbv/api-docs", file=sys.stderr)
        sys.exit(1)
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


# Filters that would change the output format or the country handling, so they are not accepted in --filters
RESERVED_FILTERS = {
    "dataFormat": "the output format is set by this script",
    "fields": "the metadata fields are set by this script",
    "limit": "all matching sequences are downloaded",
    "offset": "all matching sequences are downloaded",
    "orderBy": "ordering is not used",
    "downloadAsFile": "the output is handled by this script",
    "compression": "the output is handled by this script",
    "geoLocCountry": "use --scope to choose countries",
}


def parse_filters(text):
    """Parse extra Pathoplexus (LAPIS) filters written like a URL query string.

    'geoLocAdmin1=Ituri&hostNameScientific=Homo sapiens&lengthFrom=15000'. Separate filters with & (or
    a new line). A repeated key means either value. Values are percent-decoded and + is a space, so
    the query string of a Pathoplexus URL can be pasted. Returns a dict of key -> list of values.
    """
    filters = {}
    if not text:
        return filters
    for chunk in text.replace("\r", "\n").replace("\n", "&").split("&"):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "=" not in chunk:
            print(f"ERROR: could not read filter '{chunk}': expected key=value, "
                  f"e.g. geoLocAdmin1=Ituri", file=sys.stderr)
            sys.exit(1)
        for key, value in parse_qsl(chunk, keep_blank_values=True):
            key = key.strip()
            value = value.strip()
            if not key or value == "":
                print(f"ERROR: filter '{chunk}' needs both a name and a value", file=sys.stderr)
                sys.exit(1)
            if key in RESERVED_FILTERS:
                print(f"ERROR: filter '{key}' is not allowed in --filters ({RESERVED_FILTERS[key]})",
                      file=sys.stderr)
                sys.exit(1)
            filters.setdefault(key, []).append(value)
    return filters


def merge_filters(base, extra):
    """Add the user's filters to the built-in ones. A user filter replaces the built-in one of the
    same name, and any groupName filter replaces the default group restriction."""
    merged = dict(base)
    if any(k == "groupName" or k.startswith("groupName.") for k in extra):
        merged.pop("groupName", None)
    merged.update(extra)
    return merged


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
        "--min-date", default="2025-01-01",
        help="Earliest collection date to include (default: 2025-01-01)"
    )
    parser.add_argument(
        "--max-date", default=None,
        help="Latest collection date to include (default: no upper limit)"
    )
    parser.add_argument(
        "--filters", default=None,
        help="Extra Pathoplexus filters as a URL query string, e.g. "
             "'geoLocAdmin1=Ituri&lengthFrom=15000'. Applied to the metadata and the sequences. "
             "A filter replaces the built-in one of the same name"
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
    extra_filters = parse_filters(args.filters)
    if extra_filters:
        print("Extra filters: " + "&".join(f"{k}={v}" for k, vs in extra_filters.items() for v in vs))

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
        filters = {"geoLocCountry": country_name,
                   "sampleCollectionDateRangeLowerFrom": args.min_date}
        if args.max_date:
            filters["sampleCollectionDateRangeUpperTo"] = args.max_date
        if countries == ["drc"] and args.group:
            filters["groupName"] = args.group
        filters = merge_filters(filters, extra_filters)

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
