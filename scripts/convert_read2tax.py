#!/usr/bin/env python3
"""Convert common classifier outputs to the ring-trial read2tax TSV schema."""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path

STATUSES = {"classified", "unclassified", "excluded_preclassification", "error"}


def open_text(path: str, mode: str = "rt"):
    return gzip.open(path, mode, encoding="utf-8", newline="") if path.endswith(".gz") else open(path, mode, encoding="utf-8", newline="")


def normalize_id(value: str, strip_mate_suffix: bool) -> str:
    value = value.strip().lstrip("@").split()[0]
    if "|:|" in value:  # Kraken2 paired-output identifier convention
        value = value.split("|:|", 1)[0]
    if strip_mate_suffix and value.endswith(("/1", "/2")):
        value = value[:-2]
    return value


def load_manifest(path: str, sample_id: str):
    rows, aliases = [], {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "query_id"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Manifest requires columns: {sorted(required)}")
        for row in reader:
            if row["sample_id"] != sample_id:
                continue
            query_id = row["query_id"]
            if query_id in aliases:
                raise ValueError(f"Duplicate manifest query_id: {query_id}")
            rows.append(query_id)
            for key in ("query_id", "read1_id", "read2_id"):
                if row.get(key):
                    aliases[row[key].split()[0]] = query_id
    if not rows:
        raise ValueError(f"No manifest rows found for sample_id={sample_id}")
    return rows, aliases


def parse_kraken2(path: str, strip_suffix: bool):
    with open_text(path) as handle:
        for line_no, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or fields[0] not in {"C", "U"}:
                raise ValueError(f"Invalid Kraken2 line {line_no}")
            qid = normalize_id(fields[1], strip_suffix)
            yield qid, ("classified" if fields[0] == "C" else "unclassified"), int(fields[2]), "", ""


def parse_kaiju(path: str, strip_suffix: bool):
    with open_text(path) as handle:
        for line_no, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or fields[0] not in {"C", "U"}:
                raise ValueError(f"Invalid Kaiju line {line_no}")
            qid = normalize_id(fields[1], strip_suffix)
            taxid = int(fields[2]) if fields[0] == "C" and len(fields) > 2 else 0
            score = fields[3] if fields[0] == "C" and len(fields) > 3 else ""
            yield qid, ("classified" if fields[0] == "C" else "unclassified"), taxid, score, ("kaiju_score" if score else "")


def parse_centrifuge(path: str, strip_suffix: bool):
    best = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        needed = {"readID", "taxID", "score"}
        if not needed.issubset(reader.fieldnames or []):
            raise ValueError(f"Centrifuge output requires columns: {sorted(needed)}")
        for row in reader:
            qid = normalize_id(row["readID"], strip_suffix)
            score = float(row["score"])
            if qid not in best or score > best[qid][0]:
                best[qid] = (score, int(row["taxID"]))
    for qid, (score, taxid) in best.items():
        yield qid, "classified", taxid, str(score), "centrifuge_score"


def parse_generic(path: str, strip_suffix: bool, id_column: str, taxid_column: str):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {id_column, taxid_column}.issubset(reader.fieldnames or []):
            raise ValueError(f"Generic input requires columns {id_column!r} and {taxid_column!r}")
        for row in reader:
            qid = normalize_id(row[id_column], strip_suffix)
            taxid = int(row[taxid_column] or 0)
            yield qid, ("classified" if taxid > 0 else "unclassified"), taxid, "", ""


def write_output(args, manifest_ids, aliases, records):
    assignments = {}
    for raw_id, status, taxid, score, score_type in records:
        query_id = aliases.get(raw_id, raw_id)
        if query_id not in aliases:
            raise ValueError(f"Classifier ID is absent from manifest: {raw_id}")
        if query_id in assignments:
            raise ValueError(f"Multiple assignments map to query_id {query_id}; reconcile mates/hits explicitly")
        assignments[query_id] = (status, taxid, score, score_type)

    out_fields = ["analysis_id", "sample_id", "query_id", "status", "tax_id", "score", "score_type", "comment"]
    with open_text(args.output, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=out_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for query_id in manifest_ids:
            status, taxid, score, score_type = assignments.get(query_id, (args.missing_status, 0, "", ""))
            writer.writerow({"analysis_id": args.analysis_id, "sample_id": args.sample_id, "query_id": query_id,
                             "status": status, "tax_id": taxid, "score": score,
                             "score_type": score_type, "comment": ""})


def validate(path: str, manifest_ids, analysis_id: str, sample_id: str):
    seen, counts = set(), {s: 0 for s in STATUSES}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"analysis_id", "sample_id", "query_id", "status", "tax_id"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"read2tax file requires columns: {sorted(required)}")
        for line_no, row in enumerate(reader, 2):
            if row["analysis_id"] != analysis_id or row["sample_id"] != sample_id:
                raise ValueError(f"Unexpected analysis/sample ID at line {line_no}")
            qid, status, taxid = row["query_id"], row["status"], int(row["tax_id"])
            if qid in seen:
                raise ValueError(f"Duplicate query_id at line {line_no}: {qid}")
            if status not in STATUSES or (status == "classified") != (taxid > 0):
                raise ValueError(f"Invalid status/tax_id combination at line {line_no}")
            seen.add(qid); counts[status] += 1
    expected = set(manifest_ids)
    missing, extra = expected - seen, seen - expected
    if missing or extra:
        raise ValueError(f"Manifest mismatch: {len(missing)} missing, {len(extra)} extra query IDs")
    print(f"VALID\trows={len(seen)}\t" + "\t".join(f"{k}={v}" for k, v in sorted(counts.items())))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--format", choices=["kraken2", "centrifuge", "kaiju", "generic", "validate"], required=True)
    parser.add_argument("--input", required=True, help="Classifier output or standardized read2tax TSV[.gz]")
    parser.add_argument("--manifest", required=True, help="Organiser query manifest TSV[.gz]")
    parser.add_argument("--analysis-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output", help="Output read2tax TSV[.gz]; not used with --format validate")
    parser.add_argument("--strip-mate-suffix", action="store_true", help="Strip terminal /1 or /2 before manifest lookup")
    parser.add_argument("--id-column", default="query_id", help="Generic input ID column")
    parser.add_argument("--taxid-column", default="tax_id", help="Generic input TaxID column")
    parser.add_argument("--missing-status", choices=["unclassified", "excluded_preclassification", "error"], default="excluded_preclassification")
    args = parser.parse_args()
    manifest_ids, aliases = load_manifest(args.manifest, args.sample_id)
    if args.format == "validate":
        validate(args.input, manifest_ids, args.analysis_id, args.sample_id); return
    if not args.output:
        parser.error("--output is required for conversion")
    parsers = {
        "kraken2": lambda: parse_kraken2(args.input, args.strip_mate_suffix),
        "centrifuge": lambda: parse_centrifuge(args.input, args.strip_mate_suffix),
        "kaiju": lambda: parse_kaiju(args.input, args.strip_mate_suffix),
        "generic": lambda: parse_generic(args.input, args.strip_mate_suffix, args.id_column, args.taxid_column),
    }
    write_output(args, manifest_ids, aliases, parsers[args.format]())
    validate(args.output, manifest_ids, args.analysis_id, args.sample_id)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, csv.Error) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
