#!/usr/bin/env python3
"""Create ring-trial query manifests and SHA-256 checksums from FASTQ files.

The script supports plain or gzip-compressed FASTQ input. For paired-end data,
one manifest row represents one fragment/read pair. For single-end or Nanopore
data, one row represents one read.

Single dataset example:

    python create_query_manifest.py \
      --sample-id DATASET_01 \
      --read1 data/DATASET_01_R1.fastq.gz \
      --read2 data/DATASET_01_R2.fastq.gz \
      --output manifests/DATASET_01.query_manifest.tsv.gz \
      --checksum-file manifests/DATASET_01.checksums.sha256

Batch example:

    python create_query_manifest.py \
      --config datasets.tsv \
      --checksum-file manifests/checksums.sha256

The batch configuration requires these tab-separated columns:

    sample_id  read1  read2  output_manifest

Leave read2 empty for single-end or Nanopore data. Relative paths in the batch
configuration are resolved relative to the configuration file.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import os
import re
import sqlite3
import sys
import tempfile
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO

VERSION = "1.0.0"
MANIFEST_HEADER = ["sample_id", "query_id", "input_unit", "read1_id", "read2_id"]
SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


class ManifestError(Exception):
    """Raised when input data cannot produce a valid manifest."""


@dataclass(frozen=True)
class Dataset:
    sample_id: str
    read1: Path
    read2: Path | None
    output: Path


def fastq_text(path: Path) -> TextIO:
    """Open a plain or gzip-compressed FASTQ as strict UTF-8 text."""
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="")
    return path.open("rt", encoding="utf-8", errors="strict", newline="")


def clean_line(line: str) -> str:
    return line.rstrip("\r\n")


def fastq_ids(path: Path) -> Iterator[str]:
    """Yield the first whitespace-delimited token of every FASTQ header."""
    with fastq_text(path) as handle:
        record_number = 0
        while True:
            header = handle.readline()
            if header == "":
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            record_number += 1

            if "" in (sequence, plus, quality):
                raise ManifestError(
                    f"{path}: truncated FASTQ record {record_number}"
                )
            header = clean_line(header)
            sequence = clean_line(sequence)
            plus = clean_line(plus)
            quality = clean_line(quality)

            if not header.startswith("@"):
                raise ManifestError(
                    f"{path}: record {record_number} header does not start with '@'"
                )
            if not plus.startswith("+"):
                raise ManifestError(
                    f"{path}: record {record_number} third line does not start with '+'"
                )
            if len(sequence) != len(quality):
                raise ManifestError(
                    f"{path}: sequence and quality lengths differ in record {record_number}"
                )

            identifier = header[1:].split(maxsplit=1)[0]
            if not identifier:
                raise ManifestError(f"{path}: empty identifier in record {record_number}")
            if "\t" in identifier:
                raise ManifestError(
                    f"{path}: tab character in identifier at record {record_number}"
                )
            yield identifier


def fragment_id(identifier: str) -> str:
    """Remove a conventional terminal /1 or /2 mate suffix."""
    return identifier[:-2] if identifier.endswith(("/1", "/2")) else identifier


class DuplicateChecker:
    """Exact duplicate detection using memory, SQLite, or no check."""

    def __init__(self, mode: str, parent: Path):
        self.mode = mode
        self.seen: set[str] | None = set() if mode == "memory" else None
        self.connection: sqlite3.Connection | None = None
        self.database_path: Path | None = None
        if mode == "sqlite":
            descriptor, name = tempfile.mkstemp(
                prefix=".manifest_ids.", suffix=".sqlite", dir=parent
            )
            os.close(descriptor)
            self.database_path = Path(name)
            self.connection = sqlite3.connect(name)
            self.connection.execute("PRAGMA journal_mode=OFF")
            self.connection.execute("PRAGMA synchronous=OFF")
            self.connection.execute("PRAGMA temp_store=FILE")
            self.connection.execute(
                "CREATE TABLE seen (query_id TEXT PRIMARY KEY) WITHOUT ROWID"
            )
            self.connection.execute("BEGIN")

    def add(self, query_id: str) -> None:
        if self.mode == "none":
            return
        if self.seen is not None:
            if query_id in self.seen:
                raise ManifestError(f"duplicate query_id after normalization: {query_id}")
            self.seen.add(query_id)
            return
        assert self.connection is not None
        try:
            self.connection.execute("INSERT INTO seen VALUES (?)", (query_id,))
        except sqlite3.IntegrityError as exc:
            raise ManifestError(
                f"duplicate query_id after normalization: {query_id}"
            ) from exc

    def close(self) -> None:
        if self.connection is not None:
            self.connection.close()
        if self.database_path is not None:
            self.database_path.unlink(missing_ok=True)


@contextmanager
def deterministic_gzip_text(path: Path):
    """Write deterministic gzip output through a temporary file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        with temporary.open("wb") as raw:
            with gzip.GzipFile(
                filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0
            ) as compressed:
                with io.TextIOWrapper(
                    compressed, encoding="utf-8", newline="", write_through=True
                ) as text:
                    yield text
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def create_manifest(dataset: Dataset, duplicate_mode: str, progress_every: int) -> int:
    paired = dataset.read2 is not None
    dataset.output.parent.mkdir(parents=True, exist_ok=True)
    checker = DuplicateChecker(duplicate_mode, dataset.output.parent)
    count = 0
    try:
        with deterministic_gzip_text(dataset.output) as output_handle:
            writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
            writer.writerow(MANIFEST_HEADER)
            read1_iterator = fastq_ids(dataset.read1)

            if paired:
                assert dataset.read2 is not None
                read2_iterator = fastq_ids(dataset.read2)
                while True:
                    read1_id = next(read1_iterator, None)
                    read2_id = next(read2_iterator, None)
                    if read1_id is None and read2_id is None:
                        break
                    if read1_id is None or read2_id is None:
                        raise ManifestError(
                            f"{dataset.sample_id}: paired FASTQ files contain different numbers of records"
                        )
                    query_id = fragment_id(read1_id)
                    read2_query_id = fragment_id(read2_id)
                    if query_id != read2_query_id:
                        raise ManifestError(
                            f"{dataset.sample_id}: unsynchronized pair at record {count + 1}: "
                            f"{read1_id!r} versus {read2_id!r}"
                        )
                    checker.add(query_id)
                    writer.writerow(
                        [dataset.sample_id, query_id, "fragment", read1_id, read2_id]
                    )
                    count += 1
                    if progress_every and count % progress_every == 0:
                        print(
                            f"{dataset.sample_id}: {count:,} fragments processed",
                            file=sys.stderr,
                        )
            else:
                for read1_id in read1_iterator:
                    query_id = read1_id
                    checker.add(query_id)
                    writer.writerow(
                        [dataset.sample_id, query_id, "read", read1_id, ""]
                    )
                    count += 1
                    if progress_every and count % progress_every == 0:
                        print(
                            f"{dataset.sample_id}: {count:,} reads processed",
                            file=sys.stderr,
                        )
    except Exception:
        dataset.output.unlink(missing_ok=True)
        raise
    finally:
        checker.close()

    unit = "fragments" if paired else "reads"
    print(f"CREATED\t{dataset.sample_id}\t{count}\t{unit}\t{dataset.output}")
    return count


def sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def checksum_name(path: Path, checksum_file: Path) -> str:
    try:
        name = os.path.relpath(path, checksum_file.parent)
    except ValueError:
        name = str(path.resolve())
    name = Path(name).as_posix()
    if "\n" in name or "\r" in name:
        raise ManifestError(f"newline in path cannot be represented safely: {path}")
    return name


def write_checksums(paths: list[Path], checksum_file: Path) -> None:
    checksum_file.parent.mkdir(parents=True, exist_ok=True)
    unique_paths: list[Path] = []
    seen: set[Path] = set()
    for path in paths:
        resolved = path.resolve()
        if resolved not in seen:
            seen.add(resolved)
            unique_paths.append(path)

    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{checksum_file.name}.", suffix=".tmp", dir=checksum_file.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        with temporary.open("wt", encoding="utf-8", newline="") as handle:
            for path in unique_paths:
                handle.write(f"{sha256(path)}  {checksum_name(path, checksum_file)}\n")
        os.replace(temporary, checksum_file)
    finally:
        temporary.unlink(missing_ok=True)
    print(f"CREATED\tchecksums\t{len(unique_paths)}\tfiles\t{checksum_file}")


def resolve_config_path(value: str, base: Path) -> Path:
    candidate = Path(value).expanduser()
    return candidate if candidate.is_absolute() else base / candidate


def read_config(path: Path) -> list[Dataset]:
    base = path.resolve().parent
    datasets: list[Dataset] = []
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "read1", "read2", "output_manifest"}
        if not required.issubset(reader.fieldnames or []):
            raise ManifestError(
                "configuration requires columns: " + ", ".join(sorted(required))
            )
        for line_number, row in enumerate(reader, 2):
            if not row["sample_id"] or not row["read1"] or not row["output_manifest"]:
                raise ManifestError(
                    f"{path}: missing required value at line {line_number}"
                )
            datasets.append(
                Dataset(
                    sample_id=row["sample_id"],
                    read1=resolve_config_path(row["read1"], base),
                    read2=(
                        resolve_config_path(row["read2"], base)
                        if row["read2"]
                        else None
                    ),
                    output=resolve_config_path(row["output_manifest"], base),
                )
            )
    if not datasets:
        raise ManifestError(f"{path}: configuration contains no datasets")
    return datasets


def validate_datasets(
    datasets: list[Dataset], checksum_file: Path, force: bool
) -> None:
    sample_ids: set[str] = set()
    outputs: set[Path] = set()
    inputs: set[Path] = set()

    for dataset in datasets:
        if not SAMPLE_ID_RE.fullmatch(dataset.sample_id):
            raise ManifestError(
                f"invalid sample_id {dataset.sample_id!r}; use letters, numbers, '.', '_' or '-'"
            )
        if dataset.sample_id in sample_ids:
            raise ManifestError(f"duplicate sample_id: {dataset.sample_id}")
        sample_ids.add(dataset.sample_id)

        for input_path in [dataset.read1, dataset.read2]:
            if input_path is None:
                continue
            if not input_path.is_file():
                raise ManifestError(f"input FASTQ does not exist: {input_path}")
            inputs.add(input_path.resolve())

        if not dataset.output.name.endswith(".tsv.gz"):
            raise ManifestError(
                f"manifest output must end with .tsv.gz: {dataset.output}"
            )
        resolved_output = dataset.output.resolve()
        if resolved_output in outputs:
            raise ManifestError(f"duplicate manifest output: {dataset.output}")
        if resolved_output in inputs:
            raise ManifestError(f"manifest output overlaps an input FASTQ: {dataset.output}")
        outputs.add(resolved_output)
        if dataset.output.exists() and not force:
            raise ManifestError(
                f"output already exists: {dataset.output}; use --force to replace"
            )

    if checksum_file.resolve() in inputs or checksum_file.resolve() in outputs:
        raise ManifestError("checksum file overlaps an input or manifest output")
    if checksum_file.exists() and not force:
        raise ManifestError(
            f"checksum file already exists: {checksum_file}; use --force to replace"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--config", type=Path, help="Batch configuration TSV")
    mode.add_argument("--sample-id", help="Dataset ID for single-dataset mode")
    parser.add_argument("--read1", type=Path, help="R1, single-end, or Nanopore FASTQ")
    parser.add_argument("--read2", type=Path, help="R2 FASTQ for paired-end data")
    parser.add_argument("--output", type=Path, help="Output query manifest .tsv.gz")
    parser.add_argument(
        "--checksum-file",
        type=Path,
        help="Output SHA-256 file; default: checksums.sha256 beside the output/config",
    )
    parser.add_argument(
        "--duplicate-check",
        choices=["sqlite", "memory", "none"],
        default="sqlite",
        help="Exact query-ID duplicate check (default: sqlite; memory is faster but uses RAM)",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1_000_000,
        help="Print progress every N records; 0 disables progress messages",
    )
    parser.add_argument(
        "--force", action="store_true", help="Replace existing manifests and checksum file"
    )
    args = parser.parse_args()

    if args.progress_every < 0:
        parser.error("--progress-every must be non-negative")
    if args.config is None and (args.read1 is None or args.output is None):
        parser.error("single-dataset mode requires --read1 and --output")
    if args.config is not None and any(
        value is not None for value in (args.read1, args.read2, args.output)
    ):
        parser.error("--read1, --read2, and --output cannot be combined with --config")
    return args


def main() -> int:
    args = parse_args()
    if args.config is not None:
        datasets = read_config(args.config)
        checksum_file = args.checksum_file or args.config.resolve().parent / "checksums.sha256"
    else:
        datasets = [
            Dataset(
                sample_id=args.sample_id,
                read1=args.read1,
                read2=args.read2,
                output=args.output,
            )
        ]
        checksum_file = args.checksum_file or args.output.parent / "checksums.sha256"

    validate_datasets(datasets, checksum_file, args.force)
    checksum_paths: list[Path] = []
    for dataset in datasets:
        create_manifest(dataset, args.duplicate_check, args.progress_every)
        checksum_paths.append(dataset.read1)
        if dataset.read2 is not None:
            checksum_paths.append(dataset.read2)
        checksum_paths.append(dataset.output)
    write_checksums(checksum_paths, checksum_file)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ManifestError, OSError, UnicodeError, csv.Error, sqlite3.Error) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
