# Metagenomics Computational Ring Trial 2026

This repository contains the reporting templates, format specification, example
files, query manifests, and validation utilities for the 2026 computational
metagenomics ring trial organised by the Friedrich-Loeffler-Institut (FLI) and
the Robert Koch Institute (RKI).

The aim is to compare the taxonomic classification results produced by the
participants' usual metagenomics workflows. There are no prescribed analysis
tools. Participants should document the software, versions, parameters,
databases, and relevant filtering steps sufficiently to make the results
interpretable and reproducible.

## Important dates

| Milestone | Date |
| --- | --- |
| Start of the computational ring trial | 1 March 2026 |
| Submission of results | 31 May 2026 |
| Joint results workshop at FLI | September 2026 |

## Datasets

Four metagenomic datasets are provided separately through the download route
communicated by the organisers:

| Dataset ID | Sequencing | Query unit used for reporting |
| --- | --- | --- |
| `DATASET_01` | Illumina | One fragment for paired-end data; one read for single-end data |
| `DATASET_02` | Illumina | One fragment for paired-end data; one read for single-end data |
| `DATASET_03` | Illumina | One fragment for paired-end data; one read for single-end data |
| `DATASET_04` | Oxford Nanopore | One read |

Sequencing files are not stored in this repository. Please verify all downloaded
files against the checksums supplied by the organisers.

## Quick start

1. Download the four sequencing datasets and verify their checksums.
2. Analyse every dataset using your laboratory's usual taxonomic metagenomics
   workflow.
3. Complete
   [`templates/Metagenomics_Ring_Trial_Reporting_Template.xlsx`](templates/Metagenomics_Ring_Trial_Reporting_Template.xlsx).
4. Generate one standardized `read2tax` result file for each dataset and
   analysis. You may use
   [`scripts/convert_read2tax.py`](scripts/convert_read2tax.py) or generate the
   standardized format directly.
5. Validate every `read2tax` file against the corresponding organiser-provided
   query manifest.
6. Submit the completed workbook, `read2tax` files, and SHA-256 checksums through
   the submission route specified by the organisers.

Detailed requirements are defined in
[`docs/reporting_specification.md`](docs/reporting_specification.md).

## Repository contents

```text
.
├── README.md
├── CHANGELOG.md
├── templates/
│   └── Metagenomics_Ring_Trial_Reporting_Template.xlsx
├── examples/
│   ├── read2tax_example.tsv
│   └── query_manifest_example.tsv
├── scripts/
│   └── convert_read2tax.py
├── manifests/
│   ├── DATASET_01.query_manifest.tsv.gz
│   ├── DATASET_02.query_manifest.tsv.gz
│   ├── DATASET_03.query_manifest.tsv.gz
│   ├── DATASET_04.query_manifest.tsv.gz
│   └── checksums.sha256
└── docs/
    └── reporting_specification.md
```

## Required submission

Each laboratory should submit:

- one completed metadata workbook;
- one gzip-compressed `read2tax` TSV for every dataset and submitted analysis;
- a SHA-256 checksum for every submitted result file; and
- any optional supplementary output the participant considers helpful for
  interpretation.

Each participant must designate exactly one analysis as the **primary
analysis**. This should represent the laboratory's usual workflow. Additional
or experimental workflows may be submitted under separate `analysis_id`
values.

Recommended filenames:

```text
<lab_code>.metadata.xlsx
<analysis_id>.<sample_id>.read2tax.tsv.gz
<analysis_id>.<sample_id>.read2tax.tsv.gz.sha256
```

## `read2tax` output

The five required columns are:

```text
analysis_id    sample_id    query_id    status    tax_id
```

Optional columns are `score`, `score_type`, and `comment`. The file must contain
exactly one row for every query in the corresponding query manifest, including
queries that were unclassified or removed before classification.

For paired-end Illumina data, one `query_id` represents an entire fragment/read
pair. For Nanopore data, one `query_id` represents one read.

## Conversion examples

The conversion helper supports Kraken2, Centrifuge, Kaiju, and generic
two-column TSV output. For example:

```bash
python scripts/convert_read2tax.py \
  --format kraken2 \
  --input results/DATASET_01.kraken2.tsv \
  --manifest manifests/DATASET_01.query_manifest.tsv.gz \
  --analysis-id LAB01_A01 \
  --sample-id DATASET_01 \
  --strip-mate-suffix \
  --output LAB01_A01.DATASET_01.read2tax.tsv.gz
```

Validate a standardized file with:

```bash
python scripts/convert_read2tax.py \
  --format validate \
  --input LAB01_A01.DATASET_01.read2tax.tsv.gz \
  --manifest manifests/DATASET_01.query_manifest.tsv.gz \
  --analysis-id LAB01_A01 \
  --sample-id DATASET_01
```

Run `python scripts/convert_read2tax.py --help` for all options.

> [!IMPORTANT]
> The converter deliberately stops when several classifier records map to the
> same benchmark query. This commonly occurs when paired mates were classified
> independently. Participants must apply and document an explicit reconciliation
> rule rather than allowing the script to select one mate silently.

## Versioning

Use the repository release specified by the organisers. The template version
and repository release used for the analysis must be recorded in the reporting
workbook. Corrections or format changes after the initial release will be listed
in `CHANGELOG.md`.

## Questions and submission

General questions that may be relevant to all participants can be raised through
the repository's issue tracker, if enabled. Questions involving confidential
sample information, participant-specific results, or file transfer should be
sent directly to the organisers:

- Dirk Höper, FLI: `Dirk.Hoeper@fli.de`
- Andrea Thürmer, RKI: `thuermera@rki.de`

Files larger than 5 MB should not be sent by email. Contact the organisers for
an upload link.
