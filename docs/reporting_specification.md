# Reporting specification

## 1. Scope

This document defines the metadata and taxonomic read-classification output to
be submitted for the 2026 computational metagenomics ring trial. Its purpose is
to make workflows and results comparable without prescribing a particular
classifier, database, or pipeline.

The spreadsheet template is used to describe the workflow. Taxonomic results
are submitted separately as tab-separated `read2tax` files because these files
may contain millions of rows.

The keywords **must**, **must not**, **required**, **should**, and **optional**
indicate the strength of each requirement.

## 2. Submission package

Each laboratory must submit:

1. one completed metadata workbook;
2. one `read2tax.tsv.gz` file for each combination of `analysis_id` and
   `sample_id`;
3. a SHA-256 checksum for every `read2tax` file; and
4. optionally, supplementary reports or native classifier output.

The completed metadata workbook is the authoritative description of the
analysis. Native classifier reports do not replace the standardized `read2tax`
files.

## 3. Stable identifiers

### 3.1 Laboratory code

The organisers assign each participant a stable `lab_code`, for example
`LAB01`. The code must be used consistently in filenames and analysis IDs.

### 3.2 Sample identifier

The `sample_id` identifies one distributed dataset, for example `DATASET_01`.
Participants must not rename or reinterpret these identifiers.

### 3.3 Analysis identifier

An `analysis_id` identifies one complete analytical workflow, for example
`LAB01_A01`.

- The same `analysis_id` must be used for all datasets processed with the same
  workflow.
- A materially different workflow, database strategy, or parameterization must
  receive a different `analysis_id`.
- Each laboratory must designate exactly one analysis as its primary analysis.
- The primary analysis should represent the participant's usual workflow.

### 3.4 Query identifier

The `query_id` is provided by the organisers in a dataset-specific query
manifest. It is the authoritative identifier for one benchmark unit.

The unique key of a result row is:

```text
analysis_id + sample_id + query_id
```

## 4. Query manifests and reporting units

The organisers provide one query manifest per dataset. Participants must use
the manifest identifiers exactly as supplied.

### 4.1 Single-end and Nanopore datasets

One query represents one sequencing read. The manifest `query_id` corresponds
to the read identifier supplied by the organisers.

### 4.2 Paired-end Illumina datasets

One query represents one DNA fragment/read pair. The two mates therefore share
one benchmark `query_id`.

This definition avoids counting paired fragments twice and supports tools that
classify paired reads jointly. If a workflow classifies the mates separately,
the participant must use a deterministic reconciliation rule to produce one
final assignment for the fragment. The rule must be documented in the
`Workflow_Steps` sheet.

Examples of possible reconciliation rules include a taxonomy-aware lowest
common ancestor or a documented score-based decision. The reporting format does
not prescribe one rule because it should reflect the participant's actual
workflow. Silently selecting the first or last mate is not acceptable.

If one mate is removed but the other is classified, the fragment's status is
determined by the final classification result. A fragment is
`excluded_preclassification` only when it did not reach the taxonomic
classification stage at all.

## 5. Metadata workbook

Participants must not rename worksheets, columns, field identifiers, controlled
vocabulary values, dataset IDs, or analysis IDs after these have been assigned.
Additional explanation should be placed in the available notes or comments
columns.

Orange cells in the template are intended for participant input.

### 5.1 `Submission`

Contains laboratory identity, contact information, submission date, template
version, and general declarations. This sheet is completed once per laboratory.

### 5.2 `Analyses`

Contains one row for each complete workflow. Participants should record:

- whether the analysis is primary;
- pipeline or workflow name and version;
- public repository or persistent identifier, if available;
- execution environment and operating system;
- workflow engine and container technology;
- execution date; and
- relevant general comments.

A public pipeline URL does not replace exact version information. Prefer a
release, tag, commit identifier, container digest, or other immutable reference.

### 5.3 `Analysis_Samples`

Maps analyses to datasets. Any dataset-specific deviation from the general
workflow must be documented here.

### 5.4 `Workflow_Steps`

Contains one row per ordered analytical step and analysis. The `step_order`
must describe the actual order in which data were processed.

For every applicable step, participants should provide:

- a controlled `step_category`;
- tool name and exact version;
- command line or relevant non-default parameters;
- input and output of the step;
- whether the step was applied to all datasets; and
- explanatory notes where necessary.

Unused example rows may be removed. Missing steps should not be invented merely
to fill the template. If several tools are used within one conceptual step,
record them as separate ordered rows.

The controlled categories are:

| Category | Meaning |
| --- | --- |
| `raw_read_qc` | Inspection or quality assessment without removing reads |
| `adapter_or_primer_removal` | Adapter, barcode, or primer trimming |
| `quality_filtering` | Length, quality, or complexity filtering |
| `host_depletion` | Removal of host-associated reads |
| `read_correction` | Read correction or other preprocessing |
| `taxonomic_classification` | Primary taxonomic assignment |
| `abundance_estimation` | Abundance re-estimation following classification |
| `postclassification_filtering` | Thresholds or filters applied to assignments |
| `taxonomic_reassignment` | LCA, reassignment, or rank normalization |
| `visualization_or_reporting` | Report generation or visualization |
| `other` | Another step explained in the notes column |

### 5.5 `Databases`

Contains one row per exact database snapshot or build. Record, where applicable:

- database name and provider;
- database release;
- snapshot or download date;
- taxonomy system and release;
- source URL or accession;
- build command and relevant settings; and
- custom additions, exclusions, or modifications.

The database name alone is insufficient. For rolling databases, the download
date and build procedure are particularly important.

### 5.6 `Step_Databases`

Links workflow steps to database records. Use more than one row if a single
step uses several databases. Each `database_id` must refer to a row in the
`Databases` sheet.

### 5.7 `Sample_Results`

Contains one row per `analysis_id` and `sample_id`. It records the submitted
filename, checksum, query counts, and number of distinct positive TaxIDs.

The following count must hold:

```text
input_query_units
  = classified
  + unclassified
  + excluded_preclassification
  + error
```

The workbook calculates the reported total and difference from the input count.
A non-zero difference must be resolved before submission.

## 6. `read2tax` file format

### 6.1 Technical format

Each result file must:

- use UTF-8 encoding;
- use tabs as delimiters;
- contain exactly one header row;
- contain no comment or metadata lines before or after the table;
- contain exactly one row for every query in the corresponding manifest;
- contain no duplicate unique keys;
- preferably be compressed with gzip; and
- use the filename
  `<analysis_id>.<sample_id>.read2tax.tsv.gz`.

Rows may appear in any order, although manifest order is recommended.

### 6.2 Columns

| Column | Requirement | Type | Definition |
| --- | --- | --- | --- |
| `analysis_id` | Required | Text | Stable ID matching the metadata workbook |
| `sample_id` | Required | Text | Exact organiser-provided dataset ID |
| `query_id` | Required | Text | Exact query ID from the organiser manifest |
| `status` | Required | Controlled text | Classification outcome defined below |
| `tax_id` | Required | Integer | Positive NCBI Taxonomy ID when classified; otherwise `0` |
| `score` | Optional | Number | Native, untransformed classifier score or confidence |
| `score_type` | Conditional | Text | Required whenever `score` is present |
| `comment` | Optional | Text | Query-specific explanation |

The first five columns are mandatory. Optional columns should retain the names
shown above.

### 6.3 Classification status

| Status | Required `tax_id` | Meaning |
| --- | ---: | --- |
| `classified` | Positive integer | The primary workflow produced a final assignment |
| `unclassified` | `0` | The query reached classification but no assignment was produced |
| `excluded_preclassification` | `0` | The query was removed before classification, for example by QC or host depletion |
| `error` | `0` | The query could not be processed because of a technical error |

An absent row is not equivalent to `unclassified`. Every manifest query must be
represented explicitly so that filtering and missing output can be distinguished
from failure to classify.

### 6.4 Taxonomic assignment

Report the final assignment from the submitted workflow. Do not force every
assignment to species rank. Assignments to genus, family, or higher ranks are
valid when that is the workflow's final result.

`tax_id` refers to an NCBI Taxonomy identifier. The corresponding taxonomy
release or snapshot must be documented in the metadata workbook. Taxon names
and ranks should not be included as authoritative fields because organisers can
derive them consistently from the documented taxonomy.

Participants whose primary workflow uses a taxonomy that cannot be mapped
reliably to NCBI Taxonomy, such as a purely GTDB-based identifier scheme, should
contact the organisers before analysis. Native assignments may additionally be
submitted as supplementary output, but an undocumented or lossy mapping should
not be performed silently.

### 6.5 Scores

Scores are optional because tools calculate fundamentally different quantities.
If a score is reported:

- it must be the native, untransformed value;
- `score_type` must state what the value represents; and
- any threshold applied to the score must be documented under
  `Workflow_Steps`.

Scores from different tools must not be assumed to be directly comparable.

### 6.6 Example

```text
analysis_id\tsample_id\tquery_id\tstatus\ttax_id\tscore\tscore_type\tcomment
LAB01_A01\tDATASET_01\tA00123:42:H3V5YDSX7:1:1101:1000:1000\tclassified\t562\t0.984\tclassifier_confidence\t
LAB01_A01\tDATASET_01\tA00123:42:H3V5YDSX7:1:1101:1000:1001\tunclassified\t0\t\t\t
LAB01_A01\tDATASET_01\tA00123:42:H3V5YDSX7:1:1101:1000:1002\texcluded_preclassification\t0\t\t\thost depletion
```

The repository also contains a downloadable example under `examples/`.

## 7. Conversion and validation

Participants may create the standardized table directly or use
`scripts/convert_read2tax.py`. The helper supports:

- Kraken2 output;
- Centrifuge classification output;
- Kaiju output;
- generic two-column query-ID/TaxID tables; and
- validation of an already standardized file.

The helper fills classifier IDs missing from the raw output using the status
specified by `--missing-status`. Its default is
`excluded_preclassification`. Participants must verify that this accurately
describes their workflow. If the raw classifier omits unclassified reads, use
`--missing-status unclassified` instead.

Successful validation checks:

- required columns;
- expected analysis and sample IDs;
- allowed status values;
- consistency between `status` and `tax_id`;
- duplicate query IDs; and
- completeness relative to the query manifest.

Successful structural validation does not establish that the biological
assignment or metadata are correct.

## 8. Multiple or alternative analyses

Alternative approaches are welcome when they are clearly separated from the
primary workflow.

- Give each alternative a distinct `analysis_id`.
- Add corresponding rows to all applicable workbook sheets.
- Submit a separate `read2tax` file per dataset and alternative.
- Do not combine results from several classifiers in one analysis unless that
  combination is an explicit part of the workflow.
- Describe ensemble, consensus, or arbitration rules as dedicated workflow
  steps.

## 9. Checksums and final quality control

Generate SHA-256 checksums after compression:

```bash
sha256sum LAB01_A01.DATASET_01.read2tax.tsv.gz \
  > LAB01_A01.DATASET_01.read2tax.tsv.gz.sha256
```

On macOS, use:

```bash
shasum -a 256 LAB01_A01.DATASET_01.read2tax.tsv.gz \
  > LAB01_A01.DATASET_01.read2tax.tsv.gz.sha256
```

Before submission, confirm that:

- exactly one analysis is marked as primary;
- every included analysis is linked to the applicable datasets;
- tool and database versions are specific and reproducible;
- all result files pass validation;
- workbook counts agree with the validated files;
- filenames contain the correct identifiers; and
- checksums were calculated from the final submitted files.

## 10. Notes for organisers

Before releasing the repository, the organisers should:

1. replace provisional dataset descriptions with the final metadata;
2. generate the four authoritative query manifests;
3. verify that paired-end read identifiers map unambiguously to fragments;
4. publish checksums for datasets, manifests, and templates;
5. test conversion and validation with at least one Illumina and one Nanopore
   workflow;
6. freeze and tag the repository version used for the trial;
7. decide how host-associated assignments will be handled in evaluation; and
8. define whether the result workshop uses coded or named laboratory results.

If a biological truth set is unavailable for a dataset, results should be
described as interlaboratory agreement rather than accuracy. Evaluation should
consider several taxonomic ranks and separately report classified,
unclassified, and pre-classification-excluded fractions.

