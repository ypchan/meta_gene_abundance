# meta-gene-abundance

`meta-gene-abundance` estimates target-gene abundance from metagenomic FASTQ reads.
It maps reads to a target-gene nucleotide FASTA with `minimap2`, filters primary
alignments with `samtools`, resolves paired-end target conflicts, and writes one
gene-level long table containing raw counts, RPM/FPM, RPKM/FPKM, and coverage.

The implementation intentionally remains a single-script style tool:

- repository entry point: `gene_abundance.py`
- installable command: `gene_abundance`
- no workflow framework or package layout required

## Installation

### Install from Git

```bash
pip install git+https://github.com/ypchan/meta_gene_abundance.git
```

After installation:

```bash
gene_abundance --help
```

### Install from a Local Clone

```bash
git clone https://github.com/ypchan/meta_gene_abundance.git
cd meta_gene_abundance
pip install -e .
```

### External Command Dependencies

The Python package installs only Python dependencies. The analysis also requires
these executables in `PATH`:

- `minimap2`
- `samtools`
- `seqkit`

A typical conda setup is:

```bash
conda install -c bioconda minimap2 samtools seqkit
```

## Quick Start

Create a tab-delimited sample list with a header:

```text
sample_id	r1_path	r2_path	year	month	depth
sample01	/path/sample01.R1.fq.gz	/path/sample01.R2.fq.gz	2024	7	0-2
sample02	/path/sample02.R1.fq.gz	/path/sample02.R2.fq.gz	2024	7	2-5
```

Run the analysis:

```bash
gene_abundance \
  -i samples.tsv \
  --reference target_genes.fa \
  --prefix gene_abundance \
  --jobs 4 \
  --threads 8
```

Main outputs:

```text
gene_abundance.gene_abundance.long.tsv
gene_abundance.sample_qc.tsv
gene_abundance.reference.tsv
gene_abundance.run.log
gene_abundance.state.json
```

## Workflow

```text
target_genes.fa
      |
      | minimap2 -x sr -d
      v
reference.mmi
      |
      v
FASTQ -> seqkit stats -> minimap2 -ax sr --secondary=no
      -> samtools view filter -> retained BAM
      -> fragment target assignment
      -> samtools coverage
      -> long TSV + sample QC + checkpoint
```

The tool supports either an input FASTA or a prebuilt minimap2 index:

```text
--reference target_genes.fa  -> build a temporary .mmi index for this run
--index reference.mmi        -> reuse an existing minimap2 .mmi index
```

Recommended prebuild command:

```bash
minimap2 -x sr -d reference.mmi target_genes.fa
```

When `--index` is used, target IDs and lengths are read from the minimap2 SAM
header emitted by the index.

## Input Modes

Choose exactly one read input mode.

### Sample-List Mode

Use `-i`, `-inlist`, or `--inlist`.

The file must be tab-delimited and must have a header.

Required columns:

```text
sample_id
r1_path
```

Optional columns:

```text
r2_path
year
month
depth
```

Rules:

- `sample_id` must be non-empty and unique.
- `r1_path` must point to a non-empty FASTQ or FASTQ.GZ file.
- `r2_path` may be empty, `NA`, `None`, `-`, or `.` for single-end data.
- Relative FASTQ paths are resolved relative to the input-list directory.
- Paired-end R1 and R2 files must contain the same number of records. This is
  checked with `seqkit stats`.
- Extra columns are ignored.

### Single-Sample Mode

```bash
gene_abundance \
  --r1 sample.R1.fq.gz \
  --r2 sample.R2.fq.gz \
  --sample-id sample01 \
  --reference target_genes.fa \
  --prefix sample01.genes
```

For single-end reads, omit `--r2`.

If `--sample-id` is omitted, it is inferred from the R1 filename after removing
common FASTQ and R1 suffixes.

## Reference FASTA Requirements

- Each FASTA record represents one target gene sequence.
- The first whitespace-delimited token in each FASTA header is used as
  `gene_id`.
- `gene_id` values must be unique.
- For metagenomic DNA abundance, use genomic gene sequences when that is the
  intended target. Avoid mixing incompatible target definitions.
- Highly redundant alleles or exact duplicates should be reduced before mapping
  if you want high MAPQ and unambiguous assignment.

## Key Parameters

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `-i`, `-inlist`, `--inlist` | none | Tab-delimited sample list with header. |
| `--r1` | none | R1 FASTQ or single-end FASTQ. |
| `--r2` | none | R2 FASTQ for paired-end single-sample mode. |
| `--reference` | none | Target-gene FASTA. Mutually exclusive with `--index`. |
| `--index` | none | Existing minimap2 `.mmi` index. Mutually exclusive with `--reference`. |
| `--prefix` | required | Output prefix. The parent directory must exist. |
| `--jobs` | `1` | Number of samples processed concurrently. |
| `--threads` | `8` | Approximate CPU budget per sample. |
| `--min-mapq` | `20` | Minimum MAPQ retained by `samtools view`. |
| `--sort-memory` | `1G` | Memory per samtools sort thread. |
| `--min-fragments` | `5` | Detection threshold for assigned fragments. |
| `--min-breadth` | `50.0` | Detection threshold for coverage breadth percent. |
| `--keep-bam` | off | Keep coordinate-sorted BAM/BAI files in `PREFIX.bam/`. |
| `--force` | off | Delete previous outputs and checkpoint, then recompute. |
| `--no-resume` | off | Ignore reusable checkpoint entries and recompute samples. |
| `--tmp-dir` | system temp | Parent directory for temporary files. |

## Mapping and Filtering

For each sample, reads are mapped as:

```bash
minimap2 -ax sr --secondary=no -t THREADS reference.mmi R1 [R2]
```

The SAM stream is filtered into BAM with:

```bash
samtools view -b -q MIN_MAPQ -F 2308
```

`-F 2308` removes:

- unmapped records
- secondary alignments
- supplementary alignments

The output therefore uses retained primary alignments with MAPQ at least
`--min-mapq`.

## Paired-End Conflict Resolution

Fragment assignment is done after name collation with `samtools collate`.

```text
read pair
  |
  +-- both ends retained, same gene
  |      -> assign 1 fragment and 2 reads to that gene
  |
  +-- only one end retained
  |      -> assign 1 fragment and 1 read to that gene
  |
  +-- both ends retained, different genes
         -> compare alignment quality
```

For different-gene paired-end conflicts, the winner is selected by:

```text
AS > MAPQ > aligned_query_length > lower_NM
```

The winning target receives the whole fragment and both read-end counts. If the
quality key is exactly tied, the fragment is marked unresolved and excluded from
assigned abundance.

## Boundary and Partial Alignments

No special boundary-only rule is applied.

Reads mapping near gene starts or ends are handled by normal minimap2 SAM output:

- soft-clipped or partially aligned reads can be retained if they pass MAPQ and
  primary-alignment filters
- clipped bases outside the reference do not contribute to coverage
- `aligned_query_length` is computed from CIGAR `M`, `I`, `=`, and `X`
- coverage breadth and mean depth are computed by `samtools coverage` from the
  retained BAM

This means boundary evidence can contribute to abundance, but reliability is
controlled by MAPQ, conflict resolution, coverage breadth, and the detection
thresholds.

## Output Files

### `PREFIX.gene_abundance.long.tsv`

One row per completed sample and target gene. This is the main abundance table.

Important columns:

| Column | Meaning |
| --- | --- |
| `sample_id`, `year`, `month`, `depth` | Sample metadata. Metadata fields come from the input list when present. |
| `library_type` | `PE` or `SE`. |
| `gene_id` | Target gene ID from FASTA or minimap2 index header. |
| `gene_length` | Target gene length. |
| `total_clean_reads` | Total input reads counted by `seqkit stats`. |
| `total_clean_fragments` | Input fragments. For paired-end data this is the R1 count. |
| `aligned_reads` | Retained aligned reads reported by `samtools coverage`. |
| `assigned_reads` | Reads assigned by the fragment resolver. |
| `assigned_fragments` | Fragments assigned by the fragment resolver. |
| `covered_bases` | Covered reference bases from `samtools coverage`. |
| `breadth_pct` | Percent of gene bases covered. |
| `mean_depth` | Mean coverage depth. |
| `mean_baseq` | Mean base quality from retained alignments. |
| `mean_mapq` | Mean MAPQ from retained alignments. |
| `read_rpm` | `1e6 * assigned_reads / total_clean_reads`. |
| `fragment_fpm` | `1e6 * assigned_fragments / total_clean_fragments`. |
| `rpkm` | `1e9 * assigned_reads / (gene_length * total_clean_reads)`. |
| `fpkm` | `1e9 * assigned_fragments / (gene_length * total_clean_fragments)`. |
| `detected` | `1` if `assigned_fragments >= --min-fragments` and `breadth_pct >= --min-breadth`; otherwise `0`. |

### `PREFIX.sample_qc.tsv`

One row per completed sample. It records library size, retained alignment
records, target-assigned counts, conflict counts, timing, and thread allocation.

### `PREFIX.reference.tsv`

A two-column reference dictionary:

```text
gene_id	gene_length
```

This file is redundant with the `gene_id` and `gene_length` columns in the long
table, but it is useful for auditing the exact target set and for downstream
tools that want a separate gene-length table.

### `PREFIX.run.log`

Run log with command progress and errors.

### `PREFIX.state.json`

Atomic checkpoint used for resume. Completed samples are reused only when their
FASTQ path, size, and modification time match the checkpoint and the mapping
signature is unchanged.

### `PREFIX.failed.tsv`

Written only when one or more samples fail. Successful samples remain in the
checkpoint and can be reused after the failed inputs are corrected.

### `PREFIX.bam/`

Created only with `--keep-bam`. Contains coordinate-sorted, indexed BAM files
for retained alignments.

## Resume Behavior

The checkpoint stores completed sample results in `PREFIX.state.json`.

Reusable entries require:

- same FASTQ paths
- same FASTQ file sizes
- same FASTQ modification times
- same target IDs
- same reference/index and mapping settings

Changing display thresholds such as `--min-fragments` or `--min-breadth` rewrites
output tables without remapping. Changing reference, index, MAPQ, or read files
requires recomputation with `--force`.

## Notes

- The tool does not create an analysis directory; outputs are written as
  `PREFIX.*` files.
- Older `PREFIX.gene_abundance.rpkm.tsv` and
  `PREFIX.gene_abundance.fpkm.tsv` files from previous versions are deleted
  during output writing because the long table now contains both metrics.
