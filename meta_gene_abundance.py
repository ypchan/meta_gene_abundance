#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import logging
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, TextIO

from rich.console import Console
from rich.logging import RichHandler
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)
from rich.table import Table
from rich_argparse import RawDescriptionRichHelpFormatter


PROGRAM = "meta_gene_abundance"
CHECKPOINT_SCHEMA = 6
CONFLICT_RULE = "AS>MAPQ>aligned_query_length>lower_NM; exact ties unresolved"
MISSING_R2 = {"", "NA", "Na", "na", "N/A", "None", "none", "-", "."}
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


INPUT_HELP = r"""
Input modes
===========
Choose exactly one analysis input mode:

A. Sample-list mode
-------------------
Use -i/-inlist with a TAB-separated file. The first line must be a header.

Required columns:

    sample_id    r1_path

Optional columns:

    r2_path    year    month    depth

Example (paired-end):

    sample_id\tyear\tmonth\tdepth\tr1_path\tr2_path
    201704_MF1_0002\t2017\t4\t0-2\t/path/sample.R1.fq.gz\t/path/sample.R2.fq.gz

Example (single-end):

    sample_id\tr1_path
    SRR12345678\t/path/SRR12345678.fastq.gz

Input-list rules:
  1. sample_id must be non-empty and unique.
  2. r1_path must point to a non-empty FASTQ/FASTQ.GZ file.
  3. r2_path may be empty, NA, None, -, or . for single-end data.
  4. Relative FASTQ paths are resolved relative to the input-list directory.
  5. Paired files must contain the same number of records; this is checked with
     seqkit stats before mapping.
  6. Extra columns are ignored.

B. Single-sample mode
---------------------
Use --r1 and optionally --r2:

    meta_gene_abundance.py \
        --r1 sample.R1.fq.gz \
        --r2 sample.R2.fq.gz \
        --sample-id sample01 \
        --reference target.fa \
        --prefix sample01.genes

If --sample-id is omitted, it is inferred from the R1 filename after removing
common FASTQ and R1 suffixes.

Reference FASTA requirements
============================
  1. Each record represents one target gene sequence.
  2. The first whitespace-delimited token in each FASTA header is gene_id.
  3. gene_id values must be unique.
  4. For metagenomic DNA abundance, use genomic gene sequences containing
     introns where applicable rather than spliced transcript CDS sequences.
  5. Exact duplicates and highly similar alleles should be reduced before
     mapping to avoid low MAPQ caused by multimapping.

Index input
===========
Minimap2 uses one index file:

    reference.mmi

Choose exactly one reference source:

  1. --reference TARGET.fa
       Build a temporary minimap2 index automatically for this run.

  2. --index reference.mmi
       Reuse an existing minimap2 index.

When --index is used, target gene IDs and lengths are read from the minimap2
SAM header produced by the index. Recommended index command:

    minimap2 -x sr -d reference.mmi reference.fa
"""


DEFAULT_HELP = r"""
Default parameters and analysis definitions
===========================================
Mapper:
  mapper                    minimap2
  preset                    -ax sr
  secondary alignments      --secondary=no
  minimum MAPQ              20
  unmapped/supplementary    removed before counting

Index:
  --reference provided      build a temporary index per run
  --index provided          reuse an existing index
  index command             minimap2 -x sr -d
  index threads             --threads

Paired-end conflict resolution:
  same target on both ends  assign 1 fragment and 2 reads to that target
  only one retained end     assign 1 fragment and 1 read to that target
  different targets         compare AS, then MAPQ, then aligned query length,
                            then lower NM; winning target receives the entire
                            fragment and both read-end counts
  exact quality tie         unresolved and excluded from assigned abundance

Normalization:
  RPKM = 1e9 * assigned_reads / (gene_length * total_clean_reads)
  FPKM = 1e9 * assigned_fragments / (gene_length * total_clean_fragments)

Coverage:
  breadth_pct, mean_depth, mean_baseq and mean_mapq are calculated from actual
  retained alignments. A losing mate reassigned by the conflict rule does not
  artificially increase physical coverage.

Detection flag:
  assigned_fragments >= 5 AND breadth_pct >= 50

Parallelism:
  --jobs 1
  --threads 8 per sample as an approximate CPU budget.

Output organization:
  No analysis directory is created. Results are PREFIX.* files. Gene-level
  abundance values are written to PREFIX.meta_gene_abundance.long.tsv. A single
  PREFIX.bam directory is created only when --keep-bam is specified.

Resume behavior:
  Completed samples are stored atomically in PREFIX.state.json. A sample is
  reused only when its FASTQ path, size and modification time match. Changes to
  display/detection thresholds rewrite output tables without repeating mapping.
"""


@dataclass(frozen=True)
class Sample:
    sample_id: str
    year: str
    month: str
    depth: str
    r1_path: str
    r2_path: str

    @property
    def paired(self) -> bool:
        return self.r2_path not in MISSING_R2


@dataclass(frozen=True)
class ThreadPlan:
    minimap2_threads: int
    samtools_view_threads: int
    samtools_sort_threads: int
    samtools_collate_threads: int
    seqkit_threads: int


@dataclass(frozen=True)
class WorkerConfig:
    index_path: str
    temp_root: str
    target_ids: tuple[str, ...]
    target_lengths: dict[str, int]
    target_set: frozenset[str]
    min_mapq: int
    sort_memory: str
    keep_bam: bool
    bam_dir: str | None
    thread_plan: ThreadPlan


@dataclass(frozen=True)
class Alignment:
    reference: str
    as_score: int
    mapq: int
    aligned_query_length: int
    nm: int

    @property
    def quality_key(self) -> tuple[int, int, int, int]:
        return (self.as_score, self.mapq, self.aligned_query_length, -self.nm)


_PRINT_LOCK = threading.Lock()


def build_console() -> Console:
    width = shutil.get_terminal_size(fallback=(120, 40)).columns
    return Console(width=max(60, width // 2), highlight=False)


CONSOLE = build_console()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=PROGRAM,
        description=(
            "Map metagenomic reads to target genes with minimap2 and "
            "calculate read-level RPKM, fragment-level FPKM and coverage."
        ),
        allow_abbrev=False,
        formatter_class=RawDescriptionRichHelpFormatter,
        epilog=r"""
Examples
========
Run a sample list with a temporary minimap2 index:

meta_gene_abundance.py \
    -i samples.tsv \
    --reference target_genes.fa \
    --prefix results \
    --threads 16

Main outputs
============
PREFIX.meta_gene_abundance.long.tsv
PREFIX.sample_qc.tsv
PREFIX.reference.tsv
PREFIX.run.log
PREFIX.state.json
""",
    )

    required = parser.add_argument_group("required arguments")
    input_mode = required.add_mutually_exclusive_group(required=True)
    input_mode.add_argument(
        "-i", "-inlist", "--inlist", dest="manifest", type=Path, metavar="INLIST",
        help=(
            "TAB-separated sample list for multi-sample/list mode. Header is "
            "required; required columns: sample_id and r1_path."
        ),
    )
    input_mode.add_argument(
        "--r1", type=Path,
        help="R1 or single-end FASTQ for single-sample mode.",
    )
    reference_source = required.add_mutually_exclusive_group(required=True)
    reference_source.add_argument(
        "-r", "--reference", type=Path,
        help="Non-redundant target-gene nucleotide FASTA.",
    )
    reference_source.add_argument(
        "--index", type=Path,
        help="Existing minimap2 .mmi index file, preferably built with -x sr.",
    )
    required.add_argument(
        "-p", "--prefix", type=Path, required=True,
        help="Output prefix for analysis mode; parent directory must exist.",
    )

    single = parser.add_argument_group("single-sample arguments")
    single.add_argument(
        "--r2", type=Path, default=None,
        help="R2 FASTQ; omit for single-end data (default: none).",
    )
    single.add_argument(
        "--sample-id", default=None,
        help="Sample ID; inferred from --r1 when omitted.",
    )
    optional = parser.add_argument_group("optional arguments")
    optional.add_argument(
        "--jobs", type=int, default=1,
        help="Samples processed concurrently (default: %(default)s).",
    )
    optional.add_argument(
        "--threads", type=int, default=8,
        help="Approximate CPU budget per sample (default: %(default)s).",
    )
    optional.add_argument(
        "--min-mapq", type=int, default=20,
        help="Minimum primary-alignment MAPQ (default: %(default)s).",
    )
    optional.add_argument(
        "--sort-memory", default="1G",
        help="samtools sort memory per thread (default: %(default)s).",
    )
    optional.add_argument(
        "--min-fragments", type=int, default=5,
        help="Minimum assigned fragments for detected=1 (default: %(default)s).",
    )
    optional.add_argument(
        "--min-breadth", type=float, default=50.0,
        help="Minimum coverage breadth percent for detected=1 (default: %(default)s).",
    )
    optional.add_argument(
        "--keep-bam", action="store_true",
        help="Keep coordinate-sorted BAM/BAI files in PREFIX.bam/.",
    )
    optional.add_argument(
        "--force", action="store_true",
        help="Delete previous outputs/checkpoint and recompute all samples.",
    )
    optional.add_argument(
        "--no-resume", action="store_true",
        help="Ignore completed samples in the checkpoint and recompute all samples.",
    )
    optional.add_argument(
        "--tmp-dir", type=Path, default=None,
        help="Parent directory for temporary files [default: system TMPDIR].",
    )
    optional.add_argument(
        "--help-input", "--help_input", action="store_true",
        help="Show detailed input-format requirements and exit.",
    )
    optional.add_argument(
        "--help-default", "--help_default", action="store_true",
        help="Show detailed default parameters and counting definitions and exit.",
    )
    return parser

def configure_logging(log_path: Path, append: bool) -> logging.Logger:
    logger = logging.getLogger(PROGRAM)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    rich_handler = RichHandler(
        console=CONSOLE,
        show_time=True,
        show_level=True,
        show_path=False,
        markup=True,
        rich_tracebacks=True,
    )
    rich_handler.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(rich_handler)

    file_handler = logging.FileHandler(
        log_path,
        mode="a" if append else "w",
        encoding="utf-8",
    )
    file_handler.setFormatter(
        logging.Formatter(
            "%(asctime)s\t%(levelname)s\t%(message)s",
            "%Y-%m-%d %H:%M:%S",
        )
    )
    logger.addHandler(file_handler)
    return logger


def require_programs(programs: Iterable[str]) -> None:
    missing = [program for program in programs if shutil.which(program) is None]
    if missing:
        raise RuntimeError("Required programs not found: " + ", ".join(missing))


def run_checked(
    command: list[str],
    *,
    stdout: Any = subprocess.PIPE,
    stderr: Any = subprocess.PIPE,
) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        command,
        text=True,
        stdout=stdout,
        stderr=stderr,
        check=False,
    )
    if completed.returncode != 0:
        stderr_text = completed.stderr if isinstance(completed.stderr, str) else ""
        raise RuntimeError(
            f"Command failed ({completed.returncode}): {shlex.join(command)}\n"
            f"{stderr_text[-5000:]}"
        )
    return completed


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_fasta(path: Path) -> tuple[list[str], dict[str, int]]:
    ids: list[str] = []
    lengths: dict[str, int] = {}
    current_id: str | None = None
    current_length = 0

    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    lengths[current_id] = current_length
                current_id = line[1:].split()[0]
                if not current_id:
                    raise ValueError(f"Empty FASTA ID in {path}")
                if current_id in lengths or current_id in ids:
                    raise ValueError(f"Duplicate FASTA ID: {current_id}")
                ids.append(current_id)
                current_length = 0
            else:
                if current_id is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                current_length += len(line)

    if current_id is not None:
        lengths[current_id] = current_length

    if not ids:
        raise ValueError(f"No FASTA records found: {path}")
    empty = [seq_id for seq_id in ids if lengths.get(seq_id, 0) <= 0]
    if empty:
        raise ValueError("Empty FASTA sequences: " + ", ".join(empty[:20]))
    return ids, lengths


def append_fasta(source: Path, output: TextIO) -> None:
    final_line = ""
    with open_text(source) as handle:
        for line in handle:
            output.write(line)
            final_line = line
    if final_line and not final_line.endswith("\n"):
        output.write("\n")


def resolve_fastq_path(value: str, manifest_dir: Path) -> str:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = manifest_dir / path
    return str(path.resolve())


def read_manifest(path: Path) -> list[Sample]:
    required = {"sample_id", "r1_path"}
    samples: list[Sample] = []
    seen: set[str] = set()
    manifest_dir = path.parent

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Input list has no header")
        missing = sorted(required.difference(reader.fieldnames))
        if missing:
            raise ValueError("Input list is missing columns: " + ", ".join(missing))

        for line_number, row in enumerate(reader, start=2):
            sample_id = (row.get("sample_id") or "").strip()
            if not sample_id:
                raise ValueError(f"Empty sample_id at input-list line {line_number}")
            if sample_id in seen:
                raise ValueError(f"Duplicate sample_id: {sample_id}")
            seen.add(sample_id)

            raw_r1 = (row.get("r1_path") or "").strip()
            raw_r2 = (row.get("r2_path") or "").strip()
            if not raw_r1:
                raise ValueError(f"Empty r1_path for sample {sample_id}")

            r1_path = resolve_fastq_path(raw_r1, manifest_dir)
            r2_path = raw_r2 if raw_r2 in MISSING_R2 else resolve_fastq_path(raw_r2, manifest_dir)
            samples.append(
                Sample(
                    sample_id=sample_id,
                    year=(row.get("year") or "").strip(),
                    month=(row.get("month") or "").strip(),
                    depth=(row.get("depth") or "").strip(),
                    r1_path=r1_path,
                    r2_path=r2_path,
                )
            )

    if not samples:
        raise ValueError("Input list contains no samples")
    return samples


def infer_sample_id(r1_path: Path) -> str:
    name = r1_path.name
    for suffix in (".gz", ".bz2", ".xz"):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    for suffix in (".fastq", ".fq", ".fasta", ".fa"):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    name = re.sub(r"(?i)([._-](R?1|READ1|FWD))$", "", name)
    name = name.strip("._-")
    if not name:
        raise ValueError("Cannot infer sample ID from --r1; provide --sample-id")
    return name


def make_single_sample(args: argparse.Namespace) -> Sample:
    assert args.r1 is not None
    r1 = args.r1.expanduser().resolve()
    r2 = args.r2.expanduser().resolve() if args.r2 is not None else None
    sample_id = (args.sample_id or infer_sample_id(r1)).strip()
    if not sample_id:
        raise ValueError("--sample-id must not be empty")
    return Sample(
        sample_id=sample_id,
        year="",
        month="",
        depth="",
        r1_path=str(r1),
        r2_path=str(r2) if r2 is not None else "",
    )

def validate_sample_files(samples: list[Sample]) -> None:
    errors: list[str] = []
    for sample in samples:
        r1 = Path(sample.r1_path)
        if not r1.is_file() or r1.stat().st_size <= 0:
            errors.append(f"{sample.sample_id}: invalid R1: {r1}")
        if sample.paired:
            r2 = Path(sample.r2_path)
            if not r2.is_file() or r2.stat().st_size <= 0:
                errors.append(f"{sample.sample_id}: invalid R2: {r2}")
    if errors:
        preview = "\n".join(errors[:20])
        extra = "" if len(errors) <= 20 else f"\n... and {len(errors) - 20} more"
        raise FileNotFoundError(preview + extra)


def sample_file_signature(sample: Sample) -> dict[str, Any]:
    signature: dict[str, Any] = {
        "sample_id": sample.sample_id,
        "paired": sample.paired,
    }
    for label, value in (
        ("r1", sample.r1_path),
        ("r2", sample.r2_path if sample.paired else ""),
    ):
        if not value:
            signature[label] = None
            continue
        stat = Path(value).stat()
        signature[label] = {
            "path": str(Path(value).resolve()),
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
        }
    return signature


def atomic_json_dump(data: dict[str, Any], path: Path) -> None:
    temp_path = Path(f"{path}.tmp")
    with temp_path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temp_path, path)


def atomic_tsv_write(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    temp_path = Path(f"{path}.tmp")
    with temp_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
    os.replace(temp_path, path)


def make_thread_plan(total_threads: int) -> ThreadPlan:
    return ThreadPlan(
        minimap2_threads=max(1, total_threads - 1),
        samtools_view_threads=1,
        samtools_sort_threads=max(1, total_threads - 1),
        samtools_collate_threads=max(1, min(2, total_threads // 4)),
        seqkit_threads=max(1, total_threads),
    )


def seqkit_counts(paths: list[str], threads: int) -> dict[str, int]:
    command = ["seqkit", "stats", "-T", "-j", str(max(1, threads)), *paths]
    completed = run_checked(command)
    rows = list(csv.DictReader(completed.stdout.splitlines(), delimiter="\t"))
    if len(rows) != len(paths):
        raise RuntimeError(
            f"seqkit stats returned {len(rows)} rows for {len(paths)} input files"
        )

    counts: dict[str, int] = {}
    for path, row in zip(paths, rows):
        raw_count = row.get("num_seqs")
        if raw_count is None:
            raise RuntimeError("seqkit stats -T output lacks num_seqs")
        counts[path] = int(raw_count.replace(",", ""))
    return counts


def cigar_aligned_query_length(cigar: str) -> int:
    if cigar == "*":
        return 0
    return sum(
        int(length)
        for length, operation in CIGAR_RE.findall(cigar)
        if operation in {"M", "I", "=", "X"}
    )


def parse_alignment(fields: list[str]) -> tuple[int, Alignment] | None:
    if len(fields) < 11:
        return None
    flag = int(fields[1])
    if flag & 0x4 or flag & 0x100 or flag & 0x800:
        return None
    reference = fields[2]
    if reference == "*":
        return None

    tags: dict[str, str] = {}
    for value in fields[11:]:
        pieces = value.split(":", 2)
        if len(pieces) == 3:
            tags[pieces[0]] = pieces[2]

    alignment = Alignment(
        reference=reference,
        as_score=int(tags.get("AS", "-1000000000")),
        mapq=int(fields[4]),
        aligned_query_length=cigar_aligned_query_length(fields[5]),
        nm=int(tags.get("NM", "1000000000")),
    )

    if flag & 0x40:
        mate = 1
    elif flag & 0x80:
        mate = 2
    else:
        mate = 0
    return mate, alignment


def choose_best(current: Alignment | None, candidate: Alignment) -> Alignment:
    if current is None or candidate.quality_key > current.quality_key:
        return candidate
    return current


def assign_group(
    mate_alignments: dict[int, Alignment],
    paired: bool,
    assigned_reads: Counter[str],
    assigned_fragments: Counter[str],
    qc: Counter[str],
) -> None:
    if not paired:
        candidates = list(mate_alignments.values())
        if not candidates:
            return
        best = max(candidates, key=lambda value: value.quality_key)
        assigned_reads[best.reference] += 1
        assigned_fragments[best.reference] += 1
        qc["single_end_assigned_fragments"] += 1
        return

    first = mate_alignments.get(1)
    second = mate_alignments.get(2)
    unknown = mate_alignments.get(0)

    if first is None and second is None:
        if unknown is not None:
            assigned_reads[unknown.reference] += 1
            assigned_fragments[unknown.reference] += 1
            qc["single_retained_mate_fragments"] += 1
        return

    if first is None or second is None:
        winner = first if first is not None else second
        assert winner is not None
        assigned_reads[winner.reference] += 1
        assigned_fragments[winner.reference] += 1
        qc["single_retained_mate_fragments"] += 1
        return

    if first.reference == second.reference:
        assigned_reads[first.reference] += 2
        assigned_fragments[first.reference] += 1
        qc["same_reference_pair_fragments"] += 1
        return

    qc["conflict_fragments"] += 1
    if first.quality_key > second.quality_key:
        winner = first
    elif second.quality_key > first.quality_key:
        winner = second
    else:
        qc["unresolved_tie_fragments"] += 1
        return

    assigned_reads[winner.reference] += 2
    assigned_fragments[winner.reference] += 1
    qc["resolved_conflict_fragments"] += 1


def count_assigned_fragments(
    filtered_bam: Path,
    paired: bool,
    collate_threads: int,
) -> tuple[Counter[str], Counter[str], Counter[str], int]:
    collate_command = [
        "samtools", "collate",
        "-@", str(max(1, collate_threads)),
        "-u", "-O",
        str(filtered_bam),
    ]
    collate = subprocess.Popen(
        collate_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert collate.stdout is not None
    view = subprocess.Popen(
        ["samtools", "view", "-"],
        stdin=collate.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    collate.stdout.close()
    assert view.stdout is not None

    assigned_reads: Counter[str] = Counter()
    assigned_fragments: Counter[str] = Counter()
    qc: Counter[str] = Counter()
    current_name: str | None = None
    mate_alignments: dict[int, Alignment] = {}
    retained_alignment_records = 0

    for line in view.stdout:
        fields = line.rstrip("\n").split("\t")
        parsed = parse_alignment(fields)
        if parsed is None:
            continue
        qname = fields[0]
        mate, alignment = parsed
        retained_alignment_records += 1

        if current_name is None:
            current_name = qname
        elif qname != current_name:
            assign_group(
                mate_alignments,
                paired,
                assigned_reads,
                assigned_fragments,
                qc,
            )
            current_name = qname
            mate_alignments = {}

        mate_alignments[mate] = choose_best(mate_alignments.get(mate), alignment)

    if current_name is not None:
        assign_group(
            mate_alignments,
            paired,
            assigned_reads,
            assigned_fragments,
            qc,
        )

    view_stderr = view.stderr.read() if view.stderr is not None else ""
    view_return = view.wait()
    collate_stderr_raw = collate.stderr.read() if collate.stderr is not None else b""
    collate_return = collate.wait()
    collate_stderr = collate_stderr_raw.decode("utf-8", errors="replace")

    if view_return != 0:
        raise RuntimeError(f"samtools view failed during fragment counting:\n{view_stderr[-4000:]}")
    if collate_return != 0:
        raise RuntimeError(f"samtools collate failed:\n{collate_stderr[-4000:]}")

    return assigned_reads, assigned_fragments, qc, retained_alignment_records


def read_coverage(coord_bam: Path) -> dict[str, dict[str, float | int]]:
    completed = run_checked(["samtools", "coverage", str(coord_bam)])
    result: dict[str, dict[str, float | int]] = {}
    for line in completed.stdout.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 9:
            continue
        result[fields[0]] = {
            "aligned_reads": int(fields[3]),
            "covered_bases": int(fields[4]),
            "breadth_pct": float(fields[5]),
            "mean_depth": float(fields[6]),
            "mean_baseq": float(fields[7]),
            "mean_mapq": float(fields[8]),
        }
    return result


def run_mapping_pipeline(
    command: list[str],
    filtered_bam: Path,
    min_mapq: int,
    view_threads: int,
    log_path: Path,
) -> None:
    with log_path.open("w", encoding="utf-8") as log_handle:
        mapper = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=log_handle,
        )
        assert mapper.stdout is not None
        view_command = [
            "samtools", "view",
            "-@", str(max(1, view_threads)),
            "-b",
            "-q", str(min_mapq),
            "-F", "2308",
            "-o", str(filtered_bam),
            "-",
        ]
        view = subprocess.Popen(
            view_command,
            stdin=mapper.stdout,
            stdout=subprocess.DEVNULL,
            stderr=log_handle,
        )
        mapper.stdout.close()
        view_return = view.wait()
        mapper_return = mapper.wait()

    if mapper_return != 0 or view_return != 0:
        details = log_path.read_text(encoding="utf-8", errors="replace")
        raise RuntimeError(
            "Mapping pipeline failed\n"
            f"minimap2 return code: {mapper_return}\n"
            f"samtools view return code: {view_return}\n"
            f"{details[-6000:]}"
        )


def process_sample(sample: Sample, config: WorkerConfig) -> dict[str, Any]:
    started = time.perf_counter()
    sample_tmp = Path(
        tempfile.mkdtemp(prefix=f"{sample.sample_id}.", dir=config.temp_root)
    )
    filtered_bam = sample_tmp / "retained.unsorted.bam"
    coord_bam = sample_tmp / "retained.coord.bam"
    mapper_log = sample_tmp / "minimap2.log"

    try:
        count_started = time.perf_counter()
        fastq_paths = [sample.r1_path]
        if sample.paired:
            fastq_paths.append(sample.r2_path)
        counts = seqkit_counts(fastq_paths, config.thread_plan.seqkit_threads)
        r1_count = counts[sample.r1_path]
        if sample.paired:
            r2_count = counts[sample.r2_path]
            if r1_count != r2_count:
                raise ValueError(
                    f"R1/R2 read-count mismatch: R1={r1_count}, R2={r2_count}"
                )
            total_clean_fragments = r1_count
            total_clean_reads = r1_count + r2_count
        else:
            total_clean_fragments = r1_count
            total_clean_reads = r1_count
        library_count_seconds = time.perf_counter() - count_started

        mapping_command = [
            "minimap2",
            "-ax", "sr",
            "--secondary=no",
            "-t", str(config.thread_plan.minimap2_threads),
            config.index_path,
            sample.r1_path,
        ]
        if sample.paired:
            mapping_command.append(sample.r2_path)

        mapping_started = time.perf_counter()
        run_mapping_pipeline(
            mapping_command,
            filtered_bam,
            config.min_mapq,
            config.thread_plan.samtools_view_threads,
            mapper_log,
        )
        mapping_seconds = time.perf_counter() - mapping_started

        fragment_started = time.perf_counter()
        assigned_reads, assigned_fragments, fragment_qc, retained_records = (
            count_assigned_fragments(
                filtered_bam,
                sample.paired,
                config.thread_plan.samtools_collate_threads,
            )
        )
        fragment_counting_seconds = time.perf_counter() - fragment_started

        coverage_started = time.perf_counter()
        run_checked([
            "samtools", "sort",
            "-@", str(config.thread_plan.samtools_sort_threads),
            "-m", config.sort_memory,
            "-o", str(coord_bam),
            str(filtered_bam),
        ])
        coverage = read_coverage(coord_bam)
        coverage_seconds = time.perf_counter() - coverage_started

        if config.keep_bam:
            assert config.bam_dir is not None
            bam_dir = Path(config.bam_dir)
            destination_bam = bam_dir / f"{sample.sample_id}.q{config.min_mapq}.bam"
            run_checked([
                "samtools", "index",
                "-@", str(config.thread_plan.samtools_sort_threads),
                str(coord_bam),
            ])
            shutil.move(str(coord_bam), destination_bam)
            shutil.move(str(coord_bam) + ".bai", Path(str(destination_bam) + ".bai"))

        gene_rows: list[dict[str, Any]] = []
        for gene_id in config.target_ids:
            gene_length = config.target_lengths[gene_id]
            reads = int(assigned_reads.get(gene_id, 0))
            fragments = int(assigned_fragments.get(gene_id, 0))
            cov = coverage.get(gene_id, {})
            read_rpm = (
                1_000_000.0 * reads / total_clean_reads
                if total_clean_reads > 0 else 0.0
            )
            fragment_fpm = (
                1_000_000.0 * fragments / total_clean_fragments
                if total_clean_fragments > 0 else 0.0
            )
            rpkm = (
                1_000_000_000.0 * reads / (gene_length * total_clean_reads)
                if total_clean_reads > 0 else 0.0
            )
            fpkm = (
                1_000_000_000.0 * fragments /
                (gene_length * total_clean_fragments)
                if total_clean_fragments > 0 else 0.0
            )
            gene_rows.append({
                "gene_id": gene_id,
                "gene_length": gene_length,
                "total_clean_reads": total_clean_reads,
                "total_clean_fragments": total_clean_fragments,
                "aligned_reads": int(cov.get("aligned_reads", 0)),
                "assigned_reads": reads,
                "assigned_fragments": fragments,
                "covered_bases": int(cov.get("covered_bases", 0)),
                "breadth_pct": float(cov.get("breadth_pct", 0.0)),
                "mean_depth": float(cov.get("mean_depth", 0.0)),
                "mean_baseq": float(cov.get("mean_baseq", 0.0)),
                "mean_mapq": float(cov.get("mean_mapq", 0.0)),
                "read_rpm": read_rpm,
                "fragment_fpm": fragment_fpm,
                "rpkm": rpkm,
                "fpkm": fpkm,
            })

        target_assigned_reads = sum(
            assigned_reads.get(gene_id, 0) for gene_id in config.target_set
        )
        target_assigned_fragments = sum(
            assigned_fragments.get(gene_id, 0) for gene_id in config.target_set
        )
        target_aligned_reads = sum(
            int(coverage.get(gene_id, {}).get("aligned_reads", 0))
            for gene_id in config.target_set
        )
        all_aligned_reads = sum(
            int(values.get("aligned_reads", 0)) for values in coverage.values()
        )
        elapsed_seconds = time.perf_counter() - started

        qc = {
            "library_type": "PE" if sample.paired else "SE",
            "total_clean_reads": total_clean_reads,
            "total_clean_fragments": total_clean_fragments,
            "retained_alignment_records": retained_records,
            "all_reference_aligned_reads": all_aligned_reads,
            "target_aligned_reads": target_aligned_reads,
            "target_assigned_reads": target_assigned_reads,
            "target_assigned_fragments": target_assigned_fragments,
            "same_reference_pair_fragments": int(fragment_qc["same_reference_pair_fragments"]),
            "single_retained_mate_fragments": int(fragment_qc["single_retained_mate_fragments"]),
            "conflict_fragments": int(fragment_qc["conflict_fragments"]),
            "resolved_conflict_fragments": int(fragment_qc["resolved_conflict_fragments"]),
            "unresolved_tie_fragments": int(fragment_qc["unresolved_tie_fragments"]),
            "target_assigned_read_pct": (
                100.0 * target_assigned_reads / total_clean_reads
                if total_clean_reads > 0 else 0.0
            ),
            "library_count_seconds": library_count_seconds,
            "mapping_seconds": mapping_seconds,
            "fragment_counting_seconds": fragment_counting_seconds,
            "coverage_seconds": coverage_seconds,
            "elapsed_seconds": elapsed_seconds,
            "clean_reads_per_mapping_second": (
                total_clean_reads / mapping_seconds if mapping_seconds > 0 else 0.0
            ),
            "clean_fragments_per_second": (
                total_clean_fragments / elapsed_seconds if elapsed_seconds > 0 else 0.0
            ),
            "minimap2_threads": config.thread_plan.minimap2_threads,
            "samtools_sort_threads": config.thread_plan.samtools_sort_threads,
            "seqkit_threads": config.thread_plan.seqkit_threads,
        }
        return {
            "sample_id": sample.sample_id,
            "gene_rows": gene_rows,
            "qc": qc,
        }
    finally:
        shutil.rmtree(sample_tmp, ignore_errors=True)


def valid_completed_entry(
    entry: dict[str, Any],
    expected_signature: dict[str, Any],
    target_ids: list[str],
) -> bool:
    if entry.get("sample_signature") != expected_signature:
        return False
    result = entry.get("result")
    if not isinstance(result, dict):
        return False
    rows = result.get("gene_rows")
    if not isinstance(rows, list) or len(rows) != len(target_ids):
        return False
    row_ids = [row.get("gene_id") for row in rows if isinstance(row, dict)]
    return row_ids == target_ids and isinstance(result.get("qc"), dict)


def write_outputs(
    prefix: Path,
    samples: list[Sample],
    target_ids: list[str],
    target_lengths: dict[str, int],
    completed_results: dict[str, dict[str, Any]],
    min_fragments: int,
    min_breadth: float,
) -> None:
    long_rows: list[dict[str, Any]] = []
    qc_rows: list[dict[str, Any]] = []

    for sample in samples:
        result = completed_results.get(sample.sample_id)
        if result is None:
            continue
        for raw_row in result["gene_rows"]:
            row = {
                "sample_id": sample.sample_id,
                "year": sample.year,
                "month": sample.month,
                "depth": sample.depth,
                "library_type": result["qc"]["library_type"],
                **raw_row,
            }
            row["detected"] = int(
                int(row["assigned_fragments"]) >= min_fragments
                and float(row["breadth_pct"]) >= min_breadth
            )
            long_rows.append(row)
        qc_rows.append({
            "sample_id": sample.sample_id,
            "year": sample.year,
            "month": sample.month,
            "depth": sample.depth,
            **result["qc"],
        })

    long_fields = [
        "sample_id", "year", "month", "depth", "library_type",
        "gene_id", "gene_length",
        "total_clean_reads", "total_clean_fragments",
        "aligned_reads", "assigned_reads", "assigned_fragments",
        "covered_bases", "breadth_pct", "mean_depth",
        "mean_baseq", "mean_mapq",
        "read_rpm", "fragment_fpm", "rpkm", "fpkm", "detected",
    ]
    qc_fields = [
        "sample_id", "year", "month", "depth", "library_type",
        "total_clean_reads", "total_clean_fragments",
        "retained_alignment_records", "all_reference_aligned_reads",
        "target_aligned_reads", "target_assigned_reads",
        "target_assigned_fragments", "same_reference_pair_fragments",
        "single_retained_mate_fragments", "conflict_fragments",
        "resolved_conflict_fragments", "unresolved_tie_fragments",
        "target_assigned_read_pct", "library_count_seconds",
        "mapping_seconds", "fragment_counting_seconds", "coverage_seconds",
        "elapsed_seconds", "clean_reads_per_mapping_second",
        "clean_fragments_per_second", "minimap2_threads",
        "samtools_sort_threads", "seqkit_threads",
    ]

    atomic_tsv_write(
        Path(f"{prefix}.meta_gene_abundance.long.tsv"),
        long_rows,
        long_fields,
    )
    atomic_tsv_write(
        Path(f"{prefix}.sample_qc.tsv"),
        qc_rows,
        qc_fields,
    )
    atomic_tsv_write(
        Path(f"{prefix}.reference.tsv"),
        [
            {"gene_id": gene_id, "gene_length": target_lengths[gene_id]}
            for gene_id in target_ids
        ],
        ["gene_id", "gene_length"],
    )
    Path(f"{prefix}.meta_gene_abundance.rpkm.tsv").unlink(missing_ok=True)
    Path(f"{prefix}.meta_gene_abundance.fpkm.tsv").unlink(missing_ok=True)


def display_parameters(
    args: argparse.Namespace,
    samples: list[Sample],
    target_count: int,
    plan: ThreadPlan,
    index_description: str,
) -> None:
    table = Table(title="Analysis parameters", show_header=False, box=None)
    table.add_column("Parameter", style="bold cyan", no_wrap=True)
    table.add_column("Value", style="white")
    table.add_row("Mapper", "minimap2 -ax sr")
    table.add_row(
        "Input mode",
        f"inlist: {args.manifest}" if args.manifest else f"single sample: {samples[0].sample_id}",
    )
    table.add_row("Reference", str(args.reference) if args.reference else "from --index header")
    table.add_row("Index", index_description)
    table.add_row("Prefix", str(args.prefix))
    table.add_row("Samples", str(len(samples)))
    table.add_row("Target genes", str(target_count))
    table.add_row("Parallel samples", str(args.jobs))
    table.add_row("CPU budget/sample", str(args.threads))
    table.add_row("minimap2 threads", str(plan.minimap2_threads))
    table.add_row("samtools sort threads", str(plan.samtools_sort_threads))
    table.add_row("Minimum MAPQ", str(args.min_mapq))
    table.add_row("Conflict rule", CONFLICT_RULE)
    table.add_row("Detection", f"fragments ≥ {args.min_fragments}; breadth ≥ {args.min_breadth:g}%")
    table.add_row("Keep BAM", str(args.keep_bam))
    CONSOLE.print(table)

def output_paths(prefix: Path) -> list[Path]:
    return [
        Path(f"{prefix}.meta_gene_abundance.long.tsv"),
        Path(f"{prefix}.sample_qc.tsv"),
        Path(f"{prefix}.reference.tsv"),
        Path(f"{prefix}.failed.tsv"),
        Path(f"{prefix}.run.log"),
        Path(f"{prefix}.state.json"),
    ]



def index_files_exist(index_path: Path) -> bool:
    return index_path.is_file() and index_path.stat().st_size > 0


def index_metadata_path(index_path: Path) -> Path:
    return Path(f"{index_path}.gene_index.json")


def index_fingerprint(index_path: Path) -> dict[str, Any]:
    return {
        "path": str(index_path.resolve()),
        "sha256": file_sha256(index_path),
    }


def parse_sam_sequence_dictionary(header_text: str) -> tuple[list[str], dict[str, int]]:
    target_ids: list[str] = []
    target_lengths: dict[str, int] = {}

    for line in header_text.splitlines():
        if not line.startswith("@SQ\t"):
            continue
        tags: dict[str, str] = {}
        for field in line.split("\t")[1:]:
            if len(field) >= 4 and field[2] == ":":
                tags[field[:2]] = field[3:]
        gene_id = tags.get("SN", "")
        raw_length = tags.get("LN", "")
        if not gene_id:
            raise RuntimeError("minimap2 SAM header contains @SQ without SN")
        if gene_id in target_lengths:
            raise RuntimeError(f"Duplicate target ID in minimap2 index: {gene_id}")
        try:
            length = int(raw_length)
        except ValueError as exc:
            raise RuntimeError(
                f"Invalid target length in minimap2 SAM header: {gene_id}={raw_length}"
            ) from exc
        if length <= 0:
            raise RuntimeError(f"Non-positive target length in minimap2 SAM header: {gene_id}")
        target_ids.append(gene_id)
        target_lengths[gene_id] = length

    if not target_ids:
        raise RuntimeError("No @SQ records found in minimap2 SAM header")
    return target_ids, target_lengths


def load_index_reference_info(index_path: Path) -> tuple[list[str], dict[str, int]]:
    completed = run_checked([
        "minimap2",
        "-ax", "sr",
        "--secondary=no",
        str(index_path),
        "/dev/null",
    ])
    return parse_sam_sequence_dictionary(completed.stdout)


def build_minimap2_index(
    reference: Path,
    index_path: Path,
    threads: int,
    logger: logging.Logger | None,
) -> None:
    if not index_path.parent.is_dir():
        raise FileNotFoundError(
            f"Index parent directory does not exist: {index_path.parent}"
        )

    metadata_path = index_metadata_path(index_path)
    index_path.unlink(missing_ok=True)
    metadata_path.unlink(missing_ok=True)

    command = [
        "minimap2",
        "-x", "sr",
        "-t", str(max(1, threads)),
        "-d", str(index_path),
        str(reference),
    ]
    if logger is not None:
        logger.info("Building minimap2 index: %s", index_path)
        logger.info("Index command: %s", shlex.join(command))

    started = time.perf_counter()
    run_checked(command)
    if not index_files_exist(index_path):
        raise RuntimeError(
            "minimap2 index returned successfully but expected .mmi file "
            f"was not created: {index_path}"
        )

    target_ids, target_lengths = read_fasta(reference)
    metadata = {
        "program": PROGRAM,
        "mapper": "minimap2",
        "reference": str(reference.resolve()),
        "reference_sha256": file_sha256(reference),
        "target_ids": target_ids,
        "target_lengths": target_lengths,
        "build_seconds": time.perf_counter() - started,
        "index": index_fingerprint(index_path),
    }
    atomic_json_dump(metadata, metadata_path)
    if logger is not None:
        logger.info("Index completed in %.2f s", metadata["build_seconds"])

def main() -> int:
    parser = build_parser()

    if "--help-input" in sys.argv[1:] or "--help_input" in sys.argv[1:]:
        CONSOLE.print(INPUT_HELP)
        return 0
    if "--help-default" in sys.argv[1:] or "--help_default" in sys.argv[1:]:
        CONSOLE.print(DEFAULT_HELP)
        return 0
    if any(value in {"-m", "--manifest"} for value in sys.argv[1:]):
        parser.error("Use -i/-inlist/--inlist instead of -m/--manifest.")

    args = parser.parse_args()
    started = time.perf_counter()

    if args.reference is not None:
        args.reference = args.reference.expanduser().resolve()
    if args.index is not None:
        args.index = args.index.expanduser().resolve()
    if args.tmp_dir is not None:
        args.tmp_dir = args.tmp_dir.expanduser().resolve()

    if args.jobs < 1:
        parser.error("--jobs must be >= 1")
    if args.threads < 1:
        parser.error("--threads must be >= 1")
    if not 0 <= args.min_mapq <= 255:
        parser.error("--min-mapq must be between 0 and 255")
    if args.min_fragments < 0:
        parser.error("--min-fragments must be >= 0")
    if not 0 <= args.min_breadth <= 100:
        parser.error("--min-breadth must be between 0 and 100")
    if args.reference is not None and not args.reference.is_file():
        parser.error(f"Reference not found: {args.reference}")
    if args.index is not None and not index_files_exist(args.index):
        parser.error(f"Minimap2 index file not found for --index: {args.index}")
    if args.tmp_dir is not None and not args.tmp_dir.is_dir():
        parser.error(f"Temporary directory does not exist: {args.tmp_dir}")

    single_only_values = [args.r2, args.sample_id]
    if args.manifest is not None and any(value not in (None, "") for value in single_only_values):
        parser.error("--r2/--sample-id are only valid with --r1")
    if args.r2 is not None and args.r1 is None:
        parser.error("--r2 requires --r1")

    args.prefix = args.prefix.expanduser().resolve()
    if not args.prefix.parent.is_dir():
        parser.error(f"Prefix parent directory does not exist: {args.prefix.parent}")

    require_programs(["minimap2", "samtools", "seqkit"])

    temp_parent = str(args.tmp_dir) if args.tmp_dir is not None else None
    reference_sha256: str | None = None
    if args.reference is not None:
        target_ids, target_lengths = read_fasta(args.reference)
        reference_sha256 = file_sha256(args.reference)
    else:
        assert args.index is not None
        try:
            target_ids, target_lengths = load_index_reference_info(args.index)
        except RuntimeError as exc:
            parser.error(str(exc))

    assert args.prefix is not None
    if args.manifest is not None:
        args.manifest = args.manifest.expanduser().resolve()
        if not args.manifest.is_file():
            parser.error(f"Input list not found: {args.manifest}")
        samples = read_manifest(args.manifest)
    else:
        assert args.r1 is not None
        samples = [make_single_sample(args)]
    validate_sample_files(samples)

    paths = output_paths(args.prefix)
    log_path = Path(f"{args.prefix}.run.log")
    state_path = Path(f"{args.prefix}.state.json")

    if args.force:
        for path in paths:
            path.unlink(missing_ok=True)
        bam_dir = Path(f"{args.prefix}.bam")
        if bam_dir.exists():
            shutil.rmtree(bam_dir)

    if not state_path.exists() and not args.force:
        preexisting = [
            path for path in paths
            if path.exists() and path not in {log_path, state_path}
        ]
        if preexisting:
            parser.error(
                "Output files exist without a matching checkpoint. Use --force: "
                + ", ".join(str(path) for path in preexisting[:5])
            )

    logger = configure_logging(log_path, append=log_path.exists() and not args.force)

    try:
        plan = make_thread_plan(args.threads)
        index_description = (
            f"existing: {args.index}"
            if args.index is not None
            else "temporary, built automatically"
        )
        display_parameters(
            args,
            samples,
            len(target_ids),
            plan,
            index_description,
        )
        logger.info("Starting %s", PROGRAM)
        logger.info("Mapper: minimap2; preset: sr; conflict rule: %s", CONFLICT_RULE)
        if args.threads < 3:
            logger.warning(
                "--threads=%d is small; minimap2 and samtools run in a pipeline, "
                "so actual CPU use may briefly exceed this value.",
                args.threads,
            )

        mapping_signature = {
            "schema": CHECKPOINT_SCHEMA,
            "program": PROGRAM,
            "mapper": "minimap2",
            "mapper_preset": "-ax sr",
            "secondary_output": False,
            "suppress_unmapped": True,
            "conflict_rule": CONFLICT_RULE,
            "reference_sha256": reference_sha256,
            "target_ids": target_ids,
            "min_mapq": args.min_mapq,
            "index": index_fingerprint(args.index) if args.index else None,
        }

        state: dict[str, Any] = {
            "schema": CHECKPOINT_SCHEMA,
            "program": PROGRAM,
            "mapping_signature": mapping_signature,
            "completed": {},
        }

        if state_path.exists() and not args.force and not args.no_resume:
            with state_path.open("r", encoding="utf-8") as handle:
                previous = json.load(handle)
            if previous.get("mapping_signature") != mapping_signature:
                raise RuntimeError(
                    "Existing checkpoint was created with different references, "
                    "index or mapping settings. Use --force to recompute."
                )
            state = previous
            logger.info("Loaded checkpoint with %d stored samples", len(state.get("completed", {})))
        elif args.no_resume:
            logger.info("Resume disabled; all samples will be recomputed")

        completed_entries = state.setdefault("completed", {})
        sample_signatures = {sample.sample_id: sample_file_signature(sample) for sample in samples}

        reusable_results: dict[str, dict[str, Any]] = {}
        pending: list[Sample] = []
        for sample in samples:
            entry = completed_entries.get(sample.sample_id)
            if (
                not args.no_resume
                and isinstance(entry, dict)
                and valid_completed_entry(entry, sample_signatures[sample.sample_id], target_ids)
            ):
                reusable_results[sample.sample_id] = entry["result"]
            else:
                completed_entries.pop(sample.sample_id, None)
                pending.append(sample)

        state["completed"] = completed_entries
        atomic_json_dump(state, state_path)
        if reusable_results:
            logger.info("Resuming %d completed samples", len(reusable_results))

        if not pending:
            logger.info("All samples are complete; refreshing output tables")
            write_outputs(
                args.prefix, samples, target_ids, target_lengths,
                reusable_results, args.min_fragments, args.min_breadth,
            )
            Path(f"{args.prefix}.failed.tsv").unlink(missing_ok=True)
            elapsed = time.perf_counter() - started
            total_reads = sum(result["qc"]["total_clean_reads"] for result in reusable_results.values())
            logger.info(
                "Refreshed %d/%d samples in %.2f s; %.1f stored clean reads/s",
                len(reusable_results), len(samples), elapsed,
                total_reads / elapsed if elapsed > 0 else 0.0,
            )
            CONSOLE.print("\n[bold green]Completed successfully.[/bold green]")
            CONSOLE.print(f"Long table  : [cyan]{args.prefix}.meta_gene_abundance.long.tsv[/cyan]")
            CONSOLE.print(f"Sample QC   : [cyan]{args.prefix}.sample_qc.tsv[/cyan]")
            return 0

        failures: list[dict[str, str]] = []
        with tempfile.TemporaryDirectory(prefix="gene_minimap2.", dir=temp_parent) as temp_root:
            temp_root_path = Path(temp_root)
            if args.index is None:
                assert args.reference is not None
                index_path = temp_root_path / "reference.mmi"
                build_minimap2_index(args.reference, index_path, args.threads, logger)
            else:
                index_path = args.index

            bam_dir: Path | None = None
            if args.keep_bam:
                bam_dir = Path(f"{args.prefix}.bam")
                bam_dir.mkdir(parents=False, exist_ok=True)

            config = WorkerConfig(
                index_path=str(index_path),
                temp_root=str(temp_root_path),
                target_ids=tuple(target_ids),
                target_lengths=target_lengths,
                target_set=frozenset(target_ids),
                min_mapq=args.min_mapq,
                sort_memory=args.sort_memory,
                keep_bam=args.keep_bam,
                bam_dir=str(bam_dir) if bam_dir else None,
                thread_plan=plan,
            )

            new_results: dict[str, dict[str, Any]] = {}
            progress = Progress(
                SpinnerColumn(),
                TextColumn("[bold blue]{task.description}"),
                BarColumn(),
                MofNCompleteColumn(),
                TimeElapsedColumn(),
                TimeRemainingColumn(),
                console=CONSOLE,
            )
            with progress:
                task_id = progress.add_task("Mapping samples", total=len(pending))
                with ThreadPoolExecutor(max_workers=args.jobs) as executor:
                    future_to_sample = {
                        executor.submit(process_sample, sample, config): sample
                        for sample in pending
                    }
                    for future in as_completed(future_to_sample):
                        sample = future_to_sample[future]
                        try:
                            result = future.result()
                            new_results[sample.sample_id] = result
                            completed_entries[sample.sample_id] = {
                                "sample_signature": sample_signatures[sample.sample_id],
                                "result": result,
                            }
                            atomic_json_dump(state, state_path)
                            qc = result["qc"]
                            logger.info(
                                "[green]Completed[/green] %s: target fragments=%d; "
                                "resolved conflicts=%d; %.1f clean reads/s",
                                sample.sample_id,
                                qc["target_assigned_fragments"],
                                qc["resolved_conflict_fragments"],
                                qc["clean_reads_per_mapping_second"],
                            )
                        except Exception as exc:
                            failures.append({
                                "sample_id": sample.sample_id,
                                "error": str(exc).replace("\n", " | "),
                            })
                            logger.error("[red]Failed[/red] %s: %s", sample.sample_id, exc)
                        finally:
                            progress.advance(task_id)

            completed_results = {**reusable_results, **new_results}
            for sample in samples:
                entry = completed_entries.get(sample.sample_id)
                if isinstance(entry, dict) and valid_completed_entry(
                    entry, sample_signatures[sample.sample_id], target_ids,
                ):
                    completed_results[sample.sample_id] = entry["result"]

            write_outputs(
                args.prefix, samples, target_ids, target_lengths,
                completed_results, args.min_fragments, args.min_breadth,
            )

            failed_path = Path(f"{args.prefix}.failed.tsv")
            if failures:
                atomic_tsv_write(failed_path, failures, ["sample_id", "error"])
            else:
                failed_path.unlink(missing_ok=True)

            elapsed = time.perf_counter() - started
            total_reads = sum(result["qc"]["total_clean_reads"] for result in completed_results.values())
            logger.info(
                "Finished %d/%d samples in %.2f s; overall %.1f clean reads/s",
                len(completed_results), len(samples), elapsed,
                total_reads / elapsed if elapsed > 0 else 0.0,
            )

            if failures:
                logger.error(
                    "%d sample(s) failed. Correct the inputs and rerun; successful "
                    "samples will resume from the checkpoint.",
                    len(failures),
                )
                return 1

        CONSOLE.print("\n[bold green]Completed successfully.[/bold green]")
        CONSOLE.print(f"Long table  : [cyan]{args.prefix}.meta_gene_abundance.long.tsv[/cyan]")
        CONSOLE.print(f"Sample QC   : [cyan]{args.prefix}.sample_qc.tsv[/cyan]")
        return 0

    except Exception as exc:
        logger.exception("Analysis failed: %s", exc)
        return 1

if __name__ == "__main__":
    raise SystemExit(main())
