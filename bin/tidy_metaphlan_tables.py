#!/usr/bin/env python3
"""Filter MetaPhlAn 4 merged tables and summarise Bowtie2 host-removal logs.

Research use only. Not for clinical or diagnostic use.
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Tidy MetaPhlAn 4 merged abundance tables and Bowtie2 host-removal logs."
    )
    parser.add_argument("--merged", required=True, help="Merged MetaPhlAn abundance table (TSV).")
    parser.add_argument("--species-out", required=True, help="Species-level table (s__, excluding t__).")
    parser.add_argument("--sgb-out", required=True, help="SGB-level table (t__).")
    parser.add_argument(
        "--bowtie2-logs",
        nargs="*",
        default=[],
        help="Optional Bowtie2 alignment logs to summarise.",
    )
    parser.add_argument("--host-summary-out", help="Host-removal summary TSV (written when logs are provided).")
    return parser.parse_args()


def _is_comment(line: str) -> bool:
    return line.startswith("#")


def _last_rank_token(clade: str) -> str:
    tokens = [token for token in clade.split("|") if token]
    return tokens[-1] if tokens else ""


def is_species_clade(clade: str) -> bool:
    """True for species-level rows (last rank is s__), excluding SGB (t__) rows."""
    last = _last_rank_token(clade)
    return last.startswith("s__")


def is_sgb_clade(clade: str) -> bool:
    """True for SGB-level rows (last rank is t__)."""
    last = _last_rank_token(clade)
    return last.startswith("t__")


def read_merged_table(path: Path) -> tuple[list[str], list[list[str]]]:
    with path.open() as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip() and not _is_comment(line.strip())]
    if not lines:
        return ["clade_name"], []

    reader = csv.reader(lines, delimiter="\t")
    rows = list(reader)
    header = rows[0]
    if header and header[0] in {"clade_name", "clade"}:
        data = rows[1:]
    else:
        header = ["clade_name"] + [f"sample_{i}" for i in range(1, len(rows[0]))]
        data = rows
    return header, data


def write_table(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def parse_bowtie2_log(path: Path) -> dict[str, str]:
    text = path.read_text(errors="replace")
    sample = path.name
    for suffix in (".bowtie2.log", ".log"):
        if sample.endswith(suffix):
            sample = sample[: -len(suffix)]
            break

    total_reads = ""
    total_match = re.search(r"^(\d+) reads;", text, re.MULTILINE)
    if total_match:
        total_reads = total_match.group(1)

    rate = ""
    rate_match = re.search(r"([0-9.]+)% overall alignment rate", text)
    if rate_match:
        rate = rate_match.group(1)

    aligned_reads = ""
    unaligned_reads = ""
    if total_reads and rate:
        total = int(total_reads)
        aligned = int(round(total * float(rate) / 100.0))
        aligned = min(max(aligned, 0), total)
        aligned_reads = str(aligned)
        unaligned_reads = str(total - aligned)

    return {
        "sample": sample,
        "total_reads": total_reads,
        "aligned_reads": aligned_reads,
        "unaligned_reads": unaligned_reads,
        "overall_alignment_rate": rate,
    }


def write_host_summary(path: Path, logs: list[Path]) -> None:
    header = [
        "sample",
        "total_reads",
        "aligned_reads",
        "unaligned_reads",
        "overall_alignment_rate",
    ]
    rows = [list(parse_bowtie2_log(log).values()) for log in logs]
    write_table(path, header, rows)


def main() -> None:
    args = parse_args()
    header, data = read_merged_table(Path(args.merged))

    species_rows = [row for row in data if row and is_species_clade(row[0])]
    sgb_rows = [row for row in data if row and is_sgb_clade(row[0])]

    write_table(Path(args.species_out), header, species_rows)
    write_table(Path(args.sgb_out), header, sgb_rows)

    logs = [Path(item) for item in args.bowtie2_logs if item and Path(item).exists()]
    if args.host_summary_out:
        write_host_summary(Path(args.host_summary_out), logs)


if __name__ == "__main__":
    main()
