#!/usr/bin/env python3
"""Convert a MetaPhlAn 4 profile into Krona ktImportText input.

Research use only. Not for clinical or diagnostic use.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


RANK_PREFIXES = ("k__", "p__", "c__", "o__", "f__", "g__", "s__", "t__")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert MetaPhlAn 4 profile TSV to Krona text.")
    parser.add_argument("--profile", required=True, help="MetaPhlAn 4 profile (.txt).")
    parser.add_argument("--output", required=True, help="Krona tab-delimited text output.")
    return parser.parse_args()


def last_rank_token(clade: str) -> str:
    tokens = [token for token in clade.split("|") if token]
    return tokens[-1] if tokens else ""


def parse_lineage(clade: str) -> list[str]:
    if clade.upper().startswith("UNCLASSIFIED") or clade.upper() in {"UNKNOWN", "UNCLASSIFIED"}:
        return ["UNCLASSIFIED"]
    lineage: list[str] = []
    for token in clade.split("|"):
        stripped = token
        for prefix in RANK_PREFIXES:
            if token.startswith(prefix):
                stripped = token[len(prefix) :]
                break
        if stripped:
            lineage.append(stripped)
    return lineage or ["UNCLASSIFIED"]


def detect_abundance_index(header: list[str]) -> int:
    lowered = [col.strip().lower() for col in header]
    for name in (
        "estimated_number_of_reads_from_the_clade",
        "relative_abundance",
        "abundance",
    ):
        if name in lowered:
            return lowered.index(name)
    return 2 if len(header) > 2 else 1


def read_profile(path: Path) -> list[tuple[str, float]]:
    rows: list[tuple[str, float]] = []
    header: list[str] | None = None
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                stripped = line.lstrip("#")
                if stripped.startswith("clade_name"):
                    header = stripped.split("\t")
                continue
            parts = line.split("\t")
            if parts[0] == "clade_name":
                header = parts
                continue
            if header is None:
                header = ["clade_name", "NCBI_tax_id", "relative_abundance"]
            abundance_idx = detect_abundance_index(header)
            if abundance_idx >= len(parts):
                continue
            try:
                abundance = float(parts[abundance_idx])
            except ValueError:
                continue
            if abundance <= 0:
                continue
            rows.append((parts[0], abundance))
    return rows


def select_leaves(rows: list[tuple[str, float]]) -> list[tuple[float, list[str]]]:
    sgb: list[tuple[float, list[str]]] = []
    species: list[tuple[float, list[str]]] = []
    unclassified: list[tuple[float, list[str]]] = []
    other: list[tuple[float, list[str]]] = []
    for clade, abundance in rows:
        if clade.upper().startswith("UNCLASSIFIED") or clade.upper() in {"UNKNOWN", "UNCLASSIFIED"}:
            unclassified.append((abundance, ["UNCLASSIFIED"]))
            continue
        last = last_rank_token(clade)
        lineage = parse_lineage(clade)
        if last.startswith("t__"):
            sgb.append((abundance, lineage))
        elif last.startswith("s__"):
            species.append((abundance, lineage))
        else:
            other.append((abundance, lineage))
    if sgb:
        return unclassified + sgb
    if species:
        return unclassified + species
    return unclassified + other


def write_krona(path: Path, rows: list[tuple[float, list[str]]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for abundance, lineage in rows:
            writer.writerow([f"{abundance:g}"] + lineage)


def main() -> None:
    args = parse_args()
    leaves = select_leaves(read_profile(Path(args.profile)))
    write_krona(Path(args.output), leaves)


if __name__ == "__main__":
    main()
