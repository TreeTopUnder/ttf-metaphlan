#!/usr/bin/env python3
"""Select StrainPhlAn 4 SGB clades from a merged SGB table.

Research use only. Not for clinical or diagnostic use.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

CLADE_RE = re.compile(r"^t__SGB[A-Za-z0-9_]+$")
SUMMARY_HEADER = [
    "clade",
    "samples_present",
    "mean_abundance",
    "selection_mode",
    "status",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Select StrainPhlAn SGB clades from a merged SGB abundance table."
    )
    parser.add_argument("--sgb-table", required=True, help="Merged SGB table TSV (clade_name + samples).")
    parser.add_argument(
        "--clades",
        default="",
        help="Optional comma-separated SGB clade IDs (t__SGB...). Empty/null selects auto top-N.",
    )
    parser.add_argument("--top-n", type=int, default=5, help="Max auto-selected SGB clades.")
    parser.add_argument(
        "--min-samples",
        type=int,
        default=4,
        help="Minimum samples with abundance > 0 for auto-selection / tree building.",
    )
    parser.add_argument("--summary-out", required=True, help="Summary TSV path.")
    parser.add_argument("--clades-out", required=True, help="Plain-text list of selected clade IDs.")
    return parser.parse_args()


def sgb_token(clade: str) -> str | None:
    """Return the t__SGB token from a full lineage or a bare clade ID."""
    if not clade:
        return None
    token = clade.strip()
    if token.startswith("#") or token in {"clade_name", "clade", "UNCLASSIFIED"}:
        return None
    last = token.split("|")[-1].strip()
    if last.startswith("t__"):
        return last
    if token.startswith("t__"):
        return token
    return None


def is_valid_clade(clade: str) -> bool:
    return bool(clade) and CLADE_RE.match(clade) is not None


def parse_clade_list(raw: str | None) -> list[str]:
    if raw is None:
        return []
    text = str(raw).strip()
    if not text or text.lower() in {"null", "none"}:
        return []
    clades: list[str] = []
    seen: set[str] = set()
    for part in text.split(","):
        token = part.strip()
        if not token:
            continue
        if token not in seen:
            clades.append(token)
            seen.add(token)
    return clades


def _to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def read_sgb_stats(path: Path) -> dict[str, dict[str, float | int]]:
    """Map t__SGB ID -> samples_present, mean_abundance from a merged table."""
    stats: dict[str, dict[str, float | int]] = {}
    if not path.exists():
        return stats
    with path.open() as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip() and not line.startswith("#")]
    if not lines:
        return stats
    reader = csv.reader(lines, delimiter="\t")
    rows = list(reader)
    header = rows[0]
    data = rows[1:]
    if header and header[0] not in {"clade_name", "clade"}:
        data = rows
    for row in data:
        if not row:
            continue
        token = sgb_token(row[0])
        if not token:
            continue
        values = [_to_float(cell) for cell in row[1:]]
        n_samples = len(values)
        present = sum(1 for value in values if value > 0.0)
        mean_ab = (sum(values) / n_samples) if n_samples else 0.0
        stats[token] = {"samples_present": present, "mean_abundance": mean_ab}
    return stats


def select_clades(
    sgb_stats: dict[str, dict[str, float | int]],
    explicit_clades: list[str] | None,
    top_n: int,
    min_samples: int,
) -> tuple[list[str], list[list[str]]]:
    """Return (selected clade IDs, summary rows). Never raises for empty cohorts."""
    summary: list[list[str]] = []
    selected: list[str] = []
    explicit = list(explicit_clades or [])

    if explicit:
        for clade in explicit:
            info = sgb_stats.get(clade, {})
            present = int(info.get("samples_present", 0))
            mean_ab = float(info.get("mean_abundance", 0.0))
            if not is_valid_clade(clade):
                status = "invalid"
            else:
                status = "selected"
                selected.append(clade)
            summary.append(
                [
                    clade,
                    str(present),
                    f"{mean_ab:.6g}",
                    "explicit",
                    status,
                ]
            )
        return selected, summary

    ranked = sorted(
        (
            (clade, int(info["samples_present"]), float(info["mean_abundance"]))
            for clade, info in sgb_stats.items()
            if is_valid_clade(clade)
        ),
        key=lambda item: (-item[1], -item[2], item[0]),
    )
    qualified = [item for item in ranked if item[1] >= min_samples]
    limit = max(int(top_n), 0)
    selected_set = {item[0] for item in qualified[:limit]}
    selected = [item[0] for item in qualified[:limit]]

    for clade, present, mean_ab in ranked:
        if clade in selected_set:
            status = "selected"
        elif present < min_samples:
            status = "below_min_samples"
        else:
            status = "not_selected"
        summary.append(
            [
                clade,
                str(present),
                f"{mean_ab:.6g}",
                "auto",
                status,
            ]
        )
    return selected, summary


def write_outputs(
    selected: list[str],
    summary: list[list[str]],
    summary_out: Path,
    clades_out: Path,
) -> None:
    summary_out.parent.mkdir(parents=True, exist_ok=True)
    clades_out.parent.mkdir(parents=True, exist_ok=True)
    with summary_out.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(SUMMARY_HEADER)
        writer.writerows(summary)
    clades_out.write_text("".join(f"{clade}\n" for clade in selected))


def main() -> int:
    args = parse_args()
    if args.top_n < 0:
        print("WARNING: --top-n is negative; treating as 0 (no auto-selected clades).", file=sys.stderr)
        args.top_n = 0
    if args.min_samples < 1:
        print("WARNING: --min-samples < 1; treating as 1.", file=sys.stderr)
        args.min_samples = 1

    sgb_stats = read_sgb_stats(Path(args.sgb_table))
    explicit = parse_clade_list(args.clades)
    selected, summary = select_clades(sgb_stats, explicit, args.top_n, args.min_samples)
    write_outputs(selected, summary, Path(args.summary_out), Path(args.clades_out))

    mode = "explicit" if explicit else "auto"
    print(
        f"StrainPhlAn clade selection ({mode}): {len(selected)} selected "
        f"(top_n={args.top_n}, min_samples={args.min_samples}). "
        "Research use only — not for clinical or diagnostic use.",
        file=sys.stderr,
    )
    if not selected:
        print(
            "No SGB clades qualified for StrainPhlAn tree building "
            "(small cohort, no shared SGBs, or no valid --strainphlan_clades). "
            "Continuing without clade trees.",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
