#!/usr/bin/env python3
"""Unit tests for ttf-metaphlan helper scripts and JSON schemas."""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
FIXTURES = ROOT / "tests" / "fixtures"


class TidyTablesTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmpdir = Path(tempfile.mkdtemp(prefix="ttf-metaphlan-"))
        self.addCleanup(shutil.rmtree, self.tmpdir, True)

    def test_species_excludes_sgb_and_higher_ranks(self) -> None:
        species_out = self.tmpdir / "merged_species_table.tsv"
        sgb_out = self.tmpdir / "merged_sgb_table.tsv"
        subprocess.check_call(
            [
                sys.executable,
                str(BIN / "tidy_metaphlan_tables.py"),
                "--merged",
                str(FIXTURES / "metaphlan" / "merged_abundance_table.tsv"),
                "--species-out",
                str(species_out),
                "--sgb-out",
                str(sgb_out),
            ]
        )
        species_lines = species_out.read_text().splitlines()
        sgb_lines = sgb_out.read_text().splitlines()
        self.assertTrue(species_lines[0].startswith("clade_name"))
        species_clades = [line.split("\t")[0] for line in species_lines[1:]]
        sgb_clades = [line.split("\t")[0] for line in sgb_lines[1:]]
        self.assertTrue(all("|s__" in clade and "|t__" not in clade for clade in species_clades))
        self.assertTrue(all(clade.split("|")[-1].startswith("s__") for clade in species_clades))
        self.assertTrue(all(clade.split("|")[-1].startswith("t__") for clade in sgb_clades))
        self.assertIn(
            "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_fragilis",
            species_clades,
        )
        self.assertNotIn(
            "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_fragilis|t__SGB1855",
            species_clades,
        )
        self.assertEqual(len(species_clades), 3)
        self.assertEqual(len(sgb_clades), 3)
        fragilis = [line for line in species_lines if "s__Bacteroides_fragilis" in line][0].split("\t")
        self.assertEqual(fragilis[1], "50.0")
        self.assertEqual(fragilis[2], "30.0")

    def test_host_removal_summary_from_bowtie2_logs(self) -> None:
        host_out = self.tmpdir / "host_removal_summary.tsv"
        subprocess.check_call(
            [
                sys.executable,
                str(BIN / "tidy_metaphlan_tables.py"),
                "--merged",
                str(FIXTURES / "metaphlan" / "merged_abundance_table.tsv"),
                "--species-out",
                str(self.tmpdir / "species.tsv"),
                "--sgb-out",
                str(self.tmpdir / "sgb.tsv"),
                "--bowtie2-logs",
                str(FIXTURES / "bowtie2_logs" / "sample1.bowtie2.log"),
                str(FIXTURES / "bowtie2_logs" / "sample2.bowtie2.log"),
                "--host-summary-out",
                str(host_out),
            ]
        )
        rows = [line.split("\t") for line in host_out.read_text().splitlines()]
        self.assertEqual(
            rows[0],
            ["sample", "total_reads", "aligned_reads", "unaligned_reads", "overall_alignment_rate"],
        )
        by_sample = {row[0]: row for row in rows[1:]}
        self.assertEqual(by_sample["sample1"][1], "20000")
        self.assertEqual(by_sample["sample1"][2], "3000")
        self.assertEqual(by_sample["sample1"][3], "17000")
        self.assertEqual(by_sample["sample1"][4], "15.00")
        self.assertEqual(by_sample["sample2"][1], "18000")
        self.assertEqual(by_sample["sample2"][2], "1800")
        self.assertEqual(by_sample["sample2"][4], "10.00")


class KronaConverterTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmpdir = Path(tempfile.mkdtemp(prefix="ttf-krona-"))
        self.addCleanup(shutil.rmtree, self.tmpdir, True)

    def test_uses_sgb_leaves_and_unclassified(self) -> None:
        out = self.tmpdir / "sample1.krona.txt"
        subprocess.check_call(
            [
                sys.executable,
                str(BIN / "metaphlan_to_krona.py"),
                "--profile",
                str(FIXTURES / "profiles" / "sample1_profile.txt"),
                "--output",
                str(out),
            ]
        )
        rows = [line.split("\t") for line in out.read_text().splitlines() if line.strip()]
        self.assertEqual(rows[0][1], "UNCLASSIFIED")
        self.assertEqual(rows[0][0], "100")
        taxa = [row[-1] for row in rows[1:]]
        self.assertIn("SGB1855", taxa)
        self.assertIn("SGB15316", taxa)
        # Species-only rows must not be emitted when SGB leaves exist (avoids double counting).
        self.assertFalse(any(row[-1] == "Bacteroides_fragilis" for row in rows))
        fragilis = [row for row in rows if row[-1] == "SGB1855"][0]
        self.assertEqual(fragilis[0], "500")
        self.assertEqual(fragilis[1], "Bacteria")
        self.assertEqual(fragilis[7], "Bacteroides_fragilis")


class PhyloseqScriptTests(unittest.TestCase):
    def test_build_phyloseq_if_rscript_available(self) -> None:
        rscript = shutil.which("Rscript")
        if not rscript:
            self.skipTest("Rscript is not available in this environment")
        tmpdir = Path(tempfile.mkdtemp(prefix="ttf-phyloseq-"))
        self.addCleanup(shutil.rmtree, tmpdir, True)
        species = tmpdir / "species.tsv"
        sgb = tmpdir / "sgb.tsv"
        subprocess.check_call(
            [
                sys.executable,
                str(BIN / "tidy_metaphlan_tables.py"),
                "--merged",
                str(FIXTURES / "metaphlan" / "merged_abundance_table.tsv"),
                "--species-out",
                str(species),
                "--sgb-out",
                str(sgb),
            ]
        )
        rds = tmpdir / "metaphlan_phyloseq.rds"
        try:
            subprocess.check_call(
                [
                    rscript,
                    str(BIN / "build_phyloseq.R"),
                    "--species",
                    str(species),
                    "--sgb",
                    str(sgb),
                    "--metadata",
                    str(FIXTURES / "metadata.tsv"),
                    "--db",
                    str(FIXTURES / "metaphlan_db"),
                    "--output",
                    str(rds),
                ]
            )
        except subprocess.CalledProcessError as exc:
            self.skipTest(f"phyloseq is not installed: {exc}")
        self.assertTrue(rds.exists())
        self.assertGreater(rds.stat().st_size, 0)


class StrainphlanCladeSelectionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmpdir = Path(tempfile.mkdtemp(prefix="ttf-strainphlan-"))
        self.addCleanup(shutil.rmtree, self.tmpdir, True)
        sys.path.insert(0, str(BIN))
        import select_strainphlan_clades as selector  # noqa: WPS433

        self.selector = selector
        self.sgb = FIXTURES / "strainphlan" / "sgb_table.tsv"

    def _run(
        self,
        clades: str = "",
        top_n: int = 5,
        min_samples: int = 4,
        table: Path | None = None,
    ) -> tuple[str, str]:
        summary = self.tmpdir / "strainphlan_summary.tsv"
        clade_list = self.tmpdir / "selected_clades.txt"
        subprocess.check_call(
            [
                sys.executable,
                str(BIN / "select_strainphlan_clades.py"),
                "--sgb-table",
                str(table or self.sgb),
                "--clades",
                clades,
                "--top-n",
                str(top_n),
                "--min-samples",
                str(min_samples),
                "--summary-out",
                str(summary),
                "--clades-out",
                str(clade_list),
            ]
        )
        return summary.read_text(), clade_list.read_text()

    def test_explicit_clades_are_validated_and_selected(self) -> None:
        summary, clade_list = self._run(clades="t__SGB1877, not-a-clade, t__SGB15316")
        selected = [line for line in clade_list.splitlines() if line]
        self.assertEqual(selected, ["t__SGB1877", "t__SGB15316"])
        rows = [line.split("\t") for line in summary.splitlines()]
        self.assertEqual(rows[0], self.selector.SUMMARY_HEADER)
        by_clade = {row[0]: row for row in rows[1:]}
        self.assertEqual(by_clade["t__SGB1877"][3], "explicit")
        self.assertEqual(by_clade["t__SGB1877"][4], "selected")
        self.assertEqual(by_clade["t__SGB1877"][1], "4")
        self.assertEqual(by_clade["not-a-clade"][4], "invalid")
        self.assertEqual(by_clade["t__SGB15316"][4], "selected")
        self.assertNotIn("not-a-clade", selected)

    def test_auto_top_n_ranks_by_prevalence_then_abundance(self) -> None:
        summary, clade_list = self._run(top_n=2, min_samples=4)
        selected = [line for line in clade_list.splitlines() if line]
        self.assertEqual(selected, ["t__SGB1877", "t__SGB10068"])
        rows = [line.split("\t") for line in summary.splitlines()[1:]]
        by_clade = {row[0]: row for row in rows}
        self.assertEqual(by_clade["t__SGB1877"][4], "selected")
        self.assertEqual(by_clade["t__SGB10068"][4], "selected")
        self.assertEqual(by_clade["t__SGB15316"][4], "below_min_samples")
        self.assertEqual(by_clade["t__SGB17256"][4], "below_min_samples")
        self.assertEqual(by_clade["t__SGB1877"][3], "auto")
        extra = self._run(top_n=1, min_samples=3)
        self.assertEqual(extra[1].splitlines(), ["t__SGB1877"])

    def test_auto_zero_clades_when_below_min_samples(self) -> None:
        summary, clade_list = self._run(top_n=5, min_samples=8)
        self.assertEqual(clade_list, "")
        rows = [line.split("\t") for line in summary.splitlines()]
        self.assertEqual(rows[0], self.selector.SUMMARY_HEADER)
        self.assertTrue(all(row[4] == "below_min_samples" for row in rows[1:]))
        empty_table = self.tmpdir / "empty_sgb.tsv"
        empty_table.write_text("clade_name\tsample1\n")
        summary2, clade_list2 = self._run(min_samples=2, table=empty_table)
        self.assertEqual(clade_list2, "")
        self.assertEqual(summary2.splitlines()[0], "\t".join(self.selector.SUMMARY_HEADER))
        self.assertEqual(len(summary2.splitlines()), 1)

    def test_select_clades_helper_matches_cli(self) -> None:
        stats = self.selector.read_sgb_stats(self.sgb)
        selected, rows = self.selector.select_clades(stats, [], top_n=5, min_samples=3)
        self.assertEqual(selected, ["t__SGB1877", "t__SGB10068", "t__SGB15316"])
        statuses = {row[0]: row[4] for row in rows}
        self.assertEqual(statuses["t__SGB17256"], "below_min_samples")


class SchemaTests(unittest.TestCase):
    def test_schemas_are_valid_json_schema(self) -> None:
        try:
            from jsonschema.validators import validator_for
        except ImportError as exc:  # pragma: no cover
            self.skipTest(f"jsonschema is not installed: {exc}")

        for rel in ("nextflow_schema.json", "assets/schema_input.json"):
            path = ROOT / rel
            with path.open() as handle:
                schema = json.load(handle)
            validator_for(schema).check_schema(schema)
            self.assertIn("$schema", schema)
            self.assertEqual(schema["type"] in {"object", "array"}, True)

    def test_nextflow_schema_exposes_required_params(self) -> None:
        schema = json.loads((ROOT / "nextflow_schema.json").read_text())
        defs = schema.get("$defs") or schema.get("definitions")
        self.assertIsNotNone(defs)
        props = {}
        for group in defs.values():
            props.update(group.get("properties", {}))
        for key in [
            "input",
            "outdir",
            "metadata",
            "metaphlan_db",
            "host_index",
            "skip_host_removal",
            "analysis_type",
            "add_viruses",
            "unclassified_estimation",
            "min_read_len",
            "skip_krona",
            "skip_phyloseq",
            "run_strainphlan",
            "strainphlan_clades",
            "strainphlan_top_n",
            "strainphlan_min_samples",
            "strainphlan_phylophlan_mode",
        ]:
            self.assertIn(key, props)
        self.assertEqual(props["analysis_type"]["enum"], ["rel_ab_w_read_stats", "rel_ab"])
        self.assertFalse(props["run_strainphlan"]["default"])
        self.assertIsNone(props["strainphlan_clades"]["default"])
        self.assertEqual(props["strainphlan_top_n"]["default"], 5)
        self.assertEqual(props["strainphlan_min_samples"]["default"], 4)
        self.assertEqual(props["strainphlan_phylophlan_mode"]["enum"], ["fast", "accurate"])
        self.assertEqual(props["strainphlan_phylophlan_mode"]["default"], "fast")
        self.assertEqual(props["min_read_len"]["default"], 70)
        self.assertFalse(props["skip_host_removal"]["default"])
        self.assertTrue(props["unclassified_estimation"]["default"])
        self.assertEqual(
            props["metaphlan_db"]["default"],
            "s3://codonyat/refdb/public/metaphlan4/mpa_vJun23_CHOCOPhlAnSGB_202403",
        )
        self.assertEqual(
            props["host_index"]["default"],
            "s3://codonyat/refdb/public/host/GRCh38_noalt_as",
        )
        self.assertIn("refdb.treetopunder.com", props["metaphlan_db"]["description"])
        self.assertIn("refdb.treetopunder.com", props["host_index"]["description"])
        self.assertNotIn("pending", json.dumps(props["metaphlan_db"]).lower())
        self.assertNotIn("pending", json.dumps(props["host_index"]).lower())

    def test_nextflow_config_s3_defaults_and_offline_test_profile(self) -> None:
        config = (ROOT / "nextflow.config").read_text()
        self.assertIn(
            "s3://codonyat/refdb/public/metaphlan4/mpa_vJun23_CHOCOPhlAnSGB_202403",
            config,
        )
        self.assertIn("s3://codonyat/refdb/public/host/GRCh38_noalt_as", config)
        self.assertNotIn("will be wired once uploaded", config)
        self.assertIn("run_strainphlan", config)
        self.assertIn("strainphlan_phylophlan_mode", config)
        test_config = (ROOT / "conf" / "test.config").read_text()
        self.assertIn("tests/fixtures/metaphlan_db", test_config)
        self.assertIn("tests/fixtures/host_index", test_config)
        self.assertNotIn("s3://codonyat", test_config)


if __name__ == "__main__":
    os.chdir(ROOT)
    unittest.main(verbosity=2)
