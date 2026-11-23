# ttf-metaphlan

Shotgun Taxonomic Profiling (MetaPhlAn 4) — a lean Nextflow DSL2 pipeline for TreeTopFlow.

**Research use only. Not for clinical or diagnostic use.**

Version `0.1.0` profiles Illumina shotgun metagenomes with MetaPhlAn 4 only: fastp QC/trimming, optional Bowtie2 host removal, MetaPhlAn 4 taxonomic profiling, merged/tidy tables, optional phyloseq RDS and Krona plots, and MultiQC.

## Inputs

`--input` samplesheet CSV (required):

```csv
sample,fastq_1,fastq_2
SAMPLE_PE,/path/to/SAMPLE_PE_R1.fastq.gz,/path/to/SAMPLE_PE_R2.fastq.gz
SAMPLE_SE,/path/to/SAMPLE_SE.fastq.gz,
```

- Paired-end: set `fastq_1` and `fastq_2`
- Single-end: leave `fastq_2` blank
- Paths may be absolute or relative to the launch directory / pipeline root

Optional `--metadata` TSV for phyloseq `sample_data`. The first column must be `ID` or `sample` and must match samplesheet sample IDs:

```tsv
sample	condition	batch
SAMPLE_PE	control	1
SAMPLE_SE	treatment	1
```

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | required | Samplesheet CSV |
| `--outdir` | required | Output directory (local or `s3://`) |
| `--metadata` | `null` | Optional sample metadata TSV |
| `--metaphlan_db` | `s3://codonyat/refdb/public/metaphlan4/mpa_vJun23_CHOCOPhlAnSGB_202403` | Unpacked MetaPhlAn 4 CHOCOPhlAnSGB directory (or `.tar`/`.tar.gz`). HTTPS: `https://refdb.treetopunder.com/metaphlan4/mpa_vJun23_CHOCOPhlAnSGB_202403/` |
| `--host_index` | `s3://codonyat/refdb/public/host/GRCh38_noalt_as` | Unpacked GRCh38_noalt_as Bowtie2 directory (or `.tar`/`.tar.gz`). Required unless `--skip_host_removal`. HTTPS: `https://refdb.treetopunder.com/host/GRCh38_noalt_as/` |
| `--skip_host_removal` | `false` | Skip Bowtie2 host removal |
| `--analysis_type` | `rel_ab_w_read_stats` | MetaPhlAn `-t` (`rel_ab_w_read_stats` or `rel_ab`) |
| `--add_viruses` | `false` | Pass `--add_viruses` |
| `--unclassified_estimation` | `true` | Pass `--unclassified_estimation` |
| `--min_read_len` | `70` | fastp `--length_required` |
| `--skip_krona` | `false` | Skip Krona plots |
| `--skip_phyloseq` | `false` | Skip phyloseq RDS |

Host indexes and MetaPhlAn databases may be an unpacked directory (local or `s3://`, trailing slash optional) or a `.tar` / `.tar.gz` archive (unpacked with `nf-core/untar`). Unpacked directories skip `UNTAR` and are staged as directory channels so MetaPhlAn, Bowtie2, and phyloseq see the index files (including `*.nwk` for `phy_tree`).

## Outputs

```
outdir/
  qc/fastp/*.fastp.json
  qc/fastp/*.fastp.html
  host_removal/*.bowtie2.log
  host_removal/host_removal_summary.tsv
  metaphlan/*.profile.txt
  metaphlan/*.biom
  metaphlan/merged_abundance_table.tsv
  metaphlan/merged_species_table.tsv
  metaphlan/merged_sgb_table.tsv
  phyloseq/metaphlan_phyloseq.rds
  krona/krona_combined.html
  krona/<sample>.krona.html
  multiqc/multiqc_report.html
  multiqc/multiqc_data/multiqc_general_stats.txt
  pipeline_info/software_versions.yml
```

Host-aligned BAMs and host reads are **never** published. Only unmapped reads continue to MetaPhlAn.

`merged_species_table.tsv` is the `s__` level (excluding `t__`). `merged_sgb_table.tsv` is the SGB `t__` level. Phyloseq uses SGB abundances when present, with taxonomy ranks Kingdom through SGB, optional `sample_data`, and optional `phy_tree` if an SGB `.nwk` file is inside `--metaphlan_db`.

## Reference databases

Defaults are unpacked prefixes on S3 (`s3://codonyat/refdb/public/...`), MD5-verified and mirrored at `https://refdb.treetopunder.com`. Override with `--metaphlan_db` / `--host_index` for a local copy or a `.tar`/`.tar.gz` archive. `-profile test` points both params at `tests/fixtures/` so stub-runs stay offline.

| Resource | Default S3 prefix | HTTPS mirror | Contents |
|---|---|---|---|
| MetaPhlAn 4 (`mpa_vJun23_CHOCOPhlAnSGB_202403`) | `s3://codonyat/refdb/public/metaphlan4/mpa_vJun23_CHOCOPhlAnSGB_202403` | https://refdb.treetopunder.com/metaphlan4/mpa_vJun23_CHOCOPhlAnSGB_202403/ | Flat unpacked dir: `mpa_vJun23_CHOCOPhlAnSGB_202403.pkl`, `.{1,2,3,4,rev.1,rev.2}.bt2l`, `mpa_vJun23_CHOCOPhlAnSGB_202403.nwk`, `mpa_latest`, etc. |
| Host genome (`GRCh38_noalt_as`) | `s3://codonyat/refdb/public/host/GRCh38_noalt_as` | https://refdb.treetopunder.com/host/GRCh38_noalt_as/ | Unpacked dir: `GRCh38_noalt_as.{1,2,3,4}.bt2`, `GRCh38_noalt_as.rev.{1,2}.bt2` |

## Run

```bash
nextflow run . \
  --input samplesheet.csv \
  --outdir results \
  --metadata metadata.tsv \
  -profile docker
```

The S3 defaults above are used unless you pass `--metaphlan_db` / `--host_index`.

AWS Batch (TreeTopFlow controllers pass the worker config with `-c`; do not require `--aws_queue`):

```bash
nextflow run . \
  -profile awsbatch \
  -c /path/to/worker/awsbatch.config \
  -work-dir s3://bucket/work \
  --input samplesheet.csv \
  --outdir s3://bucket/results
```

Optional CLI overrides (`--aws_queue`, `--aws_region`, `--aws_cli_path`) are applied only when set, so they do not clobber the worker `-c` file.

Stub-run (DAG compile, no containers required):

```bash
nextflow run . -profile test -stub-run --outdir results_host
nextflow run . -profile test -stub-run --skip_host_removal true --outdir results_nohost
```

Helper-script tests:

```bash
python tests/test_scripts.py
```

## Modules

Standard `nf-core/modules` (see `modules.json`):

- `nf-core/fastp`
- `nf-core/bowtie2/align` (`save_unaligned = true`)
- `nf-core/untar`
- `nf-core/metaphlan/metaphlan` (container `metaphlan=4.1.1`)
- `nf-core/metaphlan/mergemetaphlantables`
- `nf-core/krona/ktimporttext`
- `nf-core/multiqc`

Local modules: `METAPHLAN_TIDY_TABLES`, `METAPHLAN_TO_KRONA`, `BUILD_PHYLOSEQ`.

## License

MIT. nf-core modules retain their upstream MIT license. MetaPhlAn, Bowtie2, fastp, Krona, MultiQC, and phyloseq retain their respective licenses.
