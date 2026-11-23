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
| `--metaphlan_db` | `null` | MetaPhlAn 4 DB directory or `.tar`/`.tar.gz`. Default on `refdb.treetopunder.com` will be wired once uploaded. |
| `--host_index` | `null` | Host Bowtie2 index directory or `.tar`/`.tar.gz`. Required unless `--skip_host_removal`. Default on `refdb.treetopunder.com` will be wired once uploaded. |
| `--skip_host_removal` | `false` | Skip Bowtie2 host removal |
| `--analysis_type` | `rel_ab_w_read_stats` | MetaPhlAn `-t` (`rel_ab_w_read_stats` or `rel_ab`) |
| `--add_viruses` | `false` | Pass `--add_viruses` |
| `--unclassified_estimation` | `true` | Pass `--unclassified_estimation` |
| `--min_read_len` | `70` | fastp `--length_required` |
| `--skip_krona` | `false` | Skip Krona plots |
| `--skip_phyloseq` | `false` | Skip phyloseq RDS |

Host indexes and MetaPhlAn databases may be an unpacked directory or a `.tar` / `.tar.gz` archive (unpacked with `nf-core/untar`).

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

| Resource | Target | Notes |
|---|---|---|
| MetaPhlAn 4 | `mpa_vJun23_CHOCOPhlAnSGB_202403` | CHOCOPhlAnSGB index used by MetaPhlAn 4.1.1 |
| Host genome | `GRCh38_noalt_as` | Bowtie2 index for human host removal |

Place unpacked directories (or tar archives) on local disk or S3 and pass them as `--metaphlan_db` / `--host_index`. TreeTopFlow hosted copies at `refdb.treetopunder.com` will be wired as defaults once uploaded.

## Run

```bash
nextflow run . \
  --input samplesheet.csv \
  --outdir results \
  --metaphlan_db /refs/mpa_vJun23_CHOCOPhlAnSGB_202403 \
  --host_index /refs/GRCh38_noalt_as \
  --metadata metadata.tsv \
  -profile docker
```

AWS Batch with S3 work and output directories:

```bash
nextflow run . \
  -profile awsbatch \
  --aws_queue <queue> \
  --aws_region us-east-1 \
  -work-dir s3://bucket/work \
  --input samplesheet.csv \
  --outdir s3://bucket/results \
  --metaphlan_db s3://bucket/refs/mpa_vJun23_CHOCOPhlAnSGB_202403 \
  --host_index s3://bucket/refs/GRCh38_noalt_as
```

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
