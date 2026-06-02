# DSA-phasing

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0.0-brightgreen.svg)](https://snakemake.github.io)

A Snakemake workflow for phasing long reads as well short reads using Donor-Specific genome Assemblies (DSAs). Reads are aligned to a per-sample DSA and assigned to haplotypes based on contig name tags, producing phased CRAM files with embedded epigenomic annotations (FIRE, nucleosome calls). The codes are largely based on the original [DSA-phasing workflow](https://github.com/vollgerlab/DSA-phasing)

## Features

- **Multi-sample / multi-BAM**: process many samples in one run via a manifest file
- **PacBio and ONT**: automatic per-sample detection of sequencing technology
- **Haplotype tagging**: reads tagged with `HP` (haplotype) and `oh` (original haplotype before MAPQ filtering)
- **Epigenomic annotations**: fibertools-rs FIRE calling; ONT modkit base modification calling + nucleosome detection
- **Optional shared-reference realignment**: realign DSA-phased reads to a common reference for cross-sample comparison

## Note for the Users (last updated: 05-29-2026)
- For the alignment of the PacBio HiFi long-read sequencing data (Standard PacBio or Fiber-seq) to their own genome of origin (i.e., DSA), we use mimimap2 with `--preset lr:hqae`
- For ONT alignment to the DSA, we are using `--preset lr:hq`, which is subjected to possible change in the near future after benchmarking.

## Quick start

### Prerequisites

Install [pixi](https://pixi.sh) (handles all Snakemake and tool dependencies):

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

### Installation

```bash
git clone https://github.com/vollgerlab/DSA-phasing.git
cd DSA-phasing
pixi install
```

### Run the test

```bash
pixi run test
```

### Run on your data

1. Create a config file (see [config/README.md](config/README.md) for all options):

```yaml
manifest: "path/to/manifest.tbl"
```

2. Create a manifest file (tab-separated):

```
sample	dsa	bam	h1_tag	h2_tag	platform	fiber-seq
sample1	/path/to/sample1.dsa.fa	/path/to/sample1.bam	hap1	hap2	pacbio	true
sample2	/path/to/sample2.dsa.fa	/path/to/s2_a.bam,/path/to/s2_b.bam	h1	h2	ont	true
sample3	/path/to/sample3.dsa.fa	/path/to/sample3.bam	hap1	hap2	illumina	false
```

| Column | Description |
|---|---|
| `sample` | Unique sample identifier |
| `dsa` | Path to the donor-specific assembly FASTA |
| `bam` | Input BAM/CRAM path(s), comma-separated for multiple files |
| `h1_tag` | Substring in DSA contig names identifying haplotype 1 (`NA` defaults to `h1`) |
| `h2_tag` | Substring in DSA contig names identifying haplotype 2 (`NA` defaults to `h2`) |
| `platform` | One of `pacbio`, `ont`, `illumina` — drives the aligner choice |
| `fiber-seq` | `true` or `false` — drives the modkit/FIRE branch. `illumina + fiber-seq=true` is rejected. |

3. Run:

```bash
pixi run snakemake --configfile config/config.yaml -j 64
```

Or with SLURM:

```bash
pixi run snakemake --configfile config/config.yaml --profile profiles/slurm-executor
```

## Workflow overview

```
BAM/CRAM(s) per sample
|
+-- long-read -> extract_fastq -> align (minimap2 -> DSA) ------+
+-- illumina  ------------------> align_illumina (bwa -> DSA) --+
                                                                |
                                                       haplotag_and_sort
                                                                |
                                                                +-- fiber-seq -> (modkit if ONT) -> fire -+
                                                                +-- non-fiber-seq ------------------------+
                                                                                                          |
                                                                                                     merge_sample
                                                                                                          |
                                                                                                          +-- fiber-seq -> qc_fiberseq (ft validate + ft qc)
                                                                                                          +-- long-read non-fs -> qc_longread_nofs (samtools stats)
                                                                                                          +-- illumina -> qc_illumina (stats + flagstat)
```

## Outputs

All outputs are written to `results/`:

| File | Description |
|---|---|
| `{sample}.dsa.cram` | Phased CRAM aligned to the DSA with HP tags and FIRE annotations |
| `{sample}.qc.tbl.gz` | fibertools QC table |
| `{sample}.shared.ref.cram` | (Optional) Realigned to the shared reference |

## Configuration

See [config/README.md](config/README.md) for the full list of configuration options.

## Disclaimer

The repository is entirely based on the [Mitchell's DSA-phasing repository](https://github.com/ryansohny/DSA-phasing), and re-written with the help of Claude code (Opus-4.6,-4.7 and -4.8)

## Citation

If you use this workflow, please cite this repository:

> Vollger, M.R. DSA-phasing. https://github.com/vollgerlab/DSA-phasing
or
> Sohn, M-H DSA-phasing. https://github.com/ryansohny/DSA-phasing

## Stage 2: per-tissue merge + coverage + mCG

After the phasing pipeline produces per-sample DSA CRAMs in `results/`, a second
workflow groups them by tissue, merges (single-member groups are copied), runs
mosdepth coverage (PacBio + ONT) and pb-CpG-tools mCG extraction (PacBio only),
and writes to a required `output_dir`. Grouping comes from the `decode.tsv`
produced by `make-manifest.py`; one donor + one platform per run.

```bash
# dry run
pixi run snakemake -n --snakefile workflow/merge-tissue.smk \
  --configfile config/config_merge.example.yaml

# submit (one SLURM job per tissue)
pixi run snakemake --snakefile workflow/merge-tissue.smk \
  --configfile config/config_merge.example.yaml --profile profiles/slurm-executor

# reclaim space after inspecting outputs (deletes consumed per-sample CRAMs)
pixi run snakemake --snakefile workflow/merge-tissue.smk \
  --configfile config/config_merge.example.yaml \
  --config delete_originals=true --profile profiles/slurm-executor
```

Re-runs are idempotent: once a tissue's merged CRAM exists in `output_dir`, that
group is skipped even if the per-sample inputs were deleted.
