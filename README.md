# DSA-phasing

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0.0-brightgreen.svg)](https://snakemake.github.io)

A Snakemake workflow for phasing long reads using Donor-Specific Assemblies (DSAs). Reads are aligned to a per-sample DSA and assigned to haplotypes based on contig name tags, producing phased CRAM files with embedded epigenomic annotations (FIRE, nucleosome calls).

## Features

- **Multi-sample / multi-BAM**: process many samples in one run via a manifest file
- **PacBio and ONT**: automatic per-sample detection of sequencing technology
- **Haplotype tagging**: reads tagged with `HP` (haplotype) and `oh` (original haplotype before MAPQ filtering)
- **Epigenomic annotations**: fibertools-rs FIRE calling; ONT modkit base modification calling + nucleosome detection
- **Optional shared-reference realignment**: realign DSA-phased reads to a common reference for cross-sample comparison

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
sample	dsa	bam	h1_tag	h2_tag	technology
sample1	/path/to/sample1.dsa.fa	/path/to/sample1.bam	hap1	hap2	pacbio
sample2	/path/to/sample2.dsa.fa	/path/to/s2_a.bam,/path/to/s2_b.bam	h1	h2	ont
```

| Column | Description |
|---|---|
| `sample` | Unique sample identifier |
| `dsa` | Path to the donor-specific assembly FASTA |
| `bam` | Input BAM path(s), comma-separated for multiple files |
| `h1_tag` | Substring in DSA contig names identifying haplotype 1 (`NA` defaults to `h1`) |
| `h2_tag` | Substring in DSA contig names identifying haplotype 2 (`NA` defaults to `h2`) |
| `technology` | `pacbio` or `ont` (optional; can also set globally in config) |

3. Run:

```bash
pixi run snakemake --configfile config/config.yaml -j 64
```

Or with SLURM:

```bash
pixi run snakemake --configfile config/config.yaml --profile workflow/profiles/default
```

## Workflow overview

```
BAM(s) per sample
  |
  +- extract_fastq -> align (minimap2 -> DSA) -> haplotag_and_sort (HP/oh tags)
  |                                                     |
  |                                    +----------------+----------------+
  |                                    | ONT                     PacBio  |
  |                                modkit -> fire                  fire  |
  |                                    +----------------+----------------+
  |                                                     |
  +-----------------------------------> merge_sample -> qc
                                                |
                                     (optional) +-> realign_to_shared_ref
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

## Citation

If you use this workflow, please cite this repository:

> Vollger, M.R. DSA-phasing. https://github.com/vollgerlab/DSA-phasing
