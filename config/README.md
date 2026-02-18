# DSA-phasing Configuration

This document describes all configuration options available for the DSA-phasing workflow. Configuration is primarily done through `config.yaml` and the workflow validates all options against `workflow/config.schema.yaml`.

## Required Configuration

### `manifest`

- **Type**: String (file path)
- **Description**: Path to a table containing the manifest of samples to be processed. This file should contain sample information including paths to input files.
- **Example**: `manifest: "test/manifest.tbl"`

## Optional Configuration

### `max_threads`

- **Type**: Integer
- **Default**: 64
- **Description**: Maximum number of threads for parallel processing operations throughout the workflow.
- **Example**: `max_threads: 32`

### `ont`

- **Type**: Boolean
- **Default**: false
- **Description**: Enable Oxford Nanopore Technologies (ONT) specific processing. When true, affects output requirements and final CRAM file selection to optimize for ONT data characteristics.
- **Example**: `ont: true`

### `set-sm`

- **Type**: Boolean
- **Default**: false
- **Description**: Whether to set sample names in BAM headers to match the manifest. When enabled, BAM headers will be modified to ensure consistency with the provided manifest.
- **Example**: `set-sm: true`

### `mm2_preset`

- **Type**: String
- **Default**: "lr:hq"
- **Description**: Minimap2 preset parameter for alignment. Common options include 'lr:hq' for high-quality long reads, 'map-ont' for ONT reads, or 'map-pb' for PacBio reads.
- **Example**: `mm2_preset: "map-ont"`

### `mm2_extra_options`

- **Type**: String
- **Default**: "" (empty)
- **Description**: Additional command-line options to pass to minimap2 during alignment. Allows fine-tuning of alignment parameters beyond the preset.
- **Example**: `mm2_extra_options: "-k19 -w10"`

### `min_mapq`

- **Type**: Integer
- **Default**: 1
- **Description**: Minimum mapping quality threshold for haplotype assignment. Reads below this threshold still appear in the output but have their `HP` tag cleared (set to unphased). The original assignment is preserved in the `oh` tag for debugging.
- **Example**: `min_mapq: 20`

### `reset_mapq`

- **Type**: Integer
- **Default**: disabled
- **Description**: Reset MAPQ of mapped reads to this value during haplotagging. The original MAPQ is preserved in the `om` tag. Useful because DSA-alignment MAPQ values may not be meaningful for downstream tools that filter on MAPQ. By default, resets after haplotype assignment so `min_mapq` filtering uses the original value (see `reset_mapq_before`).
- **Example**: `reset_mapq: 60`

### `reset_mapq_before`

- **Type**: Boolean
- **Default**: false
- **Description**: When true (and `reset_mapq` is set), reset MAPQ before haplotype assignment. This means `min_mapq` filtering will use the new MAPQ value instead of the original. When false (default), MAPQ is reset after assignment so filtering uses the original alignment MAPQ.
- **Example**: `reset_mapq_before: true`

### `ft_nuc_params`

- **Type**: String
- **Default**: "" (empty)
- **Description**: Additional parameters for the `ft add-nucleosomes` command in the modkit rule. Used for nucleosome detection and modification analysis.
- **Example**: `ft_nuc_params: "--nucleosome-length 60"`

### `shared_reference`

- **Type**: String (file path) or List of file paths
- **Default**: [] (empty, disabled)
- **Description**: Path to a shared reference genome. When provided, the workflow will realign the final merged CRAMs to this reference in addition to the DSA-aligned outputs. Useful for comparing samples aligned to different DSAs on a common coordinate system.
- **Example**: `shared_reference: "/path/to/hg38.fa"`

### `keep_read_assignments`

- **Type**: Boolean
- **Default**: false
- **Description**: When true, read-to-haplotype assignment TSV files are saved to `results/{sm}/` instead of being placed in `temp/` and cleaned up. Useful for downstream analysis of phasing results.
- **Example**: `keep_read_assignments: true`

## Configuration Validation

The workflow automatically validates all configuration options against the schema defined in `workflow/config.schema.yaml`. Invalid configurations will cause the workflow to fail with descriptive error messages.

## Example Configuration

```yaml
# Required
manifest: "test/manifest.tbl"

# Optional (showing non-default values)
max_threads: 32
ont: true
set-sm: true
mm2_preset: "map-ont"
mm2_extra_options: "-k19"
min_mapq: 20
ft_nuc_params: "--nucleosome-length 60"
```

## Notes

- Only `manifest` is required; all other options have sensible defaults
- Boolean values should be lowercase: `true` or `false`
- String values should be quoted if they contain special characters
- The workflow will print the loaded configuration to stderr for verification
