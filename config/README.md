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
- **Description**: Minimum mapping quality threshold for filtering reads during haplotagging and sorting. Reads with mapping quality below this value will be excluded.
- **Example**: `min_mapq: 20`

### `ft_nuc_params`

- **Type**: String
- **Default**: "" (empty)
- **Description**: Additional parameters for the `ft add-nucleosomes` command in the modkit rule. Used for nucleosome detection and modification analysis.
- **Example**: `ft_nuc_params: "--nucleosome-length 60"`

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
