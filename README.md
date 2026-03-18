# GWAS Peak Finder Batch

This repository contains a companion peak-finder based on EBISPOT `gwasAssociationFilter`, optimized for multi-trait batch processing and optional automatic tie resolution.

## Release lines

- `v1.0.0`: original script as provided
- `v1.1.0`: reviewed patch with compatibility and logic fixes
- `v2.0.0`: companion-tool naming update to distinguish this workflow from the legacy upstream tool

## Upstream relationship

- Upstream tool: `EBISPOT/gwas-utils` → `gwasAssociationFilter`
- This repository is intentionally maintained as a separate companion workflow and does not replace upstream defaults.

## v1.1.0 patch highlights

- Fixes pandas compatibility by avoiding float writes into string-typed columns.
- Fixes `--prune` so significance filtering runs correctly.
- Restricts equal p-value tie handling to variants within the current window.
- Handles invalid/blank p-values safely and treats `p=0` as most significant.
- Adds threshold argument validation (`0 < threshold < 1`).

## What the script does

`peak_finder-2.2.py` processes GWAS association tables and marks lead variants per trait/chromosome in a configurable genomic window.

Required columns:

- `trait`
- `rs_id`
- `pvalue`
- `chromosome`
- `bp_location`

It outputs the same table with an added `isTopAssociation` column.

## Quick start

```bash
python peak_finder-2.2.py -f input.tsv -o output.tsv
```

Optional arguments:

- `-w/--window` (default: `100000`)
- `-t/--threshold` (default: `1e-5`)
- `-p/--prune` (drop sub-threshold rows before peak calling)

## GitHub Pages

The project overview page is in [`docs/index.md`](docs/index.md).
