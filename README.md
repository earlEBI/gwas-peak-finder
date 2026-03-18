# Peak Finder 2.2

This repository packages `peak_finder-2.2.py` for release in two versions:

- `v1.0.0`: original script as provided
- `v1.1.0`: reviewed patch with compatibility and logic fixes

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
