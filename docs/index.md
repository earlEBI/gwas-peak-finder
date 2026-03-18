---
layout: default
title: GWAS Peak Finder Batch
---

# GWAS Peak Finder Batch

`peak_finder-2.2.py` finds lead SNP associations by trait/chromosome using a significance threshold and a configurable genomic window.

Companion tool to the legacy `gwasAssociationFilter` peak-finder in `EBISPOT/gwas-utils`.

Automatic tie resolution is enabled by default in this workflow.

## Inputs

Provide a tab-separated file containing:

- `trait`
- `rs_id`
- `pvalue`
- `chromosome`
- `bp_location`

## Output

The script adds `isTopAssociation` with:

- `true` for selected lead variants
- `false` for non-leads or sub-threshold variants
- `AUTO-REVIEWED-FALSE` for tied candidates; one in-window representative is automatically set to `true`

## Run

```bash
python peak_finder-2.2.py -f input.tsv -o output.tsv
```

Optional:

```bash
python peak_finder-2.2.py -f input.tsv -o output.tsv -w 100000 -t 1e-5 -p
```
