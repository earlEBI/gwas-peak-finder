---
layout: default
title: Peak Finder 2.2
---

# Peak Finder 2.2

`peak_finder-2.2.py` finds lead SNP associations by trait/chromosome using a significance threshold and a configurable genomic window.

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
- `AUTO-REVIEWED-FALSE` for tied candidates needing manual review

## Run

```bash
python peak_finder-2.2.py -f input.tsv -o output.tsv
```

Optional:

```bash
python peak_finder-2.2.py -f input.tsv -o output.tsv -w 100000 -t 1e-5 -p
```
