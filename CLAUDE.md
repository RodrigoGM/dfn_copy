# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Status

This repository is currently empty — no source, build files, or docs exist yet.
This CLAUDE.md is the **design/spec contract** for the first implementation,
derived from the tool spec that was provided when this repo was initialized.
It intentionally does not describe a build/lint/test workflow, because none
exists yet — update this file once real commands exist, don't let it drift
into aspirational claims about a working build.

## What the tool is

A single-purpose C++17 CLI tool that turns a **DEFND-seq** gDNA BAM into a
bin × cell copy-number count matrix.

DEFND-seq (Olsen et al.) co-sequences RNA and DNA from single nuclei using the
10x Genomics Chromium Single Cell Multiome ATAC + Gene Expression kit
"as-is," producing a gDNA library alongside the cDNA library. The gDNA BAM is
therefore a standard 10x-barcoded, coordinate-sorted BAM — the per-nucleus
cell barcode lives in a BAM tag (`CB` by default, the Cell Ranger corrected
barcode, e.g. `CB:Z:AACGTGATCCGTAAGT-1`), not in the read sequence.

The tool does one thing: stream the sorted BAM once, assign each passing read
to a user-supplied genomic bin, and emit two gzipped bin × cell matrices (raw
counts and GC-corrected counts). It does not compute bins itself, and it is
not ATAC-fragment-aware — this is plain per-read binning over gDNA alignments.

## CLI contract

This is the interface any implementation must match; do not add flags beyond
this list or change these defaults without updating this section.

| Flag | Default | Meaning |
|---|---|---|
| `--bam` | required | Coordinate-sorted, indexed DEFND-seq gDNA BAM |
| `--bins` | required | Tab-separated bin file; header-driven columns `chrom`, `start`, `end`, `gc` |
| `--out-prefix` | required | Prefix for the two output `.txt.gz` files |
| `--barcode-tag` | `CB` | BAM tag holding the cell barcode (e.g. switch to `CR` for raw) |
| `--barcodes` | none | Optional barcode allowlist file; reads with barcodes not on it are skipped |
| `--mapq` | sane default (document exact value where implemented) | Minimum MAPQ to keep a read |
| `--exclude-dups` | on | Toggle exclusion of PCR/optical duplicate reads (SAM dup flag) |
| `--position` | `end` | Which read coordinate bins on: `start`, `midpoint`, or `end` |
| `--threads` | `1` | Worker thread count |
| `--help` | — | Print usage |

Reads are always dropped if: unmapped, secondary, supplementary, QC-fail, or
(when `--exclude-dups` is on) marked duplicate — these are not configurable
beyond the MAPQ threshold and the dup toggle. Reads without a `--barcode-tag`
tag, or excluded by the `--barcodes` allowlist, are also dropped.

## Core algorithm (single streaming pass)

- One pass over the coordinate-sorted BAM. For each read passing filters,
  compute its bin coordinate from `--position` (start/midpoint/end) and look
  up the containing bin by chromosome; reads outside all provided bins are
  dropped.
- Increment a counter for `(bin, cell_barcode)`.
- Bin identifiers in output are `<chr>:<start>:<end>` (colon-separated).
- Accumulate into a bins × barcodes counts structure sized from the bin file
  and the observed/allowlisted barcode set.

## Design decisions the implementation must make explicit

These are called out in the spec as choices to be made and *stated*, not
defaults to be picked silently:

- **Matrix storage: dense vs. sparse.** A dense matrix is fine at the scale
  of tens of thousands of bins × thousands of cells, but note the memory
  footprint in the implementation notes; move to a sparse accumulator if the
  barcode set gets large, and say so explicitly rather than assuming.
- **Barcode column ordering.** Pick exactly one of: allowlist order,
  first-seen order, or sorted order — and document which one was chosen and
  why, since it determines output column order.
- **Chromosome-naming consistency.** BAM header contig names and bin-file
  chrom names must be checked for a consistent convention (`chr1` vs `1`) up
  front. A mismatch must be a clear, loud error — never a silent drop of
  every read.
- **GC correction method and scope.** Use a standard, documented approach
  (LOESS/lowess of counts vs. bin GC, correcting each bin to the
  genome-wide median — the QDNAseq/Ginkgo-style approach). State explicitly
  whether correction is per-cell (the single-cell CNV norm, but more
  compute) or pooled across cells, and the tradeoff either way. If in-C++
  LOESS adds too much complexity, the accepted alternative is emitting raw
  counts plus the bin's GC column and doing correction in a small companion
  R/Python script — state which path was taken and why.

## Outputs

Two tab-separated, gzip-compressed `.txt.gz` files, both shaped bins (rows) ×
cells (columns), bin rows named `<chr>:<start>:<end>`, cell barcodes as
column headers:

1. `<out-prefix>...` raw integer bin counts per cell.
2. `<out-prefix>...` GC-corrected bin counts per cell.

## Constraints / non-goals

- No config files, no plugin system, no output formats beyond the two
  specified gzipped TSVs.
- Barcode handling only ever reads the BAM tag (`--barcode-tag`) — never
  parse barcodes out of read sequence.
- The tool assigns reads to bins the user supplies; it does not compute or
  merge bins itself.
- Validate inputs up front (BAM index present, bin file well-formed,
  chromosome-naming consistency) rather than failing silently mid-stream.
