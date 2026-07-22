# DEFND-seq bin × cell counting tool — design spec

Date: 2026-07-21

## Purpose

Turn a coordinate-sorted, indexed DEFND-seq gDNA BAM into a bin × cell count
matrix suitable for single-cell copy-number analysis: one raw-counts matrix
and one GC-bias-corrected matrix.

## Background

DEFND-seq (Olsen et al.) co-sequences RNA and DNA from single nuclei using
the 10x Genomics Chromium Single Cell Multiome ATAC + Gene Expression kit
"as-is," producing a gDNA library alongside the cDNA library. The gDNA BAM is
therefore a standard 10x-barcoded, coordinate-sorted BAM: the per-nucleus
cell barcode lives in a BAM tag (`CB` by default — the corrected barcode,
carrying a GEM-well suffix like `-1`), not in the read sequence. This is
plain per-read binning over gDNA alignments — no ATAC fragment logic.

## Related work: cna_utils

Reviewed `rishvanth-kp/cna_utils` (`scripts/getBinCounts.py`, plus its
`docs/cnv_pipeline_manual.pdf`) as a reference bin-counting implementation.

`getBinCounts.py` itself applies **no** MAPQ or SAM-flag filtering — it
iterates every alignment record and only excludes reads landing in a
"deadzone" BED, binning on `reference_start` (alignment start). That's
because the manual documents doing all filtering *upstream*, by hand, with
plain `samtools` before `getBinCounts.py` ever runs:
`samtools markdup -r` (dedup) → `samtools view -q 30 -F 0x800`
(MAPQ + drop supplementary) → **`samtools view -f 0x40`** (keep only read 1
of each pair, discard read 2 entirely — "Remove mates," done solely so a
sequenced fragment isn't double-counted once per mate).

Two things matter for us here:

1. Their pipeline is a generic paired-end single-cell WGS protocol
   (Baslan/Kendall-Krasnitz-style, not 10x/droplet-barcoded). Both mates are
   ordinary genomic sequence, so "keep read 1, drop read 2" is an arbitrary
   but harmless choice — either mate would do.
2. **Our situation differs**: DEFND-seq's gDNA library is 10x-chemistry
   based, and in 10x reads, read 1 conventionally carries the cell
   barcode/UMI rather than genomic sequence, while read 2 carries the actual
   insert. If the gDNA BAM was produced by aligning both mates directly
   (e.g. `bwa mem` on the R1/R2 fastqs without a Cell-Ranger-style
   barcode-extraction step first), read-1 alignment records are barcode
   sequence that happened to map somewhere — noise, not signal — so the
   mate worth keeping is the *opposite* of cna_utils's choice. This can't be
   assumed to have been handled upstream the way cna_utils assumes, since a
   DEFND-seq user is unlikely to have run an equivalent manual `samtools`
   pipeline first, so the filtering is built into `dfn_bin` itself (see CLI
   contract below) rather than documented as a prerequisite step.

Also, since `getBinCounts.py` never excludes secondary (`0x100`) or
supplementary (`0x800`) alignments itself, a multi-mapped or split-aligned
read would be double-counted in their tool if the upstream `-F 0x800` step
were ever skipped. Our design filters those unconditionally (see below) so
`dfn_bin` is correct standalone, without depending on an assumed upstream
pipeline.

## Scope / non-goals

- No config files, no plugin system, no output formats beyond the two
  specified gzipped TSVs.
- The tool assigns reads to bins the user supplies; it does not compute or
  merge bins itself.
- Barcode handling only ever reads the configured BAM tag — never parsed out
  of read sequence.
- Not ATAC-fragment-aware; simple per-read alignment binning.

## Architecture

Two artifacts, each with a single responsibility:

1. **`dfn_bin`** (C++17, htslib) — one streaming pass over the BAM, applies
   read filters, assigns reads to bins, emits a raw bin × cell count matrix.
2. **`correct_gc.py`** (Python) — reads the raw matrix + bin GC column, emits
   a GC-bias-corrected matrix.

No shared library between them; they communicate only through the raw-counts
file and the original bins file. This was chosen over implementing LOESS in
C++ so the binary's only job is fast, correct BAM streaming, and the
numerically fiddly part (weighted local regression) lives in ~80 lines of
well-tested Python (`statsmodels`) rather than hand-rolled C++.

## Assumed scale

Hundreds to low-thousands of cells, 10,000–50,000 bins per run. This governs
two choices below (dense matrix storage, naive per-cell LOESS runtime). If a
future run needs tens of thousands+ of cells, the dense matrix and the
per-cell LOESS loop are the two places that would need revisiting (sparse
accumulator; batched/parallel LOESS).

## Component 1: `dfn_bin`

### CLI contract

| Flag                       | Default  | Meaning                                                                     |
|----------------------------|----------|-----------------------------------------------------------------------------|
| `--bam`                    | required | Coordinate-sorted, indexed DEFND-seq gDNA BAM                               |
| `--bins`                   | required | Tab-separated bin file; header-driven columns `chrom`, `start`, `end`, `gc` |
| `--out-prefix`             | required | Prefix for `<prefix>.raw_counts.txt.gz`                                     |
| `--barcode-tag`            | `CB`     | BAM tag holding the cell barcode (e.g. `CR` for raw)                        |
| `--barcodes`               | none     | Optional barcode allowlist file; reads with barcodes not listed are skipped |
| `--mapq`                   | `30`     | Minimum MAPQ to keep a read                                                 |
| `--exclude-dups`           | true     | Exclude reads flagged as PCR/optical duplicates (`0x400`)                   |
| `--primary-alignment-only` | true     | Exclude secondary (`0x100`) and supplementary (`0x800`) alignment records   |
| `--second-mate-only`       | true     | For paired reads, keep only read 2 (`0x80`); drop read 1 (`0x40`). 10x reads carry the cell barcode/UMI on read 1, not genomic sequence — no-op on unpaired/single-end alignments |
| `--position`               | `end`    | Which read coordinate bins on: `start`, `midpoint`, or `end`                |
| `--threads`                | `1`      | Threads passed to htslib for BGZF decompression                             |
| `--help`                   | —        | Print usage                                                                 |

Reads are always dropped, regardless of other options, if unmapped (`0x4`)
or QC-fail (`0x200`) — these two are not configurable. On top of that, by
default (`--primary-alignment-only true`) only the primary alignment record
of a read is ever counted, never its secondary/supplementary records;
(`--second-mate-only true`) drops read-1 records of a mapped pair, since in
10x chemistry read 1 conventionally carries the barcode/UMI rather than
genomic sequence (see cna_utils comparison above — this is the reverse of
their "keep read 1" convention, deliberately, because the two pipelines'
read 1/2 semantics differ); and `--exclude-dups true` drops duplicate-flagged
(`0x400`) reads. All three are exposed as toggles (not hardcoded) so they
can be turned off if inspection of a real BAM shows a different read
structure than assumed here — this assumption should be verified against an
actual DEFND-seq gDNA BAM header/records before relying on the defaults.
Reads without a `--barcode-tag` tag, or excluded by `--barcodes`, are also
dropped. Reads whose bin coordinate falls outside every provided bin are
dropped.

Implementation-wise this is mask checks against `bam1_t::core.flag` using
htslib's existing flag constants (`BAM_FUNMAP`, `BAM_FSECONDARY`,
`BAM_FQCFAIL`, `BAM_FSUPPLEMENTARY`, `BAM_FDUP`, `BAM_FPAIRED`,
`BAM_FREAD1`) rather than hand-rolled bit literals.

### Algorithm

1. Parse bin file by header (order-independent columns `chrom`, `start`,
   `end`, `gc`); build a per-chromosome sorted vector of bins for
   binary-search lookup.
2. Validate: BAM index present (`.bai`/`.csi`); bin file has all required
   columns and valid numeric start/end; chromosome-naming consistency
   between BAM header contigs and bin-file chroms (see below).
3. If `--barcodes` given, load the allowlist; barcode column order is fixed
   to allowlist order and the matrix is sized/allocated up front.
   Otherwise, columns are assigned in first-seen order as new barcodes
   appear during the streaming pass (a `unordered_map<string,size_t>` from
   barcode to column index, columns appended to a growable structure).
4. Single streaming pass over the BAM (`bam_read1`). For each read: filter
   by FLAG mask (unmapped/QC-fail always; primary-alignment-only,
   second-mate-only, dup-exclusion per their toggles, per above) and MAPQ →
   extract barcode via `--barcode-tag` → check allowlist (if any) → compute
   bin coordinate from `--position` → binary search the chrom's bin vector →
   increment `counts[bin_idx][barcode_idx]` (or drop if no containing bin).
5. Write `<prefix>.raw_counts.txt.gz`: header row = `bin<TAB>barcode1<TAB>barcode2...`;
   each row = `<chr>:<start>:<end><TAB>count...`. Bin rows in bin-file order.

### Chromosome-naming validation

Take the set of contig names from the BAM header and the set of chrom names
from the bin file. If neither is a subset of the other, compare their `chr`
prefix convention (e.g. does the BAM use `chr1`-style and the bins use
`1`-style, or vice versa) and fail with an explicit error naming both
conventions found and how many contigs/bins are affected. No automatic
renaming — the user must fix one of the two inputs.

### Storage

Dense `int32_t` matrix (`bins.size() × barcodes.size()`), row-major by bin.
At the assumed scale this is at most a few hundred MB — acceptable, avoids
the complexity of a sparse accumulator.

### Threading

`--threads` is passed to `hts_set_threads()` for multi-threaded BGZF
decompression only. The counts matrix is populated by a single thread
reading records in order — no locking/merge logic needed, since decompression
I/O (not per-read filtering) is the actual bottleneck on a coordinate-sorted
BAM.

## Component 2: `correct_gc.py`

- Inputs: `--counts <prefix>.raw_counts.txt.gz`, `--bins bins.tsv`,
  `--out-prefix <prefix>`.
- For each cell (matrix column): fit `statsmodels.nonparametric.lowess`
  of raw counts against each bin's GC fraction, then rescale the fitted
  values so the corrected column's median equals that cell's genome-wide
  median raw count (QDNAseq/Ginkgo-style: correct each bin toward the
  genome-wide median, applied per cell — the single-cell CNV norm, at the
  cost of one LOESS fit per cell rather than one fit total).
- Output: `<prefix>.gc_corrected.txt.gz`, identical shape/labels to the raw
  matrix.

## Error handling

Fail fast with a specific, actionable message (not a stack trace or silent
empty output) for: missing BAM index, bin file missing a required header
column or containing malformed rows, chromosome-naming mismatch, an empty or
unreadable barcode allowlist file. None of these produce partial output.

## Build

Minimal `Makefile`: `g++ -std=c++17 $(pkg-config --cflags htslib) ... $(pkg-config --libs htslib)`.
No other dependencies. `correct_gc.py` requires `numpy`, `statsmodels`
(both standard, pip-installable).

## Testing

TDD. A fixture script builds a small synthetic BAM (few contigs, ~100–200
reads, a handful of barcodes, including reads that should be dropped by each
filter — low MAPQ, secondary/supplementary/dup flags, read-1-of-pair,
missing CB tag, off-allowlist barcode, outside all bins) plus a matching
small `bins.tsv`, using `pysam`/`samtools`. Unit/integration tests assert
exact expected counts per (bin, barcode) cell, plus the error-handling paths (bad index,
malformed bins file, chrom-naming mismatch). `correct_gc.py` is tested
against a small synthetic counts matrix with a known, injected GC trend,
asserting the trend is flattened post-correction.

## Validation (manual sanity checks post-implementation)

- Total counts per cell should be a large majority of that cell's
  passing-filter read count (most should land in a bin).
- Plotting raw count vs. bin GC should show a visible trend; the same plot
  post-correction should be flat.
