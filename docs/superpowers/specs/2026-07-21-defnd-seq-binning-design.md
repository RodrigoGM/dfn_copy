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
carrying a GEM-well suffix like `-1`), not in the read sequence. Binning is
per-fragment (a mapped pair is consolidated to one representative alignment
before binning, see `--single-end-counting`), but there is no true
ATAC-style fragment-interval logic (no Tn5 cut-site offset correction, no
binning by the full insert span) — each counted fragment is binned by a
single coordinate (`--position`) taken from its one representative read.

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
2. **Our situation differs, confirmed against real DEFND-seq BAM records**:
   both mates of a pair carry *identical* `CR`/`CB` tags (e.g. FLAG 83/163
   pair at `chr1:9998`, both tagged `CB:Z:GGCTCACAGTAATCCA-1`). That's the
   signature of 10x **ATAC-style** chemistry (the barcode comes from a
   separate index read and is stamped onto both genomic mates), not 3′-GEX
   chemistry (where read 1 itself would be the barcode/UMI). So — unlike an
   earlier draft of this doc assumed — read 1 is *not* barcode noise here;
   both mates are real genomic sequence, symmetric in validity, and either
   one can be the correct fragment representative depending on which
   survives filtering or which alignment is longer. Observed read lengths
   also vary within a pair (e.g. 27bp vs 51bp, occasional soft-clipping),
   so a fixed "always drop read 1" rule (cna_utils's choice, harmless for
   *their* symmetric WGS mates) would be actively wrong here — it could
   discard the better-aligned mate of a pair, or the only one that passed
   filtering at all.

   What *does* carry over from cna_utils is the underlying reason they
   prune to one alignment per pair in the first place: **the count is
   meant to represent one fragment/molecule, not one alignment record** —
   counting both mates of a pair double-counts the same physical DNA
   molecule. `dfn_bin` needs to do this pruning itself (see
   `--single-end-counting` below), rather than assume it happened upstream
   the way cna_utils's manual does, since a DEFND-seq user is unlikely to
   have run an equivalent manual `samtools` pipeline first.

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
- No Tn5/ATAC-style fragment-interval logic (cut-site offsets, full-insert
  binning) — fragments are consolidated to one representative read (see
  `--single-end-counting`) and binned by that read's own coordinate.

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
| `--single-end-counting`    | `auto`   | Fragment-counting mode for paired data: `auto`, `r1`, or `r2` (see below)   |
| `--position`               | `end`    | Which read coordinate bins on: `start`, `midpoint`, or `end`                |
| `--threads`                | `1`      | Threads passed to htslib for BGZF decompression                             |
| `--help`                   | —        | Print usage                                                                 |

Reads are always dropped, regardless of other options, if unmapped (`0x4`)
or QC-fail (`0x200`) — these two are not configurable. On top of that, by
default (`--primary-alignment-only true`) only the primary alignment record
of a read is ever counted, never its secondary/supplementary records; and
`--exclude-dups true` drops duplicate-flagged (`0x400`) reads. Both are
exposed as toggles, not hardcoded, so they can be turned off if inspection
of a real BAM shows a different structure than assumed here. Reads without
a `--barcode-tag` tag, or excluded by `--barcodes`, are also dropped. Reads
whose bin coordinate falls outside every provided bin are dropped.

Implementation-wise the always-on/toggle filters above are mask checks
against `bam1_t::core.flag` using htslib's existing flag constants
(`BAM_FUNMAP`, `BAM_FSECONDARY`, `BAM_FQCFAIL`, `BAM_FSUPPLEMENTARY`,
`BAM_FDUP`) rather than hand-rolled bit literals.

### `--single-end-counting`: counting fragments, not alignment records

The count in each bin is meant to represent the number of DNA
fragments/molecules observed, not the number of alignment records — a
mapped pair must contribute exactly one count, not two, regardless of which
physical read (1 or 2) happens to be the better/surviving one. This flag
controls how a mapped pair collapses to its one representative alignment
before binning; it has no effect on already-unpaired/single-end records
(`BAM_FPAIRED` unset), which are always counted individually since there's
no pair to collapse.

- **`auto`** (default): per fragment, count whichever surviving primary
  mate has the longer reference-aligned span (`bam_endpos - pos`, from the
  CIGAR); ties break toward read 1. If only one mate of a pair survives the
  other filters (MAPQ, dup, etc.), that one is counted — the fragment is
  not dropped just because its partner didn't make it. This is the
  appropriate default for DEFND-seq's 10x ATAC-style gDNA reads, where
  either physical read number can be the informative one (see cna_utils
  comparison above) and read length can vary between mates of the same
  pair.
- **`r1`** / **`r2`**: always use a fixed mate (`BAM_FREAD1` / second-in-pair)
  as the fragment's representative, matching cna_utils's convention for a
  plain (non-10x) paired-end WGS library where either mate is equally
  valid and a fixed choice is simplest. Strict: if that specific mate
  didn't survive the other filters, the fragment is not counted at all
  (matches cna_utils's `-f 0x40` behavior exactly, for users who want
  reproducibility against that convention).

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
   dup-exclusion per its toggle, per above) and MAPQ → extract barcode via
   `--barcode-tag` → check allowlist (if any). If the read is unpaired
   (`BAM_FPAIRED` unset), proceed straight to binning. If paired, resolve it
   through the fragment-pairing buffer (below) per `--single-end-counting`
   before binning — a read is only ever binned once its fragment
   representative has been decided.
5. For a read cleared to be counted: compute bin coordinate from
   `--position` on that alignment's own coordinates → binary search the
   chrom's bin vector → increment `counts[bin_idx][barcode_idx]` (or drop
   if no containing bin).
6. Write `<prefix>.raw_counts.txt.gz`: header row = `bin<TAB>barcode1<TAB>barcode2...`;
   each row = `<chr>:<start>:<end><TAB>count...`. Bin rows in bin-file order.

### Fragment-pairing buffer

Mates of a pair aren't adjacent in a coordinate-sorted BAM, but for typical
gDNA fragment sizes (tens to a few hundred bp — well under bin size) they're
always close. Implementation: an `unordered_map<QNAME, PendingMate>` holding
reads whose partner hasn't been seen yet (`PendingMate` = alignment
coordinates, reference-aligned length, and enough state to re-derive its
bin). Each mapped, paired, primary/filter-passing read is looked up by QNAME:

- **Mate already pending** → resolve the pair now per `--single-end-counting`
  (`auto`: keep the longer span, tie→read1; `r1`/`r2`: keep the fixed mate,
  drop the fragment if the required mate isn't the one on hand), bin the
  winner, erase the map entry.
- **Mate not pending** → store this read, keyed by QNAME, and continue.
- **Flushing unresolved entries**: each stored read carries its own mate's
  expected coordinates (`RNEXT`/`PNEXT`, already in the BAM record). Once
  the stream's current position passes that expected coordinate (same
  chromosome) without the mate having appeared — it was filtered out
  upstream (MAPQ, dup, off-allowlist, etc.) — the held read is flushed as a
  singleton: counted in `auto` mode (a fragment shouldn't vanish just
  because its partner didn't survive), dropped in `r1`/`r2` mode if it
  isn't the required mate. Cross-chromosome mates (discordant pairs) are
  rare post-MAPQ-filtering in practice and are flushed as singletons the
  same way once their chromosome's records are exhausted. At end-of-file,
  any remaining pending entries are flushed the same way.

Memory footprint is bounded by the number of *unresolved* fragments within
one insert-size window at any point in the stream — small relative to the
main dense counts matrix.

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
the complexity of a sparse accumulator. The fragment-pairing buffer adds a
small, bounded amount on top (see above).

### Threading

`--threads` is passed to `hts_set_threads()` for multi-threaded BGZF
decompression only. The counts matrix and the fragment-pairing buffer are
populated by a single thread reading records in order — pairing resolution
is inherently order-dependent (it relies on coordinate-sorted proximity
between mates), so no locking/merge logic is introduced; decompression I/O,
not per-read filtering or pairing, is the actual bottleneck on a
coordinate-sorted BAM.

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
filter — low MAPQ, secondary/supplementary/dup flags, missing CB tag,
off-allowlist barcode, outside all bins) plus a matching small `bins.tsv`,
using `pysam`/`samtools`. Unit/integration tests assert exact expected
counts per (bin, barcode) cell, plus the error-handling paths (bad index,
malformed bins file, chrom-naming mismatch). Fragment-pairing gets its own
fixture cases: complete pairs with equal-length mates, complete pairs with
unequal-length mates (assert `auto` keeps the longer one), pairs where only
one mate survives other filters (assert it's still counted in `auto`, and
correctly kept/dropped per `r1`/`r2`), and a discordant cross-chromosome
pair (assert both ends are handled as singletons once flushed).
`correct_gc.py` is tested against a small synthetic counts matrix with a
known, injected GC trend, asserting the trend is flattened post-correction.

## Validation (manual sanity checks post-implementation)

- Total counts per cell should be a large majority of that cell's
  passing-filter read count (most should land in a bin).
- Plotting raw count vs. bin GC should show a visible trend; the same plot
  post-correction should be flat.
