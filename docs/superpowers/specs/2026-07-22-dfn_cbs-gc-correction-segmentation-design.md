# dfn_cbs — GC correction + CBS segmentation at scale — design spec

Date: 2026-07-22

## Purpose

`correct_gc.py` (per-cell LOWESS GC correction in Python/statsmodels) crashed
out of memory on a real run: 6,866 bins × 475,674 cells. `dfn_cbs` is a new
C++17 binary that does the same GC correction at that scale without the
memory blowup, and additionally segments each cell's corrected profile via
circular binary segmentation (CBS), producing copy-number-ready outputs
directly: a GC-corrected counts matrix, the underlying correction ratio, a
segmented version of that ratio, and an IGV-loadable `.seg` file.

`correct_gc.py` is not removed — it stays as a documented fallback for users
with small matrices who don't want to build the C++ tool chain. `dfn_cbs` is
the path for real-scale data.

## Background

Two things came together to make the crash predictable in hindsight, not a
surprise to work around blindly:

1. **The correction is inherently per-cell.** A LOWESS fit of counts against
   bin GC only makes sense within one cell's own profile — you cannot fit
   across cells. So the *computation* itself is naturally small (one
   6,866-length vector at a time) and cheap; nothing about the correction
   algorithm requires large working memory.
2. **The crash was a materialization problem, not a computation problem.**
   Python/numpy/statsmodels each hold their own copies of intermediate
   arrays (the full matrix, per-column slices, LOWESS's internal sorted
   copies, the output matrix) with no control over when those copies are
   freed. At this scale — a single `bins × cells` `float64` matrix is
   ~26GB — a handful of simultaneous incidental copies is enough to exceed
   even a generous 64GB budget.

`dfn_cbs` fixes this by making the memory model explicit and bounded
end-to-end (see Memory model below), and by doing the numerically fiddly
parts (LOWESS, CBS) in C++ where memory ownership is under our control.

## Related work: cbs+

The user's existing `cbs+` (`~/mskcc/dev/cbs+/`) is a standalone C++17
circular binary segmentation tool. Its CLI (`bin/cbs_tool`) segments **one
value series at a time** — `--in chrom,pos,value[,GC,WAVE]` — with its own
per-chromosome thread pool. It has no native multi-cell/matrix concept, so
invoking it as a subprocess once per cell (475,674 times) would mean 475k
process spawns plus a temp-file round trip per cell — a real bottleneck
independent of the segmentation work itself.

Its actual segmentation logic (`stats.cpp`'s `max_cusum_1cp`/`max_cusum_2cp`,
`cbs.cpp`'s `segment_chromosome`) is a small set of plain, dependency-free
functions with no global mutable state — safe to call directly, concurrently,
from multiple threads, each with its own `std::mt19937_64` RNG instance. This
is vendored into `dfn_cbs`'s own source tree (copied, not subprocess-called,
not built as a shared library) so `cbs+` itself stays untouched as a
general-purpose standalone tool, while `dfn_cbs` gets in-process,
zero-spawn-overhead segmentation.

`cbs+` also ships its own GC-bias hooks (`correction.cpp`:
`linear_covariate_regress`, `rolling_median_detrend`) — a linear regression
against a GC covariate, and a rolling-median detrend by genomic position.
Neither is used here: GC bias is typically non-monotonic, which is exactly
why LOWESS (not a linear fit) is the standard method (QDNAseq/Ginkgo-style),
and rolling-median detrend isn't a GC correction at all (it detrends by
position, not GC content). `dfn_cbs` reimplements LOWESS in C++ instead,
matching `correct_gc.py`'s existing, already-validated method.

One more `cbs+` behavior is deliberately **not** carried over: its `--center`
option (per-chromosome median-centering of the value series before
segmentation, on by default in `cbs+`). For single-cell CNV, this would
actively erase the signal the GC correction is designed to preserve — a
real whole-chromosome gain or loss is exactly a shift in that chromosome's
own median relative to the rest of the genome, and centering removes it
before CBS ever sees the data. `dfn_cbs` calls the vendored segmentation
core with centering forced off (not exposed as a flag — enabling it would
be actively wrong for this use case, so it isn't offered as an option).

## Scope / non-goals

- Reads only the matrix format `dfn_copy` already produces
  (`<prefix>.raw_counts.txt.gz`) and the same `bins.tsv` format `dfn_copy`
  consumes. `dfn_cbs` never touches a BAM.
- No config files, no plugin system, no output formats beyond the four
  specified below.
- Not a general-purpose segmentation tool — `cbs+` remains that; `dfn_cbs`
  vendors only the specific segmentation call it needs, with fixed
  (non-exposed) centering behavior appropriate to single-cell CNV.
- No change to `cbs+` itself.
- `correct_gc.py` is retained as-is for small-matrix use; `dfn_cbs` does not
  replace it, only supersedes it at scale.
- CBS's circular-permutation testing (default 1,000 permutations per
  candidate split) is inherently expensive per cell; at 475k cells this
  dominates runtime. `--threads` (cell-level parallelism) and `--perms` are
  the available levers — no algorithmic change to CBS's permutation testing
  is in scope here.

## Architecture

One binary, three phases, each phase's memory footprint bounded
independently so peak RSS never requires holding more than one
`bins × cells`-sized structure at a time:

**Phase 1 — load.** Read `<prefix>.raw_counts.txt.gz` fully into memory as a
dense `int32_t` `bins × cells` array (row-major by bin, matching the file's
own layout) — one read, ~13GB at this scale. Read `bins.tsv` for GC values
and (chrom, start, end) per bin, in bin-file order (matching the counts
matrix's row order). Sum each column (cell) to get its total raw count.

**Phase 2 — per-cell compute, in parallel.** Drop any cell whose column sum
is below `--min-reads`. For each surviving cell, on a worker thread from a
pool sized by `--threads`:

1. Extract that cell's column into a contiguous `bins`-length vector
   (`counts`) — trivial, ~27KB.
2. Fit LOWESS of `counts` against the bins' GC values (see LOWESS below) →
   `fitted`. Any non-positive fitted value is clamped to machine epsilon
   before use as a divisor (same safeguard as `correct_gc.py`, avoiding
   divide-by-zero or a sign-flipped ratio at a bin where the fit dips to or
   below zero).
3. Compute `lowess_ratio = counts / fitted` (mean ≈ 1, this cell's
   GC-corrected ratio track) and `gc_corrected = counts * (median(counts) /
   fitted)` (same shape, rescaled back to raw-count units — identical
   definition to `correct_gc.py`'s existing output). If `median(counts) ==
   0`, both outputs fall back to the raw `counts` unchanged (same safeguard
   as `correct_gc.py`).
4. Group `lowess_ratio` by chromosome (using the bin file's chrom column,
   in the order chromosomes first appear there) and run the vendored CBS
   core once per chromosome, centering forced off, to get segment
   boundaries and means.
5. Expand each chromosome's segments back into a `bins`-length vector
   (`segmented_lowess_ratio`) — every bin gets its segment's mean.
6. Emit this cell's three result vectors, plus one `.seg` row per segment
   (`barcode`, `chrom`, `start`, `end`, `num_bins`, `seg.mean`) — `start`/
   `end` are the genomic coordinates of the segment's first and last bin,
   `seg.mean` is the segment's mean `lowess_ratio` value.
7. Write the three result vectors into three **cell-major scratch files**
   (flat binary, `float32`, fixed record size = `bins × 4` bytes, opened in
   random-access mode so each cell writes to its own byte offset —
   `cell_index × bins × 4` — safely and independently of other threads,
   no locking needed since ranges never overlap). Append this cell's `.seg`
   rows directly to the final `.seg` output (small file; a short
   mutex-guarded append, negligible contention since it happens once per
   cell, not once per bin).
8. Free the cell's vectors immediately; nothing about this phase holds more
   than `O(bins × threads)` beyond the (already-loaded) input matrix.

Free the input matrix once Phase 2 completes — it isn't needed again.

**Phase 3 — transpose on exit.** For each of the three scratch files, in
turn (never more than one at a time): load it fully into memory (~13GB,
`float32`, cell-major), then write the final `<prefix>.*.txt.gz` in
`bins × cells` order — outer loop over bins, inner loop over surviving
cells, a strided read into the loaded scratch array, written directly to
the gzip stream (no second full-size buffer — the transpose is a read
pattern, not a copy). Free the scratch array before loading the next one.
Delete the scratch files once all three are written.

Peak RSS is therefore ≈ max(input matrix during Phase 2, one scratch array
during Phase 3) — about 13–15GB with overhead, regardless of how many of
the three output matrices there are, and does not grow if cell count grows
further (only the *scratch file sizes* grow; they're never held alongside
each other).

## LOWESS implementation

Reimplemented in C++ to match `statsmodels.nonparametric.smoothers_lowess.
lowess`'s defaults for behavioral parity with `correct_gc.py`: `frac = 2/3`
(fraction of points used in each local regression — a *fraction*, so it
scales naturally regardless of bin count, no special-casing needed), `iter =
3` robustifying (bisquare-weighted) iterations, tri-cube kernel weighting,
local linear (not constant or quadratic) fit at each point. Returns fitted
values at each original bin's GC value (`return_sorted=False` equivalent),
not resampled onto a grid.

## CLI contract

| Flag             | Default | Meaning                                                                 |
|-------------------|--------|--------------------------------------------------------------------------|
| `--counts`        | required | Path to `<prefix>.raw_counts.txt.gz` (from `dfn_copy`)                |
| `--bins`          | required | Same `bins.tsv` used to produce `--counts` (header columns `chrom`, `start`, `end`, `gc`) |
| `--out-prefix`    | required | Prefix for all four output files                                     |
| `--min-reads`     | `100000` | Minimum sum of raw binned counts for a cell to be kept; cells below this are dropped entirely from every output |
| `--threads`       | `1`      | Worker threads for cell-level parallelism in Phase 2 (not htslib decompression — `dfn_cbs` never opens a BAM) |
| `--alpha`         | `0.01`   | CBS split-acceptance p-value threshold (passed through to the vendored core) |
| `--perms`         | `1000`   | CBS circular-permutation count per candidate split                     |
| `--min-seg-len`   | `25`     | CBS minimum segment length (bins)                                      |
| `--max-depth`     | `100`    | CBS maximum recursive split depth                                      |
| `--cbs-method`    | `1cp`    | CBS split search: `1cp` (single changepoint) or `2cp` (two-changepoint scan) |
| `--seed`          | `1`      | RNG seed; per-cell RNG is seeded `seed XOR hash(barcode)` for reproducibility, mirroring `cbs+`'s own per-chromosome seeding convention |
| `--help`          | —        | Print usage                                                             |

`--center` is intentionally not exposed (see Related work: cbs+ above).

## Outputs

Four files:

1. `<prefix>.gc_corrected.txt.gz` — bins × cells, gzip TSV, raw-count-scale
   GC-corrected counts. Same shape/labels/convention as `correct_gc.py`'s
   output of the same name.
2. `<prefix>.lowess_ratio.txt.gz` — bins × cells, gzip TSV, dimensionless
   ratio (mean ≈ 1 per cell) — the value CBS segments.
3. `<prefix>.segmented_lowess_ratio.txt.gz` — bins × cells, gzip TSV, same
   units as (2) but with each bin replaced by its segment's mean.
4. `<prefix>.seg` — plain text (not gzipped — it's a row-per-segment file,
   not a matrix, and IGV loads `.seg` directly without decompression),
   columns `ID  chrom  loc.start  loc.end  num.mark  seg.mean` (standard
   IGV/DNAcopy `.seg` format), `ID` = cell barcode.

Cells dropped by `--min-reads` are absent from all four outputs — not
zeroed, not flagged, simply not present as columns (or rows, for `.seg`).

## Error handling

Fail fast with a specific, actionable message (not a stack trace or silent
empty/partial output) for: `--counts` file missing, malformed, or with a
header that doesn't match `--bins`' bin count/order; `--bins` missing a
required column or containing malformed rows; every cell filtered out by
`--min-reads` (report the threshold and the observed max column sum rather
than silently emitting empty matrices); scratch-file creation failure in
`<out-prefix>`'s directory (same class of error as `dfn_copy`'s
`--out-prefix` directory checks).

## Build

Extends the existing `Makefile` with a second binary target (`dfn_cbs`),
built from its own `src/` files plus the vendored `cbs+` core
(`stats.cpp`/`cbs.cpp`, copied under e.g. `src/vendor/cbs/`, with a comment
noting the origin and that it's a deliberate copy, not a symlink or
submodule, matching `cbs+`'s own build having no library target to link
against). No new external dependencies beyond what `dfn_copy` already
requires (zlib for gzip output); no htslib dependency for this binary since
it never touches a BAM.

## Testing

TDD, following the same pattern as `dfn_copy`. Unit tests: LOWESS against a
small synthetic (GC, count) series with a known injected trend, asserting
the fit flattens it (same style of assertion `correct_gc.py`'s existing test
uses); the vendored CBS call against a synthetic ratio series with a known
injected breakpoint, asserting a segment boundary lands at the right bin;
the cell-major scratch-file read/write round trip; the bin-major transpose
-on-exit logic against a small known matrix. Integration tests: a small
synthetic `raw_counts.txt.gz` + `bins.tsv` fixture (few chromosomes, a
handful of cells spanning both sides of `--min-reads`, at least one cell
with an injected whole-chromosome copy-number shift to verify it survives
segmentation without being centered away) asserting all four output files'
shapes, labels, and that the filtered-out cell is absent from every output.

## Validation (manual sanity checks post-implementation)

- `lowess_ratio`'s per-cell mean should be ≈1; plotting `gc_corrected` vs.
  bin GC should be flat (same check as `correct_gc.py`'s existing
  validation guidance).
- `segmented_lowess_ratio` should visibly flatten within-segment noise
  relative to `lowess_ratio` while preserving between-segment
  (whole-chromosome-scale) differences — spot-check a cell with a known
  aneuploidy in IGV against its `.seg` row.
- Peak RSS during a real run should stay near the ~13–15GB figure above,
  not scale with cell count beyond that — worth confirming directly (e.g.
  `/usr/bin/time -l` or equivalent) on the first real-scale run.
