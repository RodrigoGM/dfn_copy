# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Status

This repo ships two working C++17 binaries, `dfn_copy` and `dfn_cbs`, built from one root `Makefile`, currently at version 0.2.0 (see `VERSION`/`changelog.txt`). Both binaries, their full design rationale, and their test suites already exist — this file documents what's actually there, not a target design. Update it whenever the CLI contract, module list, or architecture changes; don't let it drift back into aspirational claims.

## Build

```bash
make            # builds both ./dfn_copy and ./dfn_cbs
make clean      # removes build/, tests/build/, and both binaries
```

`dfn_copy` requires htslib (`brew install htslib` / `apt install libhtslib-dev`); `dfn_cbs` has **no** htslib dependency by design (it never touches a BAM) — its Makefile rules deliberately use explicit compile flags instead of the shared `$(CXXFLAGS)`/`$(LDLIBS)`, which carry htslib's pkg-config output, so `dfn_cbs` alone is buildable without htslib installed at all. Preserve that separation when adding rules for either binary.

## Testing

```bash
make test                                     # C++ unit tests, one binary per module
.venv/bin/pytest tests/integration -v         # end-to-end integration tests (both binaries)
.venv/bin/pytest tests_py -v                  # correct_gc.py tests
```

Fresh-clone bootstrap for the Python side:
```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
```

Unit tests use a hand-rolled framework (`tests/test_framework.hpp`: `ASSERT_TRUE`/`ASSERT_EQ`/`ASSERT_NEAR`/`TEST_REPORT()`); BAM-specific test helpers live separately in `tests/bam_test_helpers.hpp` so `dfn_cbs`'s tests never need htslib. Integration fixtures for `dfn_copy` are pysam-generated BAMs (`tests/fixtures/`); `dfn_cbs`'s fixtures are plain gzip/text (`tests/fixtures/gen_cbs_fixture.py`), since its input is already a matrix.

## Architecture

### `dfn_copy` — BAM → bin × cell counts matrix

One streaming pass over a coordinate-sorted, indexed DEFND-seq gDNA BAM: filters reads, resolves paired fragments to one representative alignment, bins them, and emits a raw bin × cell counts matrix plus a diverted BAM of structural-variant evidence.

Modules (`src/`): `cli_args`, `bins` (bin-file parsing + coordinate lookup + chrom-naming validation), `barcode_index`, `counts_matrix` (dense growable matrix, gzip TSV writer), `read_filter`, `fragment_pairing` (QNAME-keyed lookback buffer, `auto`/`r1`/`r2` modes), `discordant_writer` (RAII BAM writer for split/discordant reads), `main.cpp` (orchestration).

CLI: `--bam`, `--bins`, `--out-prefix` (required); `--barcode-tag` (`CB`); `--barcodes`; `--mapq` (`30`); `--exclude-dups` (`true`); `--primary-alignment-only` (`true`); `--single-end-counting` (`auto`|`r1`|`r2`); `--max-insert-size` (`0`, disabled); `--position` (`start`|`midpoint`|`end`, default `end`); `--threads` (`1`, htslib BGZF decompression only); `--help`.

Outputs: `<prefix>.raw_counts.txt.gz` (bins × cells, gzip TSV), `<prefix>.discordant.bam` (+`.bai`) for split/chimeric reads and discordant pairs.

Full design rationale: `docs/superpowers/specs/2026-07-21-defnd-seq-binning-design.md`.

### `dfn_cbs` — GC correction + CBS segmentation at scale

Reads a `dfn_copy` raw-counts matrix, GC-corrects it per cell (LOWESS reimplemented in C++, matching `correct_gc.py`'s method/statsmodels defaults), and segments each cell's corrected profile via an **in-process vendored copy** of the user's separate `cbs+` project (circular binary segmentation, source at `~/mskcc/dev/cbs+/`) — copied under `src/vendor/cbs/`, never subprocess-called, `cbs+` itself never modified. Built because `correct_gc.py` OOM'd on a real 6,866-bin × 475,674-cell matrix; `dfn_cbs`'s whole point is a bounded memory model at that scale.

Three-phase, bounded-memory architecture: **Phase 1** loads the raw-counts matrix once (freed after Phase 2). **Phase 2** processes surviving cells in parallel (one worker per cell), GC-correcting and segmenting each, writing results into three cell-major scratch files via lock-free positioned I/O (`pwrite`/`pread`, disjoint byte ranges per cell — no locking needed) plus appending `.seg` rows. **Phase 3** transposes each scratch file to the final `bins × cells` gzip TSV convention, one file in memory at a time — never two full matrices resident simultaneously. Peak RSS stays near one matrix's size (~13–15GB at the motivating scale) regardless of cell count; scratch disk space needed is roughly 3× one output matrix's size, held for the duration of the run.

CBS segmentation always runs with **centering forced off** (not exposed as a flag) — normal CBS per-chromosome centering would erase real whole-chromosome copy-number differences, exactly the signal GC correction is meant to preserve. Per-cell/per-chromosome RNG seed is `--seed XOR hash(barcode) XOR hash(chrom)`, so results are independent of worker thread scheduling.

Modules (`src/`): `cbs_args`, `bin_gc` (bin/GC loader, distinct from `dfn_copy`'s `bins.hpp` — no coordinate-lookup machinery needed), `raw_counts_loader`, `cell_filter`, `lowess`, `gc_correct`, `segment_cell` (chromosome grouping + vendored-CBS call + coordinate mapping), `scratch_matrix` (positioned-I/O cell-major scratch files), `seg_file_writer`, `transpose_writer`, `progress_reporter` (stderr status line with ETA, TTY-aware), `cbs_main.cpp` (orchestration), `vendor/cbs/` (vendored CBS core — see below).

CLI: `--counts`, `--bins`, `--out-prefix` (required); `--min-reads` (`100000` — cells below this raw-count sum are dropped entirely from every output); `--threads` (`1`, cell-level parallelism — different meaning from `dfn_copy`'s `--threads`); `--alpha` (`0.01`); `--perms` (`1000`); `--min-seg-len` (`25`); `--max-depth` (`100`); `--cbs-method` (`1cp`|`2cp`); `--seed` (`1`); `--quiet`; `--help`.

Outputs: `<prefix>.gc_corrected.txt.gz`, `<prefix>.lowess_ratio.txt.gz`, `<prefix>.segmented_lowess_ratio.txt.gz` (all bins × cells gzip TSV), `<prefix>.seg` (plain text, IGV/DNAcopy format — not gzipped).

Full design rationale: `docs/superpowers/specs/2026-07-22-dfn_cbs-gc-correction-segmentation-design.md`.

### `correct_gc.py`

Python GC-correction fallback, kept for users with small matrices who don't want to build the C++ tool chain. `dfn_cbs` supersedes it at production scale but does not replace it.

## Vendored code

`src/vendor/cbs/` is a deliberate **verbatim copy** (not a submodule or symlink) of `~/mskcc/dev/cbs+/`'s segmentation core, with exactly one intentional change (`cbs.hpp`'s include of a trimmed local `io_types.hpp` instead of `cbs+`'s real `io.hpp`, which pulls in unrelated CLI/file-parsing declarations `dfn_cbs` doesn't need). See `src/vendor/cbs/README.md` for the full provenance note. Never "clean up" this code's style to match the rest of the repo, and never strip its unused paths (e.g. `correction.cpp`'s rolling/linear correction, never exercised since `dfn_cbs` always leaves `Args::correction = "none"`) — fidelity to the tested upstream source is the point. If `cbs+` changes upstream, this copy needs updating by hand.

## Constraints / non-goals

- `dfn_copy` assigns reads to bins the user supplies; it does not compute or merge bins itself, and is not ATAC-fragment-aware (no Tn5 cut-site offset correction).
- Barcode handling only ever reads the BAM tag (`--barcode-tag`) — never parsed out of read sequence.
- `dfn_cbs` never opens a BAM and has no htslib dependency, enforced at the build level.
- No config files, no plugin system, no output formats beyond what's documented above for each binary.
- Validate inputs up front (BAM index present, bin file well-formed, chromosome-naming consistency, `--counts`/`--bins` matching pair) rather than failing silently mid-stream — every fallible operation throws a specific `std::runtime_error`, never crashes or returns a sentinel.
