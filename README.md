# dfn_bin

Turns a coordinate-sorted, indexed DEFND-seq gDNA BAM into a bin×cell
fragment-count matrix for single-cell copy-number analysis, plus a
companion BAM of structural-rearrangement evidence. Full design rationale:
`docs/superpowers/specs/2026-07-21-defnd-seq-binning-design.md`.

## Build

Requires htslib (`brew install htslib` / `apt install libhtslib-dev`) and a
C++17 compiler.

```bash
make
```

Produces `./dfn_bin`. Run `make test` to run the C++ unit test suite.

## Usage

```bash
./dfn_bin --bam sample.bam --bins bins.tsv --out-prefix sample
python3 correct_gc.py --counts sample.raw_counts.txt.gz --bins bins.tsv --out-prefix sample
```

**Note:** `correct_gc.py` requires `numpy` and `statsmodels`. Install them with:
```bash
pip install numpy statsmodels
```

`bins.tsv` is a tab-separated file with a header row containing (in any
column order) `chrom`, `start`, `end`, `gc`.

### Outputs

- `sample.raw_counts.txt.gz` — bins (rows) × cells (columns), raw integer
  fragment counts. Bin rows are named `<chrom>:<start>:<end>`; columns are
  cell barcodes.
- `sample.gc_corrected.txt.gz` — same shape, GC-bias-corrected (via
  `correct_gc.py`).
- `sample.discordant.bam` (+ `.bai`) — split/chimeric reads and discordant
  read pairs (cross-chromosome, or same-chromosome with the aligner's
  proper-pair flag unset), for manual/downstream structural-variant review.
  Excluded from both counts matrices.

See `./dfn_bin --help` for the full flag list and defaults.

## Validation

After a run, sanity-check the output before trusting it:

- **Assignment rate**: total counts per cell should be a large majority of
  that cell's passing-filter read count. A low fraction suggests a bins
  file that doesn't cover most of the genome, or a chromosome-naming issue
  that silently passed validation.
- **GC-correction check**: plot raw count vs. bin GC fraction — it should
  show a visible trend. Plot the same thing post-correction (using
  `sample.gc_corrected.txt.gz`) — the trend should be flat.
- **Discordant rate**: reads landing in `sample.discordant.bam` should be a
  small minority of total reads. A large fraction suggests a library or
  alignment QC problem rather than genuine genome-wide rearrangement.
  Spot-check a handful of entries in IGV.

## Testing

```bash
make test                                    # C++ unit tests
.venv/bin/pytest tests/integration -v        # end-to-end integration tests
.venv/bin/pytest tests_py -v                 # correct_gc.py tests
```

On a fresh clone, set up the Python virtualenv with:

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
```

## License

BSD 3-Clause. See [LICENSE](LICENSE).
