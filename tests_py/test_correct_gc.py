import gzip
import subprocess
import sys
import os
import numpy as np

REPO_ROOT = os.path.join(os.path.dirname(__file__), "..")
SCRIPT = os.path.join(REPO_ROOT, "correct_gc.py")


def write_counts(path, bins, barcodes, matrix):
    with gzip.open(path, "wt") as f:
        f.write("bin\t" + "\t".join(barcodes) + "\n")
        for bin_name, row in zip(bins, matrix):
            f.write(bin_name + "\t" + "\t".join(str(v) for v in row) + "\n")


def write_bins(path, bins, gcs):
    with open(path, "w") as f:
        f.write("chrom\tstart\tend\tgc\n")
        for bin_name, gc in zip(bins, gcs):
            chrom, start, end = bin_name.split(":")
            f.write(f"{chrom}\t{start}\t{end}\t{gc}\n")


def test_gc_trend_is_flattened(tmp_path):
    rng = np.random.default_rng(0)
    n_bins = 200
    gcs = rng.uniform(0.3, 0.6, size=n_bins)
    baseline = 100
    # Inject a strong linear GC trend plus small noise.
    counts = (baseline * (1 + 3.0 * (gcs - 0.45)) + rng.normal(0, 2, size=n_bins))
    counts = np.clip(counts, 0, None).astype(int)

    bins = [f"chr1:{i*1000}:{(i+1)*1000}" for i in range(n_bins)]
    barcodes = ["AAAA-1"]

    counts_path = str(tmp_path / "raw.txt.gz")
    bins_path = str(tmp_path / "bins.tsv")
    out_prefix = str(tmp_path / "out")

    write_counts(counts_path, bins, barcodes, counts.reshape(-1, 1))
    write_bins(bins_path, bins, gcs)

    result = subprocess.run(
        [sys.executable, SCRIPT, "--counts", counts_path, "--bins", bins_path,
         "--out-prefix", out_prefix],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr

    with gzip.open(out_prefix + ".gc_corrected.txt.gz", "rt") as f:
        lines = [line.rstrip("\n").split("\t") for line in f]
    header, rows = lines[0], lines[1:]
    assert header == ["bin", "AAAA-1"]
    corrected = np.array([int(float(row[1])) for row in rows])

    raw_corr = abs(np.corrcoef(counts, gcs)[0, 1])
    corrected_corr = abs(np.corrcoef(corrected, gcs)[0, 1])
    assert raw_corr > 0.5          # confirms the injected trend is real
    assert corrected_corr < 0.15   # confirms correction flattened it


def test_shape_and_labels_match_raw_matrix(tmp_path):
    bins = ["chr1:0:1000", "chr1:1000:2000"]
    barcodes = ["AAAA-1", "BBBB-1"]
    counts = np.array([[10, 20], [15, 5]])
    gcs = [0.4, 0.5]

    counts_path = str(tmp_path / "raw.txt.gz")
    bins_path = str(tmp_path / "bins.tsv")
    out_prefix = str(tmp_path / "out")

    write_counts(counts_path, bins, barcodes, counts)
    write_bins(bins_path, bins, gcs)

    result = subprocess.run(
        [sys.executable, SCRIPT, "--counts", counts_path, "--bins", bins_path,
         "--out-prefix", out_prefix],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr

    with gzip.open(out_prefix + ".gc_corrected.txt.gz", "rt") as f:
        lines = [line.rstrip("\n").split("\t") for line in f]
    assert lines[0] == ["bin", "AAAA-1", "BBBB-1"]
    assert [row[0] for row in lines[1:]] == bins
