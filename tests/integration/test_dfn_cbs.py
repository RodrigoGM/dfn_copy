import gzip
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "fixtures"))
import gen_cbs_fixture

REPO_ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
DFN_CBS = os.path.join(REPO_ROOT, "dfn_cbs")


def read_matrix(path):
    with gzip.open(path, "rt") as f:
        lines = [line.rstrip("\n").split("\t") for line in f]
    header, rows = lines[0], lines[1:]
    return header, {row[0]: row[1:] for row in rows}


def run_dfn_cbs(tmp_path, extra_args=None):
    counts_path, bins_path = gen_cbs_fixture.main(str(tmp_path))
    out_prefix = str(tmp_path / "out")
    args = [DFN_CBS, "--counts", counts_path, "--bins", bins_path,
            "--out-prefix", out_prefix, "--min-reads", "1000",
            "--perms", "50", "--min-seg-len", "5", "--threads", "2", "--quiet"]
    if extra_args:
        args.extend(extra_args)
    result = subprocess.run(args, capture_output=True, text=True)
    return result, out_prefix


def test_all_four_outputs_are_produced_with_correct_shape(tmp_path):
    result, out_prefix = run_dfn_cbs(tmp_path)
    assert result.returncode == 0, result.stderr

    for suffix in ["gc_corrected.txt.gz", "lowess_ratio.txt.gz", "segmented_lowess_ratio.txt.gz"]:
        header, rows = read_matrix(out_prefix + "." + suffix)
        assert header == ["bin", "AAAA-1", "CCCC-1"]  # GGGG-1 filtered out
        assert len(rows) == 40  # 20 bins x 2 chromosomes

    assert os.path.exists(out_prefix + ".seg")


def test_filtered_cell_absent_from_seg_file(tmp_path):
    result, out_prefix = run_dfn_cbs(tmp_path)
    assert result.returncode == 0, result.stderr

    with open(out_prefix + ".seg") as f:
        lines = f.readlines()
    header = lines[0].rstrip("\n").split("\t")
    assert header == ["ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean"]
    ids = {line.split("\t")[0] for line in lines[1:]}
    assert "GGGG-1" not in ids
    assert ids == {"AAAA-1", "CCCC-1"}


def test_whole_chromosome_shift_survives_without_centering(tmp_path):
    result, out_prefix = run_dfn_cbs(tmp_path)
    assert result.returncode == 0, result.stderr

    header, rows = read_matrix(out_prefix + ".segmented_lowess_ratio.txt.gz")
    # header includes the leading "bin" column; rows' values have already had
    # that column sliced off, so the value index is one less than the header
    # index.
    cccc_idx = header.index("CCCC-1") - 1

    chr1_values = [float(rows[f"chr1:{i*1000}:{(i+1)*1000}"][cccc_idx]) for i in range(20)]
    chr2_values = [float(rows[f"chr2:{i*1000}:{(i+1)*1000}"][cccc_idx]) for i in range(20)]

    avg_chr1 = sum(chr1_values) / len(chr1_values)
    avg_chr2 = sum(chr2_values) / len(chr2_values)

    # lowess_ratio is mean-normalized to ~1 across the *whole* cell (not
    # per-chromosome), and with an even 20/20 bin split between chr1 (raw
    # 200/bin) and chr2 (raw 500/bin) that puts chr1 at ~0.57 and chr2 at
    # ~1.43 -- a 2.5x ratio, exactly the injected 500/200 gain. If CBS's
    # centering weren't forced off, both chromosomes' segments would
    # collapse toward the same value (ratio -> 1) instead of preserving it.
    assert avg_chr1 < 0.9
    assert avg_chr2 > 1.1
    assert avg_chr2 / avg_chr1 > 2.0


def test_every_cell_filtered_out_fails_cleanly(tmp_path):
    result, _ = run_dfn_cbs(tmp_path, extra_args=["--min-reads", "999999999"])
    assert result.returncode == 1
    assert "no cells survived" in result.stderr


def test_malformed_bins_file_fails_cleanly(tmp_path):
    counts_path, bins_path = gen_cbs_fixture.main(str(tmp_path))
    with open(bins_path, "w") as f:
        f.write("chrom\tstart\tend\n")  # no gc column
        f.write("chr1\t0\t1000\n")
    out_prefix = str(tmp_path / "out")
    result = subprocess.run(
        [DFN_CBS, "--counts", counts_path, "--bins", bins_path, "--out-prefix", out_prefix],
        capture_output=True, text=True)
    assert result.returncode == 1
    assert "missing required column" in result.stderr


def test_unwritable_out_prefix_dir_fails_cleanly(tmp_path):
    counts_path, bins_path = gen_cbs_fixture.main(str(tmp_path))
    out_prefix = str(tmp_path / "nonexistent_dir" / "out")
    result = subprocess.run(
        [DFN_CBS, "--counts", counts_path, "--bins", bins_path, "--out-prefix", out_prefix],
        capture_output=True, text=True)
    assert result.returncode == 1
    assert "libc++abi" not in result.stderr
    assert "terminating due to uncaught exception" not in result.stderr


def test_worker_exception_leaves_no_seg_or_scratch_files_behind(tmp_path):
    # A zero-bin bins.tsv (header row only) is a legitimately well-formed
    # input -- load_bin_gc and load_raw_counts both accept it -- but it
    # drives gc_correct_cell's "counts and gc must be non-empty" guard
    # inside the Phase 2 worker loop for every surviving cell, which is a
    # genuine worker-thread exception rather than a Phase-1 validation
    # error. With --min-reads 0, a cell with a trivial column-sum of 0
    # still survives cell_filter and reaches the worker loop.
    bins_path = str(tmp_path / "bins.tsv")
    counts_path = str(tmp_path / "counts.raw_counts.txt.gz")
    with open(bins_path, "w") as f:
        f.write("chrom\tstart\tend\tgc\n")
    with gzip.open(counts_path, "wt") as f:
        f.write("bin\tX-1\n")

    out_prefix = str(tmp_path / "out")
    result = subprocess.run(
        [DFN_CBS, "--counts", counts_path, "--bins", bins_path,
         "--out-prefix", out_prefix, "--min-reads", "0"],
        capture_output=True, text=True)

    assert result.returncode == 1
    assert "counts and gc must be non-empty" in result.stderr

    for suffix in [".seg", ".gc_corrected.scratch.bin",
                   ".lowess_ratio.scratch.bin",
                   ".segmented_lowess_ratio.scratch.bin"]:
        assert not os.path.exists(out_prefix + suffix), \
            f"{out_prefix + suffix} should have been cleaned up on error"
