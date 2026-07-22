import gzip
import subprocess
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "fixtures"))
import gen_smoke_fixture

REPO_ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
DFN_COPY = os.path.join(REPO_ROOT, "dfn_copy")

def read_matrix(path):
    with gzip.open(path, "rt") as f:
        lines = [line.rstrip("\n").split("\t") for line in f]
    header, rows = lines[0], lines[1:]
    return header, {row[0]: row[1:] for row in rows}

def test_clean_reads_produce_expected_counts(tmp_path):
    bam_path, bins_path = gen_smoke_fixture.main(str(tmp_path))
    out_prefix = str(tmp_path / "out")

    result = subprocess.run(
        [DFN_COPY, "--bam", bam_path, "--bins", bins_path, "--out-prefix", out_prefix],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr

    header, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    barcodes = header[1:]
    assert set(barcodes) == {"AAAA-1", "CCCC-1"}

    aaaa_idx = barcodes.index("AAAA-1")
    cccc_idx = barcodes.index("CCCC-1")

    assert int(rows["chr1:0:1000"][aaaa_idx]) == 1   # clean1
    assert int(rows["chr1:1000:2000"][aaaa_idx]) == 1  # clean2
    assert int(rows["chr1:1000:2000"][cccc_idx]) == 1  # clean3
    assert int(rows["chr2:0:1000"][cccc_idx]) == 1     # clean4

def test_last_discovered_barcode_with_no_binned_reads_does_not_crash(tmp_path):
    """Reproduces the CountsMatrix/BarcodeIndex dimension-desync bug: the
    last barcode discovered while streaming (CCCC-1) has its only read land
    outside every defined bin, so CountsMatrix's barcode dimension never
    grows to cover it via increment()'s lazy grow_to() while BarcodeIndex has
    already counted it via get_or_create(). Without the fix, this desync
    makes matrix.write_gz()'s size-mismatch check throw an uncaught
    std::runtime_error out of run() (verified separately: the process
    receives SIGABRT / returncode -6 with a bare `libc++abi: terminating due
    to uncaught exception...` message, not a clean error from the program).

    The fix (`matrix.grow_to(barcode_index.size())` immediately before
    write_gz, wrapped alongside write_gz in the same try/catch used
    elsewhere in run()) turns out to fully close the desync rather than
    merely catching it: grow_to's target is exactly barcode_index.size(),
    the same value write_gz's size-check compares against, so the mismatch
    becomes unreachable and the run succeeds cleanly. This test asserts
    that end state: a clean exit with the never-binned barcode correctly
    reported with all-zero counts, in place of a crash."""
    bam_path, bins_path = gen_smoke_fixture.gen_barcode_desync_fixture(str(tmp_path))
    out_prefix = str(tmp_path / "out")

    result = subprocess.run(
        [DFN_COPY, "--bam", bam_path, "--bins", bins_path, "--out-prefix", out_prefix],
        capture_output=True, text=True,
    )

    # A crash (e.g. SIGABRT from an uncaught exception) yields a negative
    # returncode on POSIX. The fix must produce a clean, non-negative exit
    # -- never a signal-based process death.
    assert result.returncode == 0, (
        f"expected a clean exit, got returncode={result.returncode} "
        f"(negative means the process crashed/aborted); stderr={result.stderr!r}"
    )

    header, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    barcodes = header[1:]
    assert set(barcodes) == {"AAAA-1", "CCCC-1"}

    aaaa_idx = barcodes.index("AAAA-1")
    cccc_idx = barcodes.index("CCCC-1")

    # The well-binned barcode counts normally...
    assert int(rows["chr1:0:1000"][aaaa_idx]) == 1  # binned1
    # ...and the never-binned barcode is present with all-zero counts,
    # rather than being silently dropped or corrupting the matrix.
    for bin_name in rows:
        assert int(rows[bin_name][cccc_idx]) == 0

def test_unwritable_out_prefix_dir_fails_cleanly(tmp_path):
    """Reproduces the uncaught DiscordantWriter-construction-exception crash:
    DiscordantWriter's constructor (src/discordant_writer.cpp) calls sam_open()
    on `<out-prefix>.discordant.bam` and throws std::runtime_error if that
    fails, e.g. because the parent directory of --out-prefix doesn't exist.
    Before the fix, this call in main.cpp's run() wasn't wrapped in
    try/catch, so the exception propagated uncaught out of main(), aborting
    the process (SIGABRT / negative returncode) instead of following the
    same clean "print message, return 1" contract every other fallible
    operation in run() follows (BAM open, sam_index_load, BinIndex::load,
    BarcodeIndex::load_allowlist, matrix.grow_to+write_gz).

    This is an ordinary, easily reachable user error (typo'd or missing
    --out-prefix directory), not a contrived edge case."""
    bam_path, bins_path = gen_smoke_fixture.main(str(tmp_path))
    # Parent directory "nonexistent_dir" does not exist, so sam_open() for
    # "<out_prefix>.discordant.bam" fails inside DiscordantWriter's
    # constructor.
    out_prefix = str(tmp_path / "nonexistent_dir" / "out")

    result = subprocess.run(
        [DFN_COPY, "--bam", bam_path, "--bins", bins_path, "--out-prefix", out_prefix],
        capture_output=True, text=True,
    )

    # A crash (e.g. SIGABRT from an uncaught exception) yields a negative
    # returncode on POSIX. The fix must produce a clean, non-negative,
    # non-zero exit -- never a signal-based process death.
    assert result.returncode == 1, (
        f"expected a clean returncode=1, got returncode={result.returncode} "
        f"(negative means the process crashed/aborted); stderr={result.stderr!r}"
    )
    assert "libc++abi" not in result.stderr
    assert "terminating due to uncaught exception" not in result.stderr
    # htslib itself logs an "[E::hts_open_format] ..." line to stderr before
    # DiscordantWriter's constructor throws; our own clean error message
    # (from main.cpp's catch block) must also be present.
    assert "error: cannot open discordant BAM for writing" in result.stderr
