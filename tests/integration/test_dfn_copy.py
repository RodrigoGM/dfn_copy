import gzip
import subprocess
import sys
import os
import pysam
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "fixtures"))
import gen_fixtures

REPO_ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
DFN_COPY = os.path.join(REPO_ROOT, "dfn_copy")


@pytest.fixture
def fixtures(tmp_path):
    return gen_fixtures.build_all(str(tmp_path))


def run_dfn_copy(fixtures, tmp_path, extra_args=None, out_name="out"):
    out_prefix = str(tmp_path / out_name)
    cmd = [DFN_COPY, "--bam", fixtures["bam"], "--bins", fixtures["bins"],
           "--out-prefix", out_prefix]
    if extra_args:
        cmd += extra_args
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result, out_prefix


def read_matrix(path):
    with gzip.open(path, "rt") as f:
        lines = [line.rstrip("\n").split("\t") for line in f]
    header, rows = lines[0], lines[1:]
    barcodes = header[1:]
    return barcodes, {row[0]: {bc: int(v) for bc, v in zip(barcodes, row[1:])} for row in rows}


def discordant_qnames(bam_path):
    with pysam.AlignmentFile(bam_path, "rb") as f:
        return [r.query_name for r in f]


def test_clean_reads_counted_correctly(fixtures, tmp_path):
    result, out_prefix = run_dfn_copy(fixtures, tmp_path)
    assert result.returncode == 0, result.stderr

    barcodes, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    assert rows["chr1:1000:2000"]["AAAA-1"] == 1  # clean2
    assert rows["chr1:1000:2000"]["BBBB-1"] == 1  # clean3
    assert rows["chr2:0:1000"]["BBBB-1"] == 1     # clean4


def test_decoys_do_not_inflate_bin_zero(fixtures, tmp_path):
    # qcfail1/lowmapq1/secondary1/dup1/missingcb1 all sit in chr1:0:1000
    # alongside clean1, each meant to be dropped by a different filter.
    result, out_prefix = run_dfn_copy(fixtures, tmp_path)
    assert result.returncode == 0, result.stderr
    barcodes, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    assert rows["chr1:0:1000"]["AAAA-1"] == 1  # clean1 only


def test_offbin_read_not_counted(fixtures, tmp_path):
    result, out_prefix = run_dfn_copy(fixtures, tmp_path)
    barcodes, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    total = sum(v for row in rows.values() for v in row.values())
    # offbin1 (chr1:350000) has no containing bin, so it must never
    # appear anywhere -- cross-check via the discordant BAM being the
    # only other place a read could "disappear" to (it isn't there either).
    assert "offbin1" not in discordant_qnames(out_prefix + ".discordant.bam")
    assert total > 0  # sanity: matrix isn't trivially empty


def test_properpair_normal_counted_once_not_twice(fixtures, tmp_path):
    result, out_prefix = run_dfn_copy(fixtures, tmp_path)
    assert result.returncode == 0, result.stderr
    barcodes, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    assert rows["chr1:2000:3000"]["AAAA-1"] == 1  # not 2


def test_discordant_bam_contains_split_and_cross_chrom(fixtures, tmp_path):
    result, out_prefix = run_dfn_copy(fixtures, tmp_path)
    assert result.returncode == 0, result.stderr

    qnames = discordant_qnames(out_prefix + ".discordant.bam")
    assert qnames.count("split1") == 2
    assert qnames.count("supplonoSA1") == 1
    assert qnames.count("crosschrom1") == 2
    assert qnames.count("improper1") == 2
    assert "properpair_normal" not in qnames
    assert "clean1" not in qnames


def test_max_insert_size_default_trusts_proper_pair_flag(fixtures, tmp_path):
    result, out_prefix = run_dfn_copy(fixtures, tmp_path)
    assert result.returncode == 0, result.stderr
    qnames = discordant_qnames(out_prefix + ".discordant.bam")
    assert "properpair_bigisize" not in qnames
    # Trusted as a normal pair -> resolves to read1 (tie, equal 50bp
    # lengths) -> counted once in its bin.
    _, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    assert rows["chr1:7000:8000"]["AAAA-1"] == 1
    assert rows["chr1:8000:9000"]["AAAA-1"] == 0


def test_max_insert_size_override_diverts_big_isize_pair(fixtures, tmp_path):
    result, out_prefix = run_dfn_copy(fixtures, tmp_path, extra_args=["--max-insert-size", "10000"])
    assert result.returncode == 0, result.stderr
    qnames = discordant_qnames(out_prefix + ".discordant.bam")
    assert qnames.count("properpair_bigisize") == 2
    _, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    assert rows["chr1:7000:8000"]["AAAA-1"] == 0
    assert rows["chr1:8000:9000"]["AAAA-1"] == 0


def test_auto_mode_keeps_longer_mate(fixtures, tmp_path):
    result, out_prefix = run_dfn_copy(fixtures, tmp_path)
    assert result.returncode == 0, result.stderr
    _, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    # autolen1: 27bp read1 -> chr1:3000:4000, 51bp read2 -> chr1:4000:5000.
    # auto must keep the 51bp one.
    assert rows["chr1:3000:4000"]["AAAA-1"] == 0
    assert rows["chr1:4000:5000"]["AAAA-1"] == 1


def test_r1_mode_keeps_shorter_mate(fixtures, tmp_path):
    result, out_prefix = run_dfn_copy(fixtures, tmp_path, extra_args=["--single-end-counting", "r1"])
    assert result.returncode == 0, result.stderr
    _, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    # r1 mode always keeps read1 (the 27bp one), regardless of length.
    assert rows["chr1:3000:4000"]["AAAA-1"] == 1
    assert rows["chr1:4000:5000"]["AAAA-1"] == 0


def test_auto_mode_counts_singleton_when_mate_filtered(fixtures, tmp_path):
    result, out_prefix = run_dfn_copy(fixtures, tmp_path)
    assert result.returncode == 0, result.stderr
    _, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    assert rows["chr1:5000:6000"]["AAAA-1"] == 1  # singleton1's read1


def test_r1_mode_keeps_singleton_that_is_read1(fixtures, tmp_path):
    result, out_prefix = run_dfn_copy(fixtures, tmp_path, extra_args=["--single-end-counting", "r1"])
    assert result.returncode == 0, result.stderr
    _, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    assert rows["chr1:5000:6000"]["AAAA-1"] == 1


def test_r2_mode_drops_singleton_that_is_read1(fixtures, tmp_path):
    result, out_prefix = run_dfn_copy(fixtures, tmp_path, extra_args=["--single-end-counting", "r2"])
    assert result.returncode == 0, result.stderr
    _, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    # singleton1's only survivor is read1; r2 mode requires read2 -> dropped.
    assert rows["chr1:5000:6000"]["AAAA-1"] == 0


def test_flush_trigger_read_counted_independently(fixtures, tmp_path):
    result, out_prefix = run_dfn_copy(fixtures, tmp_path)
    assert result.returncode == 0, result.stderr
    _, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    assert rows["chr1:6000:7000"]["AAAA-1"] == 1


def test_barcode_allowlist_filters_unlisted_barcodes(fixtures, tmp_path):
    result, out_prefix = run_dfn_copy(fixtures, tmp_path, extra_args=["--barcodes", fixtures["allowlist"]])
    assert result.returncode == 0, result.stderr
    barcodes, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    assert barcodes == ["AAAA-1", "BBBB-1"]  # allowlist order
    assert "ZZZZ-1" not in barcodes


def test_notallowed_barcode_counted_when_no_allowlist_given(fixtures, tmp_path):
    result, out_prefix = run_dfn_copy(fixtures, tmp_path)
    assert result.returncode == 0, result.stderr
    barcodes, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    assert "ZZZZ-1" in barcodes
    assert rows["chr1:9000:10000"]["ZZZZ-1"] == 1


def test_non_string_barcode_tag_skipped_without_crash(fixtures, tmp_path):
    # badtagtype1 carries an integer-typed CB tag (present but not type 'Z'),
    # which makes bam_aux2Z return NULL. The tool must not crash on this and
    # must simply skip the read rather than count it anywhere.
    result, out_prefix = run_dfn_copy(fixtures, tmp_path)
    assert result.returncode == 0, result.stderr
    barcodes, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    assert rows["chr1:10000:11000"] == {bc: 0 for bc in barcodes}
    assert "badtagtype1" not in discordant_qnames(out_prefix + ".discordant.bam")


def test_exclude_dups_toggle(fixtures, tmp_path):
    result_on, prefix_on = run_dfn_copy(fixtures, tmp_path, out_name="on")
    result_off, prefix_off = run_dfn_copy(fixtures, tmp_path, extra_args=["--exclude-dups", "false"],
                                          out_name="off")
    assert result_on.returncode == 0, result_on.stderr
    assert result_off.returncode == 0, result_off.stderr
    _, rows_on = read_matrix(prefix_on + ".raw_counts.txt.gz")
    _, rows_off = read_matrix(prefix_off + ".raw_counts.txt.gz")
    assert rows_on["chr1:0:1000"]["AAAA-1"] == 1    # dup1 excluded
    assert rows_off["chr1:0:1000"]["AAAA-1"] == 2   # dup1 now counted too


def test_missing_bam_index_errors_clearly(fixtures, tmp_path):
    os.remove(fixtures["bam"] + ".bai")
    result, _ = run_dfn_copy(fixtures, tmp_path)
    assert result.returncode != 0
    assert "index" in result.stderr.lower()


def test_malformed_bins_file_errors_clearly(fixtures, tmp_path):
    out_prefix = str(tmp_path / "out")
    result = subprocess.run(
        [DFN_COPY, "--bam", fixtures["bam"], "--bins", fixtures["bins_malformed"],
         "--out-prefix", out_prefix],
        capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "gc" in result.stderr.lower() or "column" in result.stderr.lower()


def test_chrom_naming_mismatch_errors_clearly(fixtures, tmp_path):
    out_prefix = str(tmp_path / "out")
    result = subprocess.run(
        [DFN_COPY, "--bam", fixtures["bam"], "--bins", fixtures["bins_mismatched"],
         "--out-prefix", out_prefix],
        capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "chromosome" in result.stderr.lower()


def test_empty_allowlist_errors_clearly(fixtures, tmp_path):
    result, _ = run_dfn_copy(fixtures, tmp_path, extra_args=["--barcodes", fixtures["allowlist_empty"]])
    assert result.returncode != 0
    assert "allowlist" in result.stderr.lower() or "empty" in result.stderr.lower()
