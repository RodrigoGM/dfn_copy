import gzip
import subprocess
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "fixtures"))
import gen_smoke_fixture

REPO_ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
DFN_BIN = os.path.join(REPO_ROOT, "dfn_bin")

def read_matrix(path):
    with gzip.open(path, "rt") as f:
        lines = [line.rstrip("\n").split("\t") for line in f]
    header, rows = lines[0], lines[1:]
    return header, {row[0]: row[1:] for row in rows}

def test_clean_reads_produce_expected_counts(tmp_path):
    bam_path, bins_path = gen_smoke_fixture.main(str(tmp_path))
    out_prefix = str(tmp_path / "out")

    result = subprocess.run(
        [DFN_BIN, "--bam", bam_path, "--bins", bins_path, "--out-prefix", out_prefix],
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
