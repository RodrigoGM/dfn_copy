#!/usr/bin/env python3
"""Generates a tiny BAM + bins.tsv covering only the 'clean read' path,
for Task 9's end-to-end smoke test. Comprehensive edge-case fixtures
(discordant reads, fragment-pairing modes, error paths) are in
gen_fixtures.py, added in Task 10."""
import pysam
import sys

def main(out_dir):
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 10000}, {"SN": "chr2", "LN": 10000}],
    }

    bam_path = f"{out_dir}/smoke.bam"
    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        def add_read(qname, tid, pos, cb, mapq=60, seq_len=50):
            a = pysam.AlignedSegment()
            a.query_name = qname
            a.query_sequence = "A" * seq_len
            a.flag = 0
            a.reference_id = tid
            a.reference_start = pos
            a.mapping_quality = mapq
            a.cigarstring = f"{seq_len}M"
            a.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
            a.tags = [("CB", cb)]
            out.write(a)

        add_read("clean1", 0, 150, "AAAA-1")
        add_read("clean2", 0, 1200, "AAAA-1")
        add_read("clean3", 0, 1300, "CCCC-1")
        add_read("clean4", 1, 150, "CCCC-1")

    pysam.index(bam_path)

    bins_path = f"{out_dir}/smoke_bins.tsv"
    with open(bins_path, "w") as f:
        f.write("chrom\tstart\tend\tgc\n")
        f.write("chr1\t0\t1000\t0.40\n")
        f.write("chr1\t1000\t2000\t0.50\n")
        f.write("chr2\t0\t1000\t0.45\n")

    return bam_path, bins_path

def gen_barcode_desync_fixture(out_dir):
    """Generates a tiny BAM + bins.tsv that reproduces the barcode/matrix
    dimension desync bug: the LAST barcode discovered while streaming
    (BarcodeIndex::get_or_create() assigns it an index the moment its first
    read is seen) has its only read land at a genomic position with no
    defined bin, so CountsMatrix::increment() (and therefore its lazy
    grow_to()) is never called for that barcode's index. CountsMatrix ends
    up one barcode short of BarcodeIndex, and matrix.write_gz()'s
    size-mismatch check should fail cleanly rather than let an uncaught
    exception abort the process."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 10000}, {"SN": "chr2", "LN": 10000}],
    }

    bam_path = f"{out_dir}/desync.bam"
    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        def add_read(qname, tid, pos, cb, mapq=60, seq_len=50):
            a = pysam.AlignedSegment()
            a.query_name = qname
            a.query_sequence = "A" * seq_len
            a.flag = 0
            a.reference_id = tid
            a.reference_start = pos
            a.mapping_quality = mapq
            a.cigarstring = f"{seq_len}M"
            a.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
            a.tags = [("CB", cb)]
            out.write(a)

        # First-discovered barcode: normal read, lands in a defined bin.
        add_read("binned1", 0, 150, "AAAA-1")
        # Last-discovered barcode: its only read sits at chr1:5000, past
        # every bin range in bins.tsv below (bins only cover chr1:0-2000
        # and chr2:0-1000), so it never lands in a bin.
        add_read("unbinned1", 0, 5000, "CCCC-1")

    pysam.index(bam_path)

    bins_path = f"{out_dir}/desync_bins.tsv"
    with open(bins_path, "w") as f:
        f.write("chrom\tstart\tend\tgc\n")
        f.write("chr1\t0\t1000\t0.40\n")
        f.write("chr1\t1000\t2000\t0.50\n")
        f.write("chr2\t0\t1000\t0.45\n")

    return bam_path, bins_path

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
