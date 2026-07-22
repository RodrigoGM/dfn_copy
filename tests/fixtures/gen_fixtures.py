#!/usr/bin/env python3
"""Builds the full fixture set for Task 10's integration tests: one BAM
covering every filtering/discordant/pairing scenario, plus supporting
bins/allowlist/malformed files for the error-handling tests.

Records are listed below in the exact order they're written, which must be
true ascending (tid, pos) order (unmapped last) to match a real
coordinate-sorted BAM -- the tool's fragment-pairing buffer assumes this.
"""
import pysam

CHR1_LEN = 400000
CHR2_LEN = 400000

# SAM FLAG bits used below (see design spec's CLI contract / discordant section):
PAIRED = 1
PROPER_PAIR = 2
UNMAP = 4
MUNMAP = 8
REVERSE = 16
MREVERSE = 32
READ1 = 64
READ2 = 128
SECONDARY = 256
QCFAIL = 512
DUP = 1024
SUPPLEMENTARY = 2048


def _add(out, qname, flag, tid, pos, mapq, mtid, mpos, isize, cb=None,
         seq_len=50, sa=None):
    a = pysam.AlignedSegment()
    a.query_name = qname
    a.query_sequence = "A" * seq_len
    a.flag = flag
    a.reference_id = tid
    a.reference_start = pos
    a.mapping_quality = mapq
    a.cigarstring = f"{seq_len}M"
    a.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
    a.next_reference_id = mtid
    a.next_reference_start = mpos
    a.template_length = isize
    tags = []
    if cb is not None:
        tags.append(("CB", cb))
    if sa is not None:
        tags.append(("SA", sa))
    a.tags = tags
    out.write(a)


def build_main_bam(out_dir):
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": CHR1_LEN}, {"SN": "chr2", "LN": CHR2_LEN}],
    }
    bam_path = f"{out_dir}/fixtures.bam"

    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        # ===== chr1 (tid=0), strictly ascending pos =====

        # --- bin chr1:0-1000: clean1 plus decoys that must NOT inflate it ---
        _add(out, "clean1", 0, 0, 150, 60, -1, -1, 0, cb="AAAA-1")
        _add(out, "qcfail1", QCFAIL, 0, 160, 60, -1, -1, 0, cb="AAAA-1")
        _add(out, "lowmapq1", 0, 0, 170, 5, -1, -1, 0, cb="AAAA-1")
        _add(out, "secondary1", SECONDARY, 0, 180, 60, -1, -1, 0, cb="AAAA-1")
        _add(out, "dup1", DUP, 0, 190, 60, -1, -1, 0, cb="AAAA-1")
        _add(out, "missingcb1", 0, 0, 200, 60, -1, -1, 0, cb=None)

        # --- bin chr1:1000-2000: two different barcodes ---
        _add(out, "clean2", 0, 0, 1200, 60, -1, -1, 0, cb="AAAA-1")
        _add(out, "clean3", 0, 0, 1300, 60, -1, -1, 0, cb="BBBB-1")

        # --- bin chr1:2000-3000: proper pair, small insert -> counted ONCE ---
        _add(out, "properpair_normal", PAIRED | PROPER_PAIR | MREVERSE | READ1,
             0, 2400, 60, 0, 2600, 250, cb="AAAA-1")
        _add(out, "properpair_normal", PAIRED | PROPER_PAIR | REVERSE | READ2,
             0, 2600, 60, 0, 2400, -250, cb="AAAA-1")

        # --- autolen1: 27bp read1 in bin chr1:3000-4000, 51bp read2 in bin
        # chr1:4000-5000 -- separate bins so "which mate won" is directly
        # observable by which bin gets the count. ---
        _add(out, "autolen1", PAIRED | PROPER_PAIR | MREVERSE | READ1,
             0, 3100, 60, 0, 4100, 1050, cb="AAAA-1", seq_len=27)

        # --- improper1: interleaved here, diverted regardless of position ---
        _add(out, "improper1", PAIRED | READ1, 0, 3500, 60, 0, 3600, 150, cb="AAAA-1")
        _add(out, "improper1", PAIRED | READ2, 0, 3600, 60, 0, 3500, -150, cb="AAAA-1")

        _add(out, "autolen1", PAIRED | PROPER_PAIR | REVERSE | READ2,
             0, 4100, 60, 0, 3100, -1050, cb="AAAA-1", seq_len=51)

        # --- bin chr1:5000-6000: singleton1 (read2 filtered by low MAPQ
        # before ever reaching the pairing buffer) ---
        _add(out, "singleton1", PAIRED | PROPER_PAIR | MREVERSE | READ1,
             0, 5100, 60, 0, 5200, 150, cb="AAAA-1")
        _add(out, "singleton1", PAIRED | PROPER_PAIR | REVERSE | READ2,
             0, 5200, 5, 0, 5100, -150, cb="AAAA-1")

        # --- bin chr1:6000-7000: unpaired read whose arrival forces the
        # pairing buffer to flush singleton1 (mate_pos=5200 < 6100) ---
        _add(out, "flush_trigger1", 0, 0, 6100, 60, -1, -1, 0, cb="AAAA-1")

        # --- properpair_bigisize: proper-pair flag set, huge isize (only
        # the isize tag is huge -- actual positions stay in separate,
        # observable bins chr1:7000-8000 / chr1:8000-9000). Default run
        # trusts the proper-pair flag; --max-insert-size overrides it. ---
        _add(out, "properpair_bigisize", PAIRED | PROPER_PAIR | MREVERSE | READ1,
             0, 7100, 60, 0, 8100, 290000, cb="AAAA-1")

        # --- crosschrom1's chr1-side mate, diverted regardless of position ---
        _add(out, "crosschrom1", PAIRED | MREVERSE | READ1, 0, 7500, 60, 1, 500, 0, cb="AAAA-1")

        _add(out, "properpair_bigisize", PAIRED | PROPER_PAIR | REVERSE | READ2,
             0, 8100, 60, 0, 7100, -290000, cb="AAAA-1")

        # --- bin chr1:9000-10000: barcode not on the allowlist (when one is used) ---
        _add(out, "notallowed1", 0, 0, 9100, 60, -1, -1, 0, cb="ZZZZ-1")

        # --- split/chimeric, diverted regardless of position ---
        _add(out, "supplonoSA1", SUPPLEMENTARY, 0, 9500, 60, -1, -1, 0, cb="AAAA-1")
        _add(out, "split1", 0, 0, 9600, 60, -1, -1, 0, cb="AAAA-1",
             sa="chr2,601,+,50M,60,0;")

        # --- off-bin: far past every defined bin ---
        _add(out, "offbin1", 0, 0, 350000, 60, -1, -1, 0, cb="AAAA-1")

        # ===== chr2 (tid=1), strictly ascending pos =====
        _add(out, "clean4", 0, 1, 150, 60, -1, -1, 0, cb="BBBB-1")
        _add(out, "crosschrom1", PAIRED | REVERSE | READ2, 1, 500, 60, 0, 7500, 0, cb="AAAA-1")
        _add(out, "split1", SUPPLEMENTARY, 1, 600, 60, -1, -1, 0, cb="AAAA-1")

        # ===== unmapped, sorted last regardless of chromosome =====
        _add(out, "unmapped1", UNMAP, -1, -1, 0, -1, -1, 0, cb="AAAA-1")

    pysam.index(bam_path)
    return bam_path


def build_bins_file(out_dir):
    path = f"{out_dir}/bins.tsv"
    with open(path, "w") as f:
        f.write("chrom\tstart\tend\tgc\n")
        f.write("chr1\t0\t1000\t0.40\n")
        f.write("chr1\t1000\t2000\t0.50\n")
        f.write("chr1\t2000\t3000\t0.45\n")
        f.write("chr1\t3000\t4000\t0.50\n")
        f.write("chr1\t4000\t5000\t0.55\n")
        f.write("chr1\t5000\t6000\t0.45\n")
        f.write("chr1\t6000\t7000\t0.50\n")
        f.write("chr1\t7000\t8000\t0.45\n")
        f.write("chr1\t8000\t9000\t0.50\n")
        f.write("chr1\t9000\t10000\t0.55\n")
        f.write("chr2\t0\t1000\t0.45\n")
    return path


def build_mismatched_bins_file(out_dir):
    # Uses "1" instead of "chr1" to trigger the naming-mismatch error.
    path = f"{out_dir}/bins_mismatched.tsv"
    with open(path, "w") as f:
        f.write("chrom\tstart\tend\tgc\n")
        f.write("1\t0\t1000\t0.40\n")
    return path


def build_malformed_bins_file(out_dir):
    path = f"{out_dir}/bins_malformed.tsv"
    with open(path, "w") as f:
        f.write("chrom\tstart\tend\n")  # missing 'gc' column
        f.write("chr1\t0\t1000\n")
    return path


def build_allowlist_file(out_dir):
    path = f"{out_dir}/allowlist.txt"
    with open(path, "w") as f:
        f.write("AAAA-1\n")
        f.write("BBBB-1\n")
    return path


def build_empty_allowlist_file(out_dir):
    path = f"{out_dir}/allowlist_empty.txt"
    open(path, "w").close()
    return path


def build_all(out_dir):
    return {
        "bam": build_main_bam(out_dir),
        "bins": build_bins_file(out_dir),
        "bins_mismatched": build_mismatched_bins_file(out_dir),
        "bins_malformed": build_malformed_bins_file(out_dir),
        "allowlist": build_allowlist_file(out_dir),
        "allowlist_empty": build_empty_allowlist_file(out_dir),
    }
