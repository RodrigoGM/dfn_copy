"""Generates a small bins.tsv + raw_counts.txt.gz pair for dfn_cbs
integration tests -- no BAM/pysam involved, since dfn_cbs's input is
already a matrix."""
import gzip
import os


def main(out_dir):
    bins_path = os.path.join(out_dir, "bins.tsv")
    counts_path = os.path.join(out_dir, "counts.raw_counts.txt.gz")

    bins = []
    for i in range(20):
        gc = 0.4 + 0.01 * (i % 5)
        bins.append(("chr1", i * 1000, (i + 1) * 1000, gc))
    for i in range(20):
        gc = 0.4 + 0.01 * (i % 5)
        bins.append(("chr2", i * 1000, (i + 1) * 1000, gc))

    with open(bins_path, "w") as f:
        f.write("chrom\tstart\tend\tgc\n")
        for chrom, start, end, gc in bins:
            f.write(f"{chrom}\t{start}\t{end}\t{gc}\n")

    # AAAA-1: diploid, ~200 reads/bin on both chromosomes, sum=8000 -- survives --min-reads 1000.
    # CCCC-1: ~200 reads/bin on chr1, ~500 reads/bin on chr2 (a whole-chromosome
    #         gain, ~2.5x) -- also survives, and its segmented output must
    #         preserve the chr1 vs chr2 difference (no centering).
    # GGGG-1: 1 read/bin everywhere, sum=40 -- filtered out by --min-reads 1000,
    #         must be absent from every output file.
    barcodes = ["AAAA-1", "CCCC-1", "GGGG-1"]
    counts = {"AAAA-1": [], "CCCC-1": [], "GGGG-1": []}
    for chrom, start, end, gc in bins:
        counts["AAAA-1"].append(200)
        counts["CCCC-1"].append(200 if chrom == "chr1" else 500)
        counts["GGGG-1"].append(1)

    with gzip.open(counts_path, "wt") as f:
        f.write("bin\t" + "\t".join(barcodes) + "\n")
        for i, (chrom, start, end, gc) in enumerate(bins):
            row = [f"{chrom}:{start}:{end}"] + [str(counts[b][i]) for b in barcodes]
            f.write("\t".join(row) + "\n")

    return counts_path, bins_path


if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
