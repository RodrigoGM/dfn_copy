#!/usr/bin/env python3
"""GC-bias correction for dfn_copy's raw counts matrix.

Method: for each cell (column), fit a LOWESS curve of raw counts against
each bin's GC fraction (statsmodels.nonparametric.lowess), then rescale
the fitted values so the corrected column's median equals that cell's own
genome-wide median raw count. This is the QDNAseq/Ginkgo-style approach
(correct each bin toward the genome-wide median), applied per cell -- the
single-cell CNV norm, at the cost of one LOWESS fit per cell rather than
one fit pooled across all cells (see design spec, Component 2).
"""
import argparse
import gzip
import sys

import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess


def read_counts_matrix(path):
    with gzip.open(path, "rt") as f:
        lines = [line.rstrip("\n").split("\t") for line in f]
    header = lines[0]
    barcodes = header[1:]
    bins = []
    matrix = []
    for row in lines[1:]:
        bins.append(row[0])
        matrix.append([float(v) for v in row[1:]])
    return bins, barcodes, np.array(matrix)


def read_bin_gc(path, bins):
    gc_by_bin = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        i_chrom = header.index("chrom")
        i_start = header.index("start")
        i_end = header.index("end")
        i_gc = header.index("gc")
        for line in f:
            fields = line.rstrip("\n").split("\t")
            name = f"{fields[i_chrom]}:{fields[i_start]}:{fields[i_end]}"
            gc_by_bin[name] = float(fields[i_gc])
    return np.array([gc_by_bin[b] for b in bins])


def correct_gc(matrix, gc):
    """matrix: (n_bins, n_cells) raw counts. gc: (n_bins,) GC fraction.
    Returns a matrix of the same shape, GC-corrected per cell."""
    n_bins, n_cells = matrix.shape
    corrected = np.zeros_like(matrix)
    for cell in range(n_cells):
        counts = matrix[:, cell]
        median = np.median(counts)
        if median == 0:
            corrected[:, cell] = counts
            continue
        fitted = lowess(counts, gc, return_sorted=False)
        fitted = np.where(fitted <= 0, np.finfo(float).eps, fitted)
        corrected[:, cell] = counts * (median / fitted)
    return corrected


def main():
    parser = argparse.ArgumentParser(description="GC-correct a dfn_copy raw counts matrix")
    parser.add_argument("--counts", required=True, help="Path to <prefix>.raw_counts.txt.gz")
    parser.add_argument("--bins", required=True, help="Path to the bins.tsv used by dfn_copy")
    parser.add_argument("--out-prefix", required=True, help="Prefix for the output file")
    args = parser.parse_args()

    bins, barcodes, matrix = read_counts_matrix(args.counts)
    gc = read_bin_gc(args.bins, bins)
    corrected = correct_gc(matrix, gc)

    out_path = args.out_prefix + ".gc_corrected.txt.gz"
    with gzip.open(out_path, "wt") as f:
        f.write("bin\t" + "\t".join(barcodes) + "\n")
        for bin_name, row in zip(bins, corrected):
            f.write(bin_name + "\t" + "\t".join(f"{v:.4f}" for v in row) + "\n")


if __name__ == "__main__":
    main()
