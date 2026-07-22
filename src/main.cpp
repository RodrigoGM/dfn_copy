#include "cli_args.hpp"
#include "bins.hpp"
#include "barcode_index.hpp"
#include "counts_matrix.hpp"
#include "read_filter.hpp"
#include "fragment_pairing.hpp"
#include "discordant_writer.hpp"

#include <htslib/sam.h>
#include <cstdio>
#include <memory>
#include <vector>
#include <string>

namespace {

void bin_candidate(const FragmentCandidate& c, const BinIndex& bin_index,
                    const std::vector<std::string>& tid_to_chrom,
                    const Args& args, CountsMatrix& matrix) {
    int64_t coord = compute_bin_coordinate(c.align_start, c.align_end, args.position);
    auto bin_idx = bin_index.find(tid_to_chrom[c.tid], coord);
    if (bin_idx) {
        matrix.increment(*bin_idx, c.barcode_idx);
    }
}

int run(const Args& args) {
    // --- Validate inputs up front ---
    samFile* bam_fp = sam_open(args.bam_path.c_str(), "r");
    if (!bam_fp) {
        std::fprintf(stderr, "error: cannot open BAM file: %s\n", args.bam_path.c_str());
        return 1;
    }
    sam_hdr_t* hdr = sam_hdr_read(bam_fp);
    if (!hdr) {
        std::fprintf(stderr, "error: cannot read BAM header: %s\n", args.bam_path.c_str());
        return 1;
    }

    hts_idx_t* idx = sam_index_load(bam_fp, args.bam_path.c_str());
    if (!idx) {
        std::fprintf(stderr,
            "error: no BAM index found for %s (run `samtools index %s`)\n",
            args.bam_path.c_str(), args.bam_path.c_str());
        return 1;
    }
    hts_idx_destroy(idx);

    BinIndex bin_index;
    try {
        bin_index = BinIndex::load(args.bins_path);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }

    std::vector<std::string> bam_contigs;
    std::vector<std::string> tid_to_chrom(sam_hdr_nref(hdr));
    for (int i = 0; i < sam_hdr_nref(hdr); ++i) {
        std::string name = sam_hdr_tid2name(hdr, i);
        bam_contigs.push_back(name);
        tid_to_chrom[i] = name;
    }
    std::vector<std::string> bin_chroms;
    for (const auto& b : bin_index.bins()) bin_chroms.push_back(b.chrom);
    std::string naming_error = validate_chrom_naming(bam_contigs, bin_chroms);
    if (!naming_error.empty()) {
        std::fprintf(stderr, "error: %s\n", naming_error.c_str());
        return 1;
    }

    BarcodeIndex barcode_index;
    if (!args.barcodes_path.empty()) {
        try {
            barcode_index = BarcodeIndex::load_allowlist(args.barcodes_path);
        } catch (const std::exception& e) {
            std::fprintf(stderr, "error: %s\n", e.what());
            return 1;
        }
    }

    // --- Set up outputs ---
    std::string discordant_path = args.out_prefix + ".discordant.bam";
    std::unique_ptr<DiscordantWriter> discordant_writer;
    try {
        discordant_writer = std::make_unique<DiscordantWriter>(discordant_path, hdr);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }

    CountsMatrix matrix(bin_index.size(), barcode_index.size());

    if (args.threads > 1) {
        hts_set_threads(bam_fp, args.threads);
    }

    // --- Stream the BAM once ---
    FragmentPairingBuffer pairing_buffer(args.single_end_counting);
    bam1_t* rec = bam_init1();

    while (sam_read1(bam_fp, hdr, rec) >= 0) {
        if (is_unconditionally_dropped(rec)) continue;

        DiscordantReason reason = classify_discordant(rec, args);
        if (reason != DiscordantReason::None) {
            try {
                discordant_writer->write(rec);
            } catch (const std::exception& e) {
                std::fprintf(stderr, "error: failed writing to discordant BAM: %s\n", e.what());
                return 1;
            }
            for (const auto& flushed : pairing_buffer.flush_up_to(rec->core.tid, rec->core.pos)) {
                bin_candidate(flushed, bin_index, tid_to_chrom, args, matrix);
            }
            continue;
        }

        for (const auto& flushed : pairing_buffer.flush_up_to(rec->core.tid, rec->core.pos)) {
            bin_candidate(flushed, bin_index, tid_to_chrom, args, matrix);
        }

        if (!passes_quality_filters(rec, args)) continue;

        uint8_t* bc_aux = bam_aux_get(rec, args.barcode_tag.c_str());
        if (!bc_aux) continue;
        char* bc_str = bam_aux2Z(bc_aux);
        if (!bc_str) continue;
        std::string barcode = bc_str;

        auto barcode_idx = barcode_index.get_or_create(barcode);
        if (!barcode_idx) continue;

        int64_t align_start = rec->core.pos;
        int64_t align_end = bam_endpos(rec);

        bool treat_as_unpaired =
            !(rec->core.flag & BAM_FPAIRED) || (rec->core.flag & BAM_FMUNMAP);

        if (treat_as_unpaired) {
            int64_t coord = compute_bin_coordinate(align_start, align_end, args.position);
            auto bin_idx = bin_index.find(tid_to_chrom[rec->core.tid], coord);
            if (bin_idx) matrix.increment(*bin_idx, *barcode_idx);
            continue;
        }

        FragmentCandidate candidate;
        candidate.qname = bam_get_qname(rec);
        candidate.tid = rec->core.tid;
        candidate.align_start = align_start;
        candidate.align_end = align_end;
        candidate.mate_pos = rec->core.mpos;
        candidate.barcode_idx = *barcode_idx;
        candidate.is_read1 = (rec->core.flag & BAM_FREAD1) != 0;

        auto resolved = pairing_buffer.process(candidate);
        if (resolved) {
            bin_candidate(*resolved, bin_index, tid_to_chrom, args, matrix);
        }
    }

    for (const auto& flushed : pairing_buffer.flush_all()) {
        bin_candidate(flushed, bin_index, tid_to_chrom, args, matrix);
    }

    bam_destroy1(rec);
    try {
        discordant_writer->close();
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: failed closing discordant BAM: %s\n", e.what());
        return 1;
    }
    sam_hdr_destroy(hdr);
    sam_close(bam_fp);

    // --- Write output matrix ---
    std::vector<std::string> bin_names;
    for (size_t i = 0; i < bin_index.size(); ++i) bin_names.push_back(bin_index.name(i));

    try {
        matrix.grow_to(barcode_index.size());
        matrix.write_gz(args.out_prefix + ".raw_counts.txt.gz", bin_names, barcode_index.names());
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }

    return 0;
}

} // namespace

int main(int argc, char** argv) {
    Args args;
    bool help_requested = false;

    if (!parse_args(argc, argv, args, help_requested)) {
        return 1;
    }
    if (help_requested) {
        std::fputs(usage_text().c_str(), stdout);
        return 0;
    }

    return run(args);
}
