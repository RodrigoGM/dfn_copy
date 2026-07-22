#include "cli_args.hpp"
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

std::string usage_text() {
    return
        "dfn_bin --bam <path> --bins <path> --out-prefix <prefix> [options]\n"
        "\n"
        "Required:\n"
        "  --bam <path>            Coordinate-sorted, indexed DEFND-seq gDNA BAM\n"
        "  --bins <path>           Tab-separated bin file (header columns: chrom,start,end,gc)\n"
        "  --out-prefix <prefix>   Prefix for output files\n"
        "\n"
        "Options:\n"
        "  --barcode-tag <tag>            BAM tag holding the cell barcode (default: CB)\n"
        "  --barcodes <path>              Optional barcode allowlist file (default: none)\n"
        "  --mapq <int>                   Minimum MAPQ to keep a read (default: 30)\n"
        "  --exclude-dups <true|false>    Exclude PCR/optical duplicates (default: true)\n"
        "  --primary-alignment-only <true|false>  Exclude secondary alignments (default: true).\n"
        "                                  Supplementary alignments are always diverted to the\n"
        "                                  discordant BAM regardless of this flag; it only\n"
        "                                  controls secondary-alignment exclusion.\n"
        "  --single-end-counting <auto|r1|r2>     Fragment-counting mode (default: auto)\n"
        "  --max-insert-size <int>        Absolute discordance cutoff, 0=disabled (default: 0)\n"
        "  --position <start|midpoint|end> Bin-assignment coordinate (default: end)\n"
        "  --threads <int>                htslib decompression threads (default: 1)\n"
        "  --help                         Print this message\n";
}

namespace {

bool parse_bool(const std::string& s, bool& out) {
    if (s == "true") { out = true; return true; }
    if (s == "false") { out = false; return true; }
    return false;
}

bool parse_int(const std::string& s, int& out) {
    try {
        size_t consumed = 0;
        out = std::stoi(s, &consumed);
        return consumed == s.size();
    } catch (const std::exception&) {
        return false;
    }
}

bool parse_int64(const std::string& s, int64_t& out) {
    try {
        size_t consumed = 0;
        out = std::stoll(s, &consumed);
        return consumed == s.size();
    } catch (const std::exception&) {
        return false;
    }
}

} // namespace

bool parse_args(int argc, char** argv, Args& args, bool& help_requested) {
    help_requested = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--help") {
            help_requested = true;
            return true;
        }

        auto next_value = [&](const char* flag_name) -> const char* {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "error: %s requires a value\n", flag_name);
                return nullptr;
            }
            return argv[++i];
        };

        if (arg == "--bam") {
            const char* v = next_value("--bam"); if (!v) return false;
            args.bam_path = v;
        } else if (arg == "--bins") {
            const char* v = next_value("--bins"); if (!v) return false;
            args.bins_path = v;
        } else if (arg == "--out-prefix") {
            const char* v = next_value("--out-prefix"); if (!v) return false;
            args.out_prefix = v;
        } else if (arg == "--barcode-tag") {
            const char* v = next_value("--barcode-tag"); if (!v) return false;
            args.barcode_tag = v;
        } else if (arg == "--barcodes") {
            const char* v = next_value("--barcodes"); if (!v) return false;
            args.barcodes_path = v;
        } else if (arg == "--mapq") {
            const char* v = next_value("--mapq"); if (!v) return false;
            if (!parse_int(v, args.mapq)) {
                std::fprintf(stderr, "error: --mapq requires an integer, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--exclude-dups") {
            const char* v = next_value("--exclude-dups"); if (!v) return false;
            if (!parse_bool(v, args.exclude_dups)) {
                std::fprintf(stderr, "error: --exclude-dups requires true|false, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--primary-alignment-only") {
            const char* v = next_value("--primary-alignment-only"); if (!v) return false;
            if (!parse_bool(v, args.primary_alignment_only)) {
                std::fprintf(stderr, "error: --primary-alignment-only requires true|false, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--single-end-counting") {
            const char* v = next_value("--single-end-counting"); if (!v) return false;
            std::string sv = v;
            if (sv == "auto") args.single_end_counting = SingleEndCountingMode::Auto;
            else if (sv == "r1") args.single_end_counting = SingleEndCountingMode::R1;
            else if (sv == "r2") args.single_end_counting = SingleEndCountingMode::R2;
            else {
                std::fprintf(stderr, "error: --single-end-counting requires auto|r1|r2, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--max-insert-size") {
            const char* v = next_value("--max-insert-size"); if (!v) return false;
            if (!parse_int64(v, args.max_insert_size)) {
                std::fprintf(stderr, "error: --max-insert-size requires an integer, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--position") {
            const char* v = next_value("--position"); if (!v) return false;
            std::string sv = v;
            if (sv == "start") args.position = PositionMode::Start;
            else if (sv == "midpoint") args.position = PositionMode::Midpoint;
            else if (sv == "end") args.position = PositionMode::End;
            else {
                std::fprintf(stderr, "error: --position requires start|midpoint|end, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--threads") {
            const char* v = next_value("--threads"); if (!v) return false;
            if (!parse_int(v, args.threads)) {
                std::fprintf(stderr, "error: --threads requires an integer, got '%s'\n", v);
                return false;
            }
        } else {
            std::fprintf(stderr, "error: unknown flag '%s'\n", arg.c_str());
            return false;
        }
    }

    if (args.bam_path.empty()) {
        std::fprintf(stderr, "error: --bam is required\n");
        return false;
    }
    if (args.bins_path.empty()) {
        std::fprintf(stderr, "error: --bins is required\n");
        return false;
    }
    if (args.out_prefix.empty()) {
        std::fprintf(stderr, "error: --out-prefix is required\n");
        return false;
    }
    if (args.barcode_tag.size() != 2) {
        std::fprintf(stderr, "error: --barcode-tag must be exactly 2 characters, got '%s'\n",
                      args.barcode_tag.c_str());
        return false;
    }

    return true;
}
