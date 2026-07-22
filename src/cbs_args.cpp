#include "cbs_args.hpp"
#include <cstdio>
#include <cstdlib>

std::string cbs_usage_text() {
    return
        "dfn_cbs --counts <path> --bins <path> --out-prefix <prefix> [options]\n"
        "\n"
        "Required:\n"
        "  --counts <path>         <prefix>.raw_counts.txt.gz produced by dfn_copy\n"
        "  --bins <path>           Same bins.tsv used to produce --counts\n"
        "  --out-prefix <prefix>   Prefix for the four output files\n"
        "\n"
        "Options:\n"
        "  --min-reads <int>              Minimum raw binned-count sum to keep a cell (default: 100000)\n"
        "  --threads <int>                Worker threads for per-cell parallelism (default: 1)\n"
        "  --alpha <float>                CBS split-acceptance p-value threshold (default: 0.01)\n"
        "  --perms <int>                  CBS circular-permutation count per split (default: 1000)\n"
        "  --min-seg-len <int>            CBS minimum segment length in bins (default: 25)\n"
        "  --max-depth <int>              CBS maximum recursive split depth (default: 100)\n"
        "  --cbs-method <1cp|2cp>         CBS split search mode (default: 1cp)\n"
        "  --seed <int>                   RNG seed; per-cell RNG uses seed XOR hash(barcode) (default: 1)\n"
        "  --quiet                        Suppress progress messages on stderr\n"
        "  --help                         Print this message\n";
}

namespace {

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

bool parse_uint64(const std::string& s, uint64_t& out) {
    try {
        size_t consumed = 0;
        out = std::stoull(s, &consumed);
        return consumed == s.size();
    } catch (const std::exception&) {
        return false;
    }
}

bool parse_double(const std::string& s, double& out) {
    try {
        size_t consumed = 0;
        out = std::stod(s, &consumed);
        return consumed == s.size();
    } catch (const std::exception&) {
        return false;
    }
}

} // namespace

bool parse_cbs_args(int argc, char** argv, CbsArgs& args, bool& help_requested) {
    help_requested = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--help") {
            help_requested = true;
            return true;
        }
        if (arg == "--quiet") {
            args.quiet = true;
            continue;
        }

        auto next_value = [&](const char* flag_name) -> const char* {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "error: %s requires a value\n", flag_name);
                return nullptr;
            }
            return argv[++i];
        };

        if (arg == "--counts") {
            const char* v = next_value("--counts"); if (!v) return false;
            args.counts_path = v;
        } else if (arg == "--bins") {
            const char* v = next_value("--bins"); if (!v) return false;
            args.bins_path = v;
        } else if (arg == "--out-prefix") {
            const char* v = next_value("--out-prefix"); if (!v) return false;
            args.out_prefix = v;
        } else if (arg == "--min-reads") {
            const char* v = next_value("--min-reads"); if (!v) return false;
            if (!parse_int64(v, args.min_reads)) {
                std::fprintf(stderr, "error: --min-reads requires an integer, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--threads") {
            const char* v = next_value("--threads"); if (!v) return false;
            if (!parse_int(v, args.threads)) {
                std::fprintf(stderr, "error: --threads requires an integer, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--alpha") {
            const char* v = next_value("--alpha"); if (!v) return false;
            if (!parse_double(v, args.alpha)) {
                std::fprintf(stderr, "error: --alpha requires a number, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--perms") {
            const char* v = next_value("--perms"); if (!v) return false;
            if (!parse_int(v, args.perms)) {
                std::fprintf(stderr, "error: --perms requires an integer, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--min-seg-len") {
            const char* v = next_value("--min-seg-len"); if (!v) return false;
            if (!parse_int(v, args.min_seg_len)) {
                std::fprintf(stderr, "error: --min-seg-len requires an integer, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--max-depth") {
            const char* v = next_value("--max-depth"); if (!v) return false;
            if (!parse_int(v, args.max_depth)) {
                std::fprintf(stderr, "error: --max-depth requires an integer, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--cbs-method") {
            const char* v = next_value("--cbs-method"); if (!v) return false;
            std::string sv = v;
            if (sv == "1cp") args.cbs_method = CbsMethod::OneCp;
            else if (sv == "2cp") args.cbs_method = CbsMethod::TwoCp;
            else {
                std::fprintf(stderr, "error: --cbs-method requires 1cp|2cp, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--seed") {
            const char* v = next_value("--seed"); if (!v) return false;
            if (!parse_uint64(v, args.seed)) {
                std::fprintf(stderr, "error: --seed requires an integer, got '%s'\n", v);
                return false;
            }
        } else {
            std::fprintf(stderr, "error: unknown flag '%s'\n", arg.c_str());
            return false;
        }
    }

    if (args.counts_path.empty()) {
        std::fprintf(stderr, "error: --counts is required\n");
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

    return true;
}
