#pragma once
#include <cstdint>
#include <string>

enum class CbsMethod { OneCp, TwoCp };

struct CbsArgs {
    std::string counts_path;
    std::string bins_path;
    std::string out_prefix;
    int64_t min_reads = 100000;
    int threads = 1;
    double alpha = 0.01;
    int perms = 1000;
    int min_seg_len = 25;
    int max_depth = 100;
    CbsMethod cbs_method = CbsMethod::OneCp;
    uint64_t seed = 1;
    bool quiet = false;
};

// Parses argv into args. Returns true on success (including a bare
// --help, which sets help_requested and leaves other fields unset).
// On error, prints a specific message to stderr and returns false.
bool parse_cbs_args(int argc, char** argv, CbsArgs& args, bool& help_requested);

std::string cbs_usage_text();
