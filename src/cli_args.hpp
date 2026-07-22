#pragma once
#include <cstdint>
#include <string>

enum class PositionMode { Start, Midpoint, End };
enum class SingleEndCountingMode { Auto, R1, R2 };

struct Args {
    std::string bam_path;
    std::string bins_path;
    std::string out_prefix;
    std::string barcode_tag = "CB";
    std::string barcodes_path;
    int mapq = 30;
    bool exclude_dups = true;
    bool primary_alignment_only = true;
    SingleEndCountingMode single_end_counting = SingleEndCountingMode::Auto;
    int64_t max_insert_size = 0;
    PositionMode position = PositionMode::End;
    int threads = 1;
};

// Parses argv into args. Returns true on success (including a bare
// --help, which sets help_requested and leaves other fields unset).
// On error, prints a specific message to stderr and returns false.
bool parse_args(int argc, char** argv, Args& args, bool& help_requested);

std::string usage_text();
