#include "test_framework.hpp"
#include "cli_args.hpp"
#include <vector>
#include <cstring>

std::pair<int, char**> make_argv(std::vector<std::string> args) {
    static std::vector<std::string> storage;
    static std::vector<char*> ptrs;
    storage = args;
    ptrs.clear();
    for (auto& s : storage) ptrs.push_back(const_cast<char*>(s.c_str()));
    return {static_cast<int>(ptrs.size()), ptrs.data()};
}

void test_defaults_with_required_only() {
    auto [argc, argv] = make_argv({"dfn_copy", "--bam", "in.bam", "--bins", "bins.tsv",
                                    "--out-prefix", "out"});
    Args args;
    bool help = false;
    ASSERT_TRUE(parse_args(argc, argv, args, help));
    ASSERT_TRUE(!help);
    ASSERT_EQ(args.bam_path, std::string("in.bam"));
    ASSERT_EQ(args.bins_path, std::string("bins.tsv"));
    ASSERT_EQ(args.out_prefix, std::string("out"));
    ASSERT_EQ(args.barcode_tag, std::string("CB"));
    ASSERT_EQ(args.barcodes_path, std::string(""));
    ASSERT_EQ(args.mapq, 30);
    ASSERT_TRUE(args.exclude_dups);
    ASSERT_TRUE(args.primary_alignment_only);
    ASSERT_TRUE(args.single_end_counting == SingleEndCountingMode::Auto);
    ASSERT_EQ(args.max_insert_size, static_cast<int64_t>(0));
    ASSERT_TRUE(args.position == PositionMode::End);
    ASSERT_EQ(args.threads, 1);
}

void test_overrides_all_flags() {
    auto [argc, argv] = make_argv({
        "dfn_copy", "--bam", "in.bam", "--bins", "bins.tsv", "--out-prefix", "out",
        "--barcode-tag", "CR", "--barcodes", "allow.txt", "--mapq", "20",
        "--exclude-dups", "false", "--primary-alignment-only", "false",
        "--single-end-counting", "r2", "--max-insert-size", "5000",
        "--position", "start", "--threads", "4"});
    Args args;
    bool help = false;
    ASSERT_TRUE(parse_args(argc, argv, args, help));
    ASSERT_EQ(args.barcode_tag, std::string("CR"));
    ASSERT_EQ(args.barcodes_path, std::string("allow.txt"));
    ASSERT_EQ(args.mapq, 20);
    ASSERT_TRUE(!args.exclude_dups);
    ASSERT_TRUE(!args.primary_alignment_only);
    ASSERT_TRUE(args.single_end_counting == SingleEndCountingMode::R2);
    ASSERT_EQ(args.max_insert_size, static_cast<int64_t>(5000));
    ASSERT_TRUE(args.position == PositionMode::Start);
    ASSERT_EQ(args.threads, 4);
}

void test_help_flag() {
    auto [argc, argv] = make_argv({"dfn_copy", "--help"});
    Args args;
    bool help = false;
    ASSERT_TRUE(parse_args(argc, argv, args, help));
    ASSERT_TRUE(help);
}

void test_missing_required_fails() {
    auto [argc, argv] = make_argv({"dfn_copy", "--bins", "bins.tsv", "--out-prefix", "out"});
    Args args;
    bool help = false;
    ASSERT_TRUE(!parse_args(argc, argv, args, help));
}

void test_unknown_flag_fails() {
    auto [argc, argv] = make_argv({"dfn_copy", "--bam", "in.bam", "--bins", "bins.tsv",
                                    "--out-prefix", "out", "--nonsense", "1"});
    Args args;
    bool help = false;
    ASSERT_TRUE(!parse_args(argc, argv, args, help));
}

void test_invalid_enum_value_fails() {
    auto [argc, argv] = make_argv({"dfn_copy", "--bam", "in.bam", "--bins", "bins.tsv",
                                    "--out-prefix", "out", "--position", "middle"});
    Args args;
    bool help = false;
    ASSERT_TRUE(!parse_args(argc, argv, args, help));
}

void test_barcode_tag_must_be_two_chars() {
    auto [argc, argv] = make_argv({"dfn_copy", "--bam", "in.bam", "--bins", "bins.tsv",
                                    "--out-prefix", "out", "--barcode-tag", "TOOLONG"});
    Args args;
    bool help = false;
    ASSERT_TRUE(!parse_args(argc, argv, args, help));
}

int main() {
    test_defaults_with_required_only();
    test_overrides_all_flags();
    test_help_flag();
    test_missing_required_fails();
    test_unknown_flag_fails();
    test_invalid_enum_value_fails();
    test_barcode_tag_must_be_two_chars();
    TEST_REPORT();
}
