#include "test_framework.hpp"
#include "cbs_args.hpp"
#include <vector>

std::pair<int, char**> make_argv(std::vector<std::string> args) {
    static std::vector<std::string> storage;
    static std::vector<char*> ptrs;
    storage = args;
    ptrs.clear();
    for (auto& s : storage) ptrs.push_back(const_cast<char*>(s.c_str()));
    return {static_cast<int>(ptrs.size()), ptrs.data()};
}

void test_defaults_with_required_only() {
    auto [argc, argv] = make_argv({"dfn_cbs", "--counts", "in.txt.gz",
                                    "--bins", "bins.tsv", "--out-prefix", "out"});
    CbsArgs args;
    bool help = false;
    ASSERT_TRUE(parse_cbs_args(argc, argv, args, help));
    ASSERT_TRUE(!help);
    ASSERT_EQ(args.counts_path, std::string("in.txt.gz"));
    ASSERT_EQ(args.bins_path, std::string("bins.tsv"));
    ASSERT_EQ(args.out_prefix, std::string("out"));
    ASSERT_EQ(args.min_reads, static_cast<int64_t>(100000));
    ASSERT_EQ(args.threads, 1);
    ASSERT_NEAR(args.alpha, 0.01, 1e-12);
    ASSERT_EQ(args.perms, 1000);
    ASSERT_EQ(args.min_seg_len, 25);
    ASSERT_EQ(args.max_depth, 100);
    ASSERT_TRUE(args.cbs_method == CbsMethod::OneCp);
    ASSERT_EQ(args.seed, static_cast<uint64_t>(1));
    ASSERT_TRUE(!args.quiet);
}

void test_overrides_all_flags() {
    auto [argc, argv] = make_argv({
        "dfn_cbs", "--counts", "in.txt.gz", "--bins", "bins.tsv", "--out-prefix", "out",
        "--min-reads", "5000", "--threads", "8", "--alpha", "0.05", "--perms", "200",
        "--min-seg-len", "10", "--max-depth", "50", "--cbs-method", "2cp",
        "--seed", "42", "--quiet"});
    CbsArgs args;
    bool help = false;
    ASSERT_TRUE(parse_cbs_args(argc, argv, args, help));
    ASSERT_EQ(args.min_reads, static_cast<int64_t>(5000));
    ASSERT_EQ(args.threads, 8);
    ASSERT_NEAR(args.alpha, 0.05, 1e-12);
    ASSERT_EQ(args.perms, 200);
    ASSERT_EQ(args.min_seg_len, 10);
    ASSERT_EQ(args.max_depth, 50);
    ASSERT_TRUE(args.cbs_method == CbsMethod::TwoCp);
    ASSERT_EQ(args.seed, static_cast<uint64_t>(42));
    ASSERT_TRUE(args.quiet);
}

void test_help_flag() {
    auto [argc, argv] = make_argv({"dfn_cbs", "--help"});
    CbsArgs args;
    bool help = false;
    ASSERT_TRUE(parse_cbs_args(argc, argv, args, help));
    ASSERT_TRUE(help);
}

void test_missing_required_fails() {
    auto [argc, argv] = make_argv({"dfn_cbs", "--bins", "bins.tsv", "--out-prefix", "out"});
    CbsArgs args;
    bool help = false;
    ASSERT_TRUE(!parse_cbs_args(argc, argv, args, help));
}

void test_unknown_flag_fails() {
    auto [argc, argv] = make_argv({"dfn_cbs", "--counts", "in.txt.gz", "--bins", "bins.tsv",
                                    "--out-prefix", "out", "--nonsense", "1"});
    CbsArgs args;
    bool help = false;
    ASSERT_TRUE(!parse_cbs_args(argc, argv, args, help));
}

void test_invalid_cbs_method_fails() {
    auto [argc, argv] = make_argv({"dfn_cbs", "--counts", "in.txt.gz", "--bins", "bins.tsv",
                                    "--out-prefix", "out", "--cbs-method", "3cp"});
    CbsArgs args;
    bool help = false;
    ASSERT_TRUE(!parse_cbs_args(argc, argv, args, help));
}

int main() {
    test_defaults_with_required_only();
    test_overrides_all_flags();
    test_help_flag();
    test_missing_required_fails();
    test_unknown_flag_fails();
    test_invalid_cbs_method_fails();
    TEST_REPORT();
}
