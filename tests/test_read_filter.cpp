#include "test_framework.hpp"
#include "bam_test_helpers.hpp"
#include "read_filter.hpp"

void test_unconditional_drops() {
    bam1_t* unmapped = make_test_record("r1", BAM_FUNMAP, -1, -1, 0, -1, -1, 0, 50);
    ASSERT_TRUE(is_unconditionally_dropped(unmapped));
    bam_destroy1(unmapped);

    bam1_t* qcfail = make_test_record("r2", BAM_FQCFAIL, 0, 100, 60, -1, -1, 0, 50);
    ASSERT_TRUE(is_unconditionally_dropped(qcfail));
    bam_destroy1(qcfail);

    bam1_t* clean = make_test_record("r3", 0, 0, 100, 60, -1, -1, 0, 50);
    ASSERT_TRUE(!is_unconditionally_dropped(clean));
    bam_destroy1(clean);
}

void test_classify_supplementary_is_split() {
    bam1_t* rec = make_test_record("r1", BAM_FSUPPLEMENTARY, 1, 100, 60, -1, -1, 0, 50);
    Args args;
    ASSERT_TRUE(classify_discordant(rec, args) == DiscordantReason::SplitChimeric);
    bam_destroy1(rec);
}

void test_classify_sa_tag_is_split() {
    bam1_t* rec = make_test_record("r1", 0, 0, 500, 60, -1, -1, 0, 50);
    std::string sa_value = "chr2,101,+,50M,60,0;";
    bam_aux_append(rec, "SA", 'Z', static_cast<int>(sa_value.size() + 1),
                    reinterpret_cast<const uint8_t*>(sa_value.c_str()));
    Args args;
    ASSERT_TRUE(classify_discordant(rec, args) == DiscordantReason::SplitChimeric);
    bam_destroy1(rec);
}

void test_classify_cross_chromosome() {
    bam1_t* rec = make_test_record("r1", BAM_FPAIRED | BAM_FREAD1, 0, 100, 60, 1, 200, 0, 50);
    Args args;
    ASSERT_TRUE(classify_discordant(rec, args) == DiscordantReason::CrossChromosome);
    bam_destroy1(rec);
}

void test_classify_improper_pair_flag_unset() {
    bam1_t* rec = make_test_record("r1", BAM_FPAIRED | BAM_FREAD1, 0, 100, 60, 0, 200, 150, 50);
    Args args;
    ASSERT_TRUE(classify_discordant(rec, args) == DiscordantReason::ImproperPair);
    bam_destroy1(rec);
}

void test_classify_proper_pair_is_clean() {
    bam1_t* rec = make_test_record("r1", BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREAD1,
                                    0, 100, 60, 0, 300, 250, 50);
    Args args;
    ASSERT_TRUE(classify_discordant(rec, args) == DiscordantReason::None);
    bam_destroy1(rec);
}

void test_classify_max_insert_size_override() {
    bam1_t* rec = make_test_record("r1", BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREAD1,
                                    0, 100, 60, 0, 50100, 50000, 50);
    Args args;
    args.max_insert_size = 0;
    ASSERT_TRUE(classify_discordant(rec, args) == DiscordantReason::None);

    args.max_insert_size = 10000;
    ASSERT_TRUE(classify_discordant(rec, args) == DiscordantReason::ImproperPair);
    bam_destroy1(rec);
}

void test_classify_unpaired_is_clean() {
    bam1_t* rec = make_test_record("r1", 0, 0, 100, 60, -1, -1, 0, 50);
    Args args;
    ASSERT_TRUE(classify_discordant(rec, args) == DiscordantReason::None);
    bam_destroy1(rec);
}

void test_classify_mate_unmapped_is_clean() {
    bam1_t* rec = make_test_record("r1", BAM_FPAIRED | BAM_FMUNMAP | BAM_FREAD1,
                                    0, 100, 60, -1, -1, 0, 50);
    Args args;
    ASSERT_TRUE(classify_discordant(rec, args) == DiscordantReason::None);
    bam_destroy1(rec);
}

void test_quality_filters_mapq() {
    bam1_t* rec = make_test_record("r1", 0, 0, 100, 10, -1, -1, 0, 50);
    Args args;
    args.mapq = 30;
    ASSERT_TRUE(!passes_quality_filters(rec, args));
    args.mapq = 5;
    ASSERT_TRUE(passes_quality_filters(rec, args));
    bam_destroy1(rec);
}

void test_quality_filters_secondary() {
    bam1_t* rec = make_test_record("r1", BAM_FSECONDARY, 0, 100, 60, -1, -1, 0, 50);
    Args args;
    args.primary_alignment_only = true;
    ASSERT_TRUE(!passes_quality_filters(rec, args));
    args.primary_alignment_only = false;
    ASSERT_TRUE(passes_quality_filters(rec, args));
    bam_destroy1(rec);
}

void test_quality_filters_dup() {
    bam1_t* rec = make_test_record("r1", BAM_FDUP, 0, 100, 60, -1, -1, 0, 50);
    Args args;
    args.exclude_dups = true;
    ASSERT_TRUE(!passes_quality_filters(rec, args));
    args.exclude_dups = false;
    ASSERT_TRUE(passes_quality_filters(rec, args));
    bam_destroy1(rec);
}

int main() {
    test_unconditional_drops();
    test_classify_supplementary_is_split();
    test_classify_sa_tag_is_split();
    test_classify_cross_chromosome();
    test_classify_improper_pair_flag_unset();
    test_classify_proper_pair_is_clean();
    test_classify_max_insert_size_override();
    test_classify_unpaired_is_clean();
    test_classify_mate_unmapped_is_clean();
    test_quality_filters_mapq();
    test_quality_filters_secondary();
    test_quality_filters_dup();
    TEST_REPORT();
}
