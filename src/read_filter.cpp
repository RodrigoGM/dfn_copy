#include "read_filter.hpp"
#include <cstdlib>

bool is_unconditionally_dropped(const bam1_t* rec) {
    return (rec->core.flag & (BAM_FUNMAP | BAM_FQCFAIL)) != 0;
}

DiscordantReason classify_discordant(const bam1_t* rec, const Args& args) {
    uint16_t flag = rec->core.flag;

    if ((flag & BAM_FSUPPLEMENTARY) || bam_aux_get(rec, "SA") != nullptr) {
        return DiscordantReason::SplitChimeric;
    }

    bool paired_and_mate_mapped = (flag & BAM_FPAIRED) && !(flag & BAM_FMUNMAP);
    if (paired_and_mate_mapped) {
        if (rec->core.mtid != rec->core.tid) {
            return DiscordantReason::CrossChromosome;
        }
        bool proper_pair = (flag & BAM_FPROPER_PAIR) != 0;
        bool oversized = args.max_insert_size > 0 &&
            std::llabs(static_cast<long long>(rec->core.isize)) > args.max_insert_size;
        if (!proper_pair || oversized) {
            return DiscordantReason::ImproperPair;
        }
    }

    return DiscordantReason::None;
}

bool passes_quality_filters(const bam1_t* rec, const Args& args) {
    if (rec->core.qual < args.mapq) return false;
    if (args.primary_alignment_only && (rec->core.flag & BAM_FSECONDARY)) return false;
    if (args.exclude_dups && (rec->core.flag & BAM_FDUP)) return false;
    return true;
}
