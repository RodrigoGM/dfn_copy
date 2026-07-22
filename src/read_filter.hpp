#pragma once
#include <htslib/sam.h>
#include "cli_args.hpp"

enum class DiscordantReason {
    None,
    SplitChimeric,
    CrossChromosome,
    ImproperPair
};

bool is_unconditionally_dropped(const bam1_t* rec);

// Must only be called on records for which is_unconditionally_dropped
// returned false.
DiscordantReason classify_discordant(const bam1_t* rec, const Args& args);

// Must only be called on records for which classify_discordant returned
// DiscordantReason::None.
bool passes_quality_filters(const bam1_t* rec, const Args& args);
