#pragma once
#include <string>
#include <htslib/sam.h>

class DiscordantWriter {
public:
    DiscordantWriter(const std::string& path, sam_hdr_t* hdr);
    ~DiscordantWriter();

    DiscordantWriter(const DiscordantWriter&) = delete;
    DiscordantWriter& operator=(const DiscordantWriter&) = delete;

    void write(const bam1_t* rec);
    void close();

private:
    std::string path_;
    samFile* fp_ = nullptr;
    sam_hdr_t* hdr_ = nullptr;
    bool closed_ = false;
};
