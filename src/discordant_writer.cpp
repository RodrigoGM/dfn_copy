#include "discordant_writer.hpp"
#include <stdexcept>

DiscordantWriter::DiscordantWriter(const std::string& path, sam_hdr_t* hdr)
    : path_(path), hdr_(hdr) {
    fp_ = sam_open(path.c_str(), "wb");
    if (!fp_) {
        throw std::runtime_error("cannot open discordant BAM for writing: " + path);
    }
    if (sam_hdr_write(fp_, hdr_) < 0) {
        sam_close(fp_);
        throw std::runtime_error("cannot write header to discordant BAM: " + path);
    }
}

void DiscordantWriter::write(const bam1_t* rec) {
    if (sam_write1(fp_, hdr_, rec) < 0) {
        throw std::runtime_error("failed writing record to discordant BAM: " + path_);
    }
}

void DiscordantWriter::close() {
    if (closed_) return;
    sam_close(fp_);
    fp_ = nullptr;
    if (sam_index_build(path_.c_str(), 0) < 0) {
        throw std::runtime_error("failed to index discordant BAM: " + path_);
    }
    closed_ = true;
}

DiscordantWriter::~DiscordantWriter() {
    if (!closed_ && fp_) {
        sam_close(fp_);
    }
}
