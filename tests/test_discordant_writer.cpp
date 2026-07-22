#include "test_framework.hpp"
#include "discordant_writer.hpp"
#include <htslib/sam.h>
#include <sys/stat.h>

bool file_exists(const std::string& path) {
    struct stat buf;
    return stat(path.c_str(), &buf) == 0;
}

void test_write_and_read_back() {
    sam_hdr_t* hdr = sam_hdr_init();
    sam_hdr_add_line(hdr, "SQ", "SN", "chr1", "LN", "1000", NULL);

    std::string path = "/tmp/test_discordant_writer_out.bam";
    {
        DiscordantWriter writer(path, hdr);
        bam1_t* rec1 = make_test_record("split1", BAM_FSUPPLEMENTARY, 0, 100, 60, -1, -1, 0, 50);
        bam1_t* rec2 = make_test_record("split2", BAM_FSUPPLEMENTARY, 0, 200, 60, -1, -1, 0, 50);
        writer.write(rec1);
        writer.write(rec2);
        bam_destroy1(rec1);
        bam_destroy1(rec2);
        writer.close();
    }

    ASSERT_TRUE(file_exists(path));
    ASSERT_TRUE(file_exists(path + ".bai"));

    samFile* in = sam_open(path.c_str(), "r");
    sam_hdr_t* read_hdr = sam_hdr_read(in);
    bam1_t* rec = bam_init1();

    int count = 0;
    while (sam_read1(in, read_hdr, rec) >= 0) {
        count++;
    }
    ASSERT_EQ(count, 2);

    bam_destroy1(rec);
    sam_hdr_destroy(read_hdr);
    sam_close(in);
    sam_hdr_destroy(hdr);
}

int main() {
    test_write_and_read_back();
    TEST_REPORT();
}
