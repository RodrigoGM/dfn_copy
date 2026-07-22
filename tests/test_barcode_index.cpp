#include "test_framework.hpp"
#include "barcode_index.hpp"
#include <fstream>

std::string write_temp_file(const std::string& name, const std::string& content) {
    std::string path = "/tmp/" + name;
    std::ofstream out(path);
    out << content;
    out.close();
    return path;
}

void test_first_seen_order() {
    BarcodeIndex idx;
    ASSERT_TRUE(!idx.has_allowlist());

    auto a = idx.get_or_create("AAAA-1");
    auto b = idx.get_or_create("CCCC-1");
    auto a2 = idx.get_or_create("AAAA-1");

    ASSERT_TRUE(a.has_value());
    ASSERT_TRUE(b.has_value());
    ASSERT_EQ(*a, static_cast<size_t>(0));
    ASSERT_EQ(*b, static_cast<size_t>(1));
    ASSERT_EQ(*a2, static_cast<size_t>(0));
    ASSERT_EQ(idx.size(), static_cast<size_t>(2));
    ASSERT_EQ(idx.names()[0], std::string("AAAA-1"));
    ASSERT_EQ(idx.names()[1], std::string("CCCC-1"));
}

void test_allowlist_order_and_rejection() {
    std::string path = write_temp_file("allow.txt", "CCCC-1\nAAAA-1\n");
    BarcodeIndex idx = BarcodeIndex::load_allowlist(path);
    ASSERT_TRUE(idx.has_allowlist());
    ASSERT_EQ(idx.size(), static_cast<size_t>(2));
    ASSERT_EQ(idx.names()[0], std::string("CCCC-1"));
    ASSERT_EQ(idx.names()[1], std::string("AAAA-1"));

    auto c = idx.get_or_create("CCCC-1");
    ASSERT_TRUE(c.has_value());
    ASSERT_EQ(*c, static_cast<size_t>(0));

    auto rejected = idx.get_or_create("GGGG-1");
    ASSERT_TRUE(!rejected.has_value());
    ASSERT_EQ(idx.size(), static_cast<size_t>(2));
}

void test_empty_allowlist_throws() {
    std::string path = write_temp_file("allow_empty.txt", "");
    bool threw = false;
    try {
        BarcodeIndex::load_allowlist(path);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);
}

void test_unreadable_allowlist_throws() {
    bool threw = false;
    try {
        BarcodeIndex::load_allowlist("/tmp/does_not_exist_allowlist.txt");
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);
}

int main() {
    test_first_seen_order();
    test_allowlist_order_and_rejection();
    test_empty_allowlist_throws();
    test_unreadable_allowlist_throws();
    TEST_REPORT();
}
