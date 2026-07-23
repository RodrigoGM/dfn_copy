#pragma once
#include <string>
#include <vector>
#include <unordered_map>

// Trimmed extract of cbs+'s include/io.hpp -- see README.md in this
// directory. Only the plain structs segment_chromosome() needs.

struct Series {
    std::string chrom;
    std::vector<long> pos;
    std::vector<double> x;
    std::unordered_map<std::string, std::vector<double>> covars;
};

struct Segment { int s, e; double mean; int level; };

struct Args {
    std::string in = "";
    std::string delimiter = "auto";
    double alpha = 0.01;
    int perms = 1000;
    int min_seg_len = 25;
    int max_depth = 100;
    unsigned long long seed = 1;
    bool center = true;
    int threads = 1;
    std::string method = "1cp";
    std::string out_format = "tsv";
    std::string out = "";
    std::string gc_col = "";
    std::string wave_col = "";
    std::string correction = "none";
    int roll_window = 101;
};
