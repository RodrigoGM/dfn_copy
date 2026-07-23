#pragma once
#include "io_types.hpp"
#include <random>
#include <vector>

void segment_chromosome(const Series& S, const Args& a, std::vector<Segment>& out, std::mt19937_64& rng);
