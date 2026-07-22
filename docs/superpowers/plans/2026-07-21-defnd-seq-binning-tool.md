# DEFND-seq Binning Tool Implementation Plan

> Note (2026-07-22): the binary described here as `dfn_bin` was renamed to
> `dfn_copy` after this plan was written. Left as-is below as a historical
> record; see README.md and CHANGELOG.txt for the current name.

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `dfn_bin`, a C++17/htslib command-line tool that streams a coordinate-sorted DEFND-seq gDNA BAM once and emits a raw bin×cell fragment-count matrix plus a discordant/split-read BAM, and `correct_gc.py`, a companion script that GC-corrects the raw matrix.

**Architecture:** `dfn_bin` is one streaming pass over the BAM: classify each record (unconditional drop / discordant-diverted / normal) → for normal paired reads, resolve to one fragment representative via a QNAME lookback buffer → bin the survivor into a dense `bins × barcodes` `int32` matrix → write it as a gzipped TSV. Discordant/split reads are diverted to a second output BAM as encountered, never entering the matrix. `correct_gc.py` is a separate, independent script (no shared code with the C++ tool) that reads the raw matrix + the bins file's GC column and writes a GC-corrected matrix.

**Tech Stack:** C++17, htslib (via pkg-config, confirmed installed: 1.22.1), zlib, GNU Make; Python 3 (pysam for test fixtures, pytest for integration tests, numpy/statsmodels for `correct_gc.py`).

## Global Constraints

- Design spec: `docs/superpowers/specs/2026-07-21-defnd-seq-binning-design.md` — every task below implements a specific section of it; re-read it if a detail here seems ambiguous.
- C++17, htslib linked via `pkg-config --cflags --libs htslib`. No other C++ dependencies besides zlib (`-lz`).
- No config files, no plugin system, no output formats beyond the three specified: `<prefix>.raw_counts.txt.gz`, `<prefix>.gc_corrected.txt.gz`, `<prefix>.discordant.bam` (+`.bai`).
- Bin identifiers in output are `<chrom>:<start>:<end>` (colon-separated).
- Barcode comes only from the configured BAM tag (default `CB`) — never parsed from read sequence.
- All CLI flags and defaults are fixed by the spec's CLI contract table — do not add, rename, or change defaults.
- htslib API calls in this plan are verified against the installed htslib 1.22.1 headers (`/opt/Homebrew/Cellar/htslib/1.22.1/include/htslib/sam.h`) — signatures are exact, not from memory.
- Python test/companion-script dependencies (`pysam`, `pytest`, `numpy`, `statsmodels`) are not installed system-wide; they're installed into a project-local `.venv` the first time they're needed (Task 9).

---

## File Structure

```
dfn_copy/
├── Makefile
├── .gitignore
├── src/
│   ├── main.cpp              # CLI entry point + full streaming pipeline
│   ├── cli_args.hpp/.cpp      # Args struct, parse_args(), usage_text()
│   ├── bins.hpp/.cpp          # Bin, BinIndex, compute_bin_coordinate, validate_chrom_naming
│   ├── barcode_index.hpp/.cpp # BarcodeIndex
│   ├── counts_matrix.hpp/.cpp # CountsMatrix (dense, growable, gzip TSV writer)
│   ├── read_filter.hpp/.cpp   # is_unconditionally_dropped, classify_discordant, passes_quality_filters
│   ├── fragment_pairing.hpp/.cpp # FragmentCandidate, FragmentPairingBuffer
│   └── discordant_writer.hpp/.cpp # DiscordantWriter
├── tests/
│   ├── test_framework.hpp     # hand-rolled ASSERT macros + make_test_record() helper
│   ├── test_cli_args.cpp
│   ├── test_bins.cpp
│   ├── test_barcode_index.cpp
│   ├── test_counts_matrix.cpp
│   ├── test_read_filter.cpp
│   ├── test_fragment_pairing.cpp
│   ├── test_discordant_writer.cpp
│   ├── fixtures/
│   │   └── gen_fixtures.py    # pysam script building all integration-test BAMs/bins files
│   └── integration/
│       ├── test_smoke.py      # Task 9's minimal end-to-end check
│       └── test_dfn_bin.py    # Task 10's full scenario coverage
├── correct_gc.py
├── tests_py/
│   └── test_correct_gc.py
└── README.md
```

---

### Task 1: Project scaffold, Makefile, build sanity check

**Files:**
- Create: `.gitignore`
- Create: `Makefile`
- Create: `src/main.cpp` (trivial stub, replaced in Task 2 and rewritten in Task 9)

**Interfaces:**
- Produces: a working `make` build producing binary `./dfn_bin`, and a `make test` target (initially a no-op, extended by every later task).

- [ ] **Step 1: Create directory structure and `.gitignore`**

```bash
mkdir -p src tests/fixtures tests/integration tests_py
```

`.gitignore`:
```
.venv/
build/
tests/build/
*.o
dfn_bin
__pycache__/
*.pyc
tests/fixtures/*.bam
tests/fixtures/*.bai
tests/fixtures/*.bam.bai
```

- [ ] **Step 2: Write the Makefile**

`Makefile`:
```makefile
CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -O2 -Isrc $(shell pkg-config --cflags htslib)
LDLIBS := $(shell pkg-config --libs htslib) -lz

BUILD_DIR := build
TEST_BUILD := tests/build
BIN := dfn_bin

.PHONY: all clean test

all: $(BIN)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(TEST_BUILD):
	mkdir -p $(TEST_BUILD)

$(BUILD_DIR)/main.o: src/main.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN): $(BUILD_DIR)/main.o
	$(CXX) $(BUILD_DIR)/main.o -o $@ $(LDLIBS)

test:
	@echo "no tests yet"

clean:
	rm -rf $(BUILD_DIR) $(TEST_BUILD) $(BIN)
```

- [ ] **Step 3: Write a trivial stub `main.cpp`**

`src/main.cpp`:
```cpp
int main() {
    return 0;
}
```

- [ ] **Step 4: Build and verify**

Run: `make`
Expected: compiles cleanly, produces `./dfn_bin`.

Run: `./dfn_bin; echo $?`
Expected: `0`

- [ ] **Step 5: Commit**

```bash
git add .gitignore Makefile src/main.cpp
git commit -m "Scaffold dfn_bin project: directories, Makefile, build sanity check"
```

---

### Task 2: CLI argument parsing

**Files:**
- Create: `src/cli_args.hpp`
- Create: `src/cli_args.cpp`
- Test: `tests/test_cli_args.cpp`
- Create: `tests/test_framework.hpp` (assertion macros only in this task; `make_test_record()` is added in Task 6)
- Modify: `Makefile` (add `cli_args.o` build rule, `test_cli_args` test rule)
- Modify: `src/main.cpp` (call `parse_args`, handle `--help` and errors)

**Interfaces:**
- Produces: `enum class PositionMode { Start, Midpoint, End }`, `enum class SingleEndCountingMode { Auto, R1, R2 }`, `struct Args { ... }` (full field list below), `bool parse_args(int argc, char** argv, Args& args, bool& help_requested)`, `std::string usage_text()`. Every later task that touches `bam1_t`/filters/pairing consumes `Args` and these two enums.

- [ ] **Step 1: Write the failing test**

`tests/test_framework.hpp`:
```cpp
#pragma once
#include <cstdio>

inline int g_tests_run = 0;
inline int g_tests_failed = 0;

#define ASSERT_TRUE(cond) do { \
    g_tests_run++; \
    if (!(cond)) { \
        std::printf("FAIL %s:%d: ASSERT_TRUE(%s)\n", __FILE__, __LINE__, #cond); \
        g_tests_failed++; \
    } \
} while (0)

#define ASSERT_EQ(a, b) do { \
    g_tests_run++; \
    auto _a = (a); auto _b = (b); \
    if (!(_a == _b)) { \
        std::printf("FAIL %s:%d: ASSERT_EQ(%s, %s)\n", __FILE__, __LINE__, #a, #b); \
        g_tests_failed++; \
    } \
} while (0)

#define TEST_REPORT() do { \
    std::printf("%d assertions, %d failed\n", g_tests_run, g_tests_failed); \
    return g_tests_failed == 0 ? 0 : 1; \
} while (0)
```

`tests/test_cli_args.cpp`:
```cpp
#include "test_framework.hpp"
#include "cli_args.hpp"
#include <vector>
#include <cstring>

std::pair<int, char**> make_argv(std::vector<std::string> args) {
    static std::vector<std::string> storage;
    static std::vector<char*> ptrs;
    storage = args;
    ptrs.clear();
    for (auto& s : storage) ptrs.push_back(const_cast<char*>(s.c_str()));
    return {static_cast<int>(ptrs.size()), ptrs.data()};
}

void test_defaults_with_required_only() {
    auto [argc, argv] = make_argv({"dfn_bin", "--bam", "in.bam", "--bins", "bins.tsv",
                                    "--out-prefix", "out"});
    Args args;
    bool help = false;
    ASSERT_TRUE(parse_args(argc, argv, args, help));
    ASSERT_TRUE(!help);
    ASSERT_EQ(args.bam_path, std::string("in.bam"));
    ASSERT_EQ(args.bins_path, std::string("bins.tsv"));
    ASSERT_EQ(args.out_prefix, std::string("out"));
    ASSERT_EQ(args.barcode_tag, std::string("CB"));
    ASSERT_EQ(args.barcodes_path, std::string(""));
    ASSERT_EQ(args.mapq, 30);
    ASSERT_TRUE(args.exclude_dups);
    ASSERT_TRUE(args.primary_alignment_only);
    ASSERT_TRUE(args.single_end_counting == SingleEndCountingMode::Auto);
    ASSERT_EQ(args.max_insert_size, static_cast<int64_t>(0));
    ASSERT_TRUE(args.position == PositionMode::End);
    ASSERT_EQ(args.threads, 1);
}

void test_overrides_all_flags() {
    auto [argc, argv] = make_argv({
        "dfn_bin", "--bam", "in.bam", "--bins", "bins.tsv", "--out-prefix", "out",
        "--barcode-tag", "CR", "--barcodes", "allow.txt", "--mapq", "20",
        "--exclude-dups", "false", "--primary-alignment-only", "false",
        "--single-end-counting", "r2", "--max-insert-size", "5000",
        "--position", "start", "--threads", "4"});
    Args args;
    bool help = false;
    ASSERT_TRUE(parse_args(argc, argv, args, help));
    ASSERT_EQ(args.barcode_tag, std::string("CR"));
    ASSERT_EQ(args.barcodes_path, std::string("allow.txt"));
    ASSERT_EQ(args.mapq, 20);
    ASSERT_TRUE(!args.exclude_dups);
    ASSERT_TRUE(!args.primary_alignment_only);
    ASSERT_TRUE(args.single_end_counting == SingleEndCountingMode::R2);
    ASSERT_EQ(args.max_insert_size, static_cast<int64_t>(5000));
    ASSERT_TRUE(args.position == PositionMode::Start);
    ASSERT_EQ(args.threads, 4);
}

void test_help_flag() {
    auto [argc, argv] = make_argv({"dfn_bin", "--help"});
    Args args;
    bool help = false;
    ASSERT_TRUE(parse_args(argc, argv, args, help));
    ASSERT_TRUE(help);
}

void test_missing_required_fails() {
    auto [argc, argv] = make_argv({"dfn_bin", "--bins", "bins.tsv", "--out-prefix", "out"});
    Args args;
    bool help = false;
    ASSERT_TRUE(!parse_args(argc, argv, args, help));
}

void test_unknown_flag_fails() {
    auto [argc, argv] = make_argv({"dfn_bin", "--bam", "in.bam", "--bins", "bins.tsv",
                                    "--out-prefix", "out", "--nonsense", "1"});
    Args args;
    bool help = false;
    ASSERT_TRUE(!parse_args(argc, argv, args, help));
}

void test_invalid_enum_value_fails() {
    auto [argc, argv] = make_argv({"dfn_bin", "--bam", "in.bam", "--bins", "bins.tsv",
                                    "--out-prefix", "out", "--position", "middle"});
    Args args;
    bool help = false;
    ASSERT_TRUE(!parse_args(argc, argv, args, help));
}

void test_barcode_tag_must_be_two_chars() {
    auto [argc, argv] = make_argv({"dfn_bin", "--bam", "in.bam", "--bins", "bins.tsv",
                                    "--out-prefix", "out", "--barcode-tag", "TOOLONG"});
    Args args;
    bool help = false;
    ASSERT_TRUE(!parse_args(argc, argv, args, help));
}

int main() {
    test_defaults_with_required_only();
    test_overrides_all_flags();
    test_help_flag();
    test_missing_required_fails();
    test_unknown_flag_fails();
    test_invalid_enum_value_fails();
    test_barcode_tag_must_be_two_chars();
    TEST_REPORT();
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++17 -Isrc -c tests/test_cli_args.cpp -o /tmp/test_cli_args.o`
Expected: FAIL — `cli_args.hpp` does not exist yet.

- [ ] **Step 3: Write `src/cli_args.hpp`**

```cpp
#pragma once
#include <cstdint>
#include <string>

enum class PositionMode { Start, Midpoint, End };
enum class SingleEndCountingMode { Auto, R1, R2 };

struct Args {
    std::string bam_path;
    std::string bins_path;
    std::string out_prefix;
    std::string barcode_tag = "CB";
    std::string barcodes_path;
    int mapq = 30;
    bool exclude_dups = true;
    bool primary_alignment_only = true;
    SingleEndCountingMode single_end_counting = SingleEndCountingMode::Auto;
    int64_t max_insert_size = 0;
    PositionMode position = PositionMode::End;
    int threads = 1;
};

// Parses argv into args. Returns true on success (including a bare
// --help, which sets help_requested and leaves other fields unset).
// On error, prints a specific message to stderr and returns false.
bool parse_args(int argc, char** argv, Args& args, bool& help_requested);

std::string usage_text();
```

- [ ] **Step 4: Write `src/cli_args.cpp`**

```cpp
#include "cli_args.hpp"
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

std::string usage_text() {
    return
        "dfn_bin --bam <path> --bins <path> --out-prefix <prefix> [options]\n"
        "\n"
        "Required:\n"
        "  --bam <path>            Coordinate-sorted, indexed DEFND-seq gDNA BAM\n"
        "  --bins <path>           Tab-separated bin file (header columns: chrom,start,end,gc)\n"
        "  --out-prefix <prefix>   Prefix for output files\n"
        "\n"
        "Options:\n"
        "  --barcode-tag <tag>            BAM tag holding the cell barcode (default: CB)\n"
        "  --barcodes <path>              Optional barcode allowlist file (default: none)\n"
        "  --mapq <int>                   Minimum MAPQ to keep a read (default: 30)\n"
        "  --exclude-dups <true|false>    Exclude PCR/optical duplicates (default: true)\n"
        "  --primary-alignment-only <true|false>  Exclude secondary alignments (default: true)\n"
        "  --single-end-counting <auto|r1|r2>     Fragment-counting mode (default: auto)\n"
        "  --max-insert-size <int>        Absolute discordance cutoff, 0=disabled (default: 0)\n"
        "  --position <start|midpoint|end> Bin-assignment coordinate (default: end)\n"
        "  --threads <int>                htslib decompression threads (default: 1)\n"
        "  --help                         Print this message\n";
}

namespace {

bool parse_bool(const std::string& s, bool& out) {
    if (s == "true") { out = true; return true; }
    if (s == "false") { out = false; return true; }
    return false;
}

bool parse_int(const std::string& s, int& out) {
    try {
        size_t consumed = 0;
        out = std::stoi(s, &consumed);
        return consumed == s.size();
    } catch (const std::exception&) {
        return false;
    }
}

bool parse_int64(const std::string& s, int64_t& out) {
    try {
        size_t consumed = 0;
        out = std::stoll(s, &consumed);
        return consumed == s.size();
    } catch (const std::exception&) {
        return false;
    }
}

} // namespace

bool parse_args(int argc, char** argv, Args& args, bool& help_requested) {
    help_requested = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--help") {
            help_requested = true;
            return true;
        }

        auto next_value = [&](const char* flag_name) -> const char* {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "error: %s requires a value\n", flag_name);
                return nullptr;
            }
            return argv[++i];
        };

        if (arg == "--bam") {
            const char* v = next_value("--bam"); if (!v) return false;
            args.bam_path = v;
        } else if (arg == "--bins") {
            const char* v = next_value("--bins"); if (!v) return false;
            args.bins_path = v;
        } else if (arg == "--out-prefix") {
            const char* v = next_value("--out-prefix"); if (!v) return false;
            args.out_prefix = v;
        } else if (arg == "--barcode-tag") {
            const char* v = next_value("--barcode-tag"); if (!v) return false;
            args.barcode_tag = v;
        } else if (arg == "--barcodes") {
            const char* v = next_value("--barcodes"); if (!v) return false;
            args.barcodes_path = v;
        } else if (arg == "--mapq") {
            const char* v = next_value("--mapq"); if (!v) return false;
            if (!parse_int(v, args.mapq)) {
                std::fprintf(stderr, "error: --mapq requires an integer, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--exclude-dups") {
            const char* v = next_value("--exclude-dups"); if (!v) return false;
            if (!parse_bool(v, args.exclude_dups)) {
                std::fprintf(stderr, "error: --exclude-dups requires true|false, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--primary-alignment-only") {
            const char* v = next_value("--primary-alignment-only"); if (!v) return false;
            if (!parse_bool(v, args.primary_alignment_only)) {
                std::fprintf(stderr, "error: --primary-alignment-only requires true|false, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--single-end-counting") {
            const char* v = next_value("--single-end-counting"); if (!v) return false;
            std::string sv = v;
            if (sv == "auto") args.single_end_counting = SingleEndCountingMode::Auto;
            else if (sv == "r1") args.single_end_counting = SingleEndCountingMode::R1;
            else if (sv == "r2") args.single_end_counting = SingleEndCountingMode::R2;
            else {
                std::fprintf(stderr, "error: --single-end-counting requires auto|r1|r2, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--max-insert-size") {
            const char* v = next_value("--max-insert-size"); if (!v) return false;
            if (!parse_int64(v, args.max_insert_size)) {
                std::fprintf(stderr, "error: --max-insert-size requires an integer, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--position") {
            const char* v = next_value("--position"); if (!v) return false;
            std::string sv = v;
            if (sv == "start") args.position = PositionMode::Start;
            else if (sv == "midpoint") args.position = PositionMode::Midpoint;
            else if (sv == "end") args.position = PositionMode::End;
            else {
                std::fprintf(stderr, "error: --position requires start|midpoint|end, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--threads") {
            const char* v = next_value("--threads"); if (!v) return false;
            if (!parse_int(v, args.threads)) {
                std::fprintf(stderr, "error: --threads requires an integer, got '%s'\n", v);
                return false;
            }
        } else {
            std::fprintf(stderr, "error: unknown flag '%s'\n", arg.c_str());
            return false;
        }
    }

    if (args.bam_path.empty()) {
        std::fprintf(stderr, "error: --bam is required\n");
        return false;
    }
    if (args.bins_path.empty()) {
        std::fprintf(stderr, "error: --bins is required\n");
        return false;
    }
    if (args.out_prefix.empty()) {
        std::fprintf(stderr, "error: --out-prefix is required\n");
        return false;
    }
    if (args.barcode_tag.size() != 2) {
        std::fprintf(stderr, "error: --barcode-tag must be exactly 2 characters, got '%s'\n",
                      args.barcode_tag.c_str());
        return false;
    }

    return true;
}
```

- [ ] **Step 5: Run test to verify it passes**

Run: `g++ -std=c++17 -Isrc tests/test_cli_args.cpp src/cli_args.cpp -o /tmp/test_cli_args && /tmp/test_cli_args`
Expected: `7 tests` worth of assertions, `0 failed`, exit code 0.

- [ ] **Step 6: Wire the Makefile**

Append to `Makefile` (before the `clean:` target):
```makefile
$(BUILD_DIR)/cli_args.o: src/cli_args.cpp src/cli_args.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_BUILD)/test_cli_args: tests/test_cli_args.cpp $(BUILD_DIR)/cli_args.o | $(TEST_BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)
```

Replace the `$(BIN):` rule and `test:` target:
```makefile
$(BIN): $(BUILD_DIR)/main.o $(BUILD_DIR)/cli_args.o
	$(CXX) $(BUILD_DIR)/main.o $(BUILD_DIR)/cli_args.o -o $@ $(LDLIBS)

TESTS := $(TEST_BUILD)/test_cli_args

.PHONY: test
test: $(TESTS)
	@for t in $(TESTS); do echo "=== $$t ==="; $$t || exit 1; done
```

- [ ] **Step 7: Wire `main.cpp` to use `parse_args`**

Replace `src/main.cpp`:
```cpp
#include "cli_args.hpp"
#include <cstdio>

int main(int argc, char** argv) {
    Args args;
    bool help_requested = false;

    if (!parse_args(argc, argv, args, help_requested)) {
        return 1;
    }
    if (help_requested) {
        std::fputs(usage_text().c_str(), stdout);
        return 0;
    }

    // Full pipeline is wired in Task 9.
    return 0;
}
```

- [ ] **Step 8: Run full build and test**

Run: `make clean && make && make test`
Expected: builds cleanly; test output shows `0 failed`.

Run: `./dfn_bin --help`
Expected: prints the usage text.

- [ ] **Step 9: Commit**

```bash
git add src/cli_args.hpp src/cli_args.cpp src/main.cpp tests/test_framework.hpp tests/test_cli_args.cpp Makefile
git commit -m "Add CLI argument parsing with full flag set and --help"
```

---

### Task 3: Bin file parsing, bin index, position computation, chromosome-naming validation

**Files:**
- Create: `src/bins.hpp`
- Create: `src/bins.cpp`
- Test: `tests/test_bins.cpp`
- Modify: `Makefile`

**Interfaces:**
- Consumes: `PositionMode` from `cli_args.hpp` (Task 2).
- Produces: `struct Bin { std::string chrom; int64_t start; int64_t end; double gc; }`, `class BinIndex` with `static BinIndex load(const std::string& path)` (throws `std::runtime_error` on malformed input), `size_t size() const`, `std::string name(size_t i) const`, `std::optional<size_t> find(const std::string& chrom, int64_t coord) const`, `const std::vector<Bin>& bins() const`; free functions `int64_t compute_bin_coordinate(int64_t align_start, int64_t align_end, PositionMode mode)` and `std::string validate_chrom_naming(const std::vector<std::string>& bam_contigs, const std::vector<std::string>& bin_chroms)` (returns `""` if consistent). Task 9 (main orchestration) and Task 10 (integration tests) rely on all of these exact names.

- [ ] **Step 1: Write the failing test**

`tests/test_bins.cpp`:
```cpp
#include "test_framework.hpp"
#include "bins.hpp"
#include <fstream>
#include <cstdio>

std::string write_temp_file(const std::string& name, const std::string& content) {
    std::string path = "/tmp/" + name;
    std::ofstream out(path);
    out << content;
    out.close();
    return path;
}

void test_load_basic() {
    std::string path = write_temp_file("bins_basic.tsv",
        "chrom\tstart\tend\tgc\n"
        "chr1\t0\t1000\t0.4\n"
        "chr1\t1000\t2000\t0.5\n"
        "chr2\t0\t1000\t0.45\n");
    BinIndex idx = BinIndex::load(path);
    ASSERT_EQ(idx.size(), static_cast<size_t>(3));
    ASSERT_EQ(idx.name(0), std::string("chr1:0:1000"));
    ASSERT_EQ(idx.name(2), std::string("chr2:0:1000"));
}

void test_load_columns_out_of_order() {
    std::string path = write_temp_file("bins_reorder.tsv",
        "gc\tend\tstart\tchrom\n"
        "0.4\t1000\t0\tchr1\n");
    BinIndex idx = BinIndex::load(path);
    ASSERT_EQ(idx.size(), static_cast<size_t>(1));
    ASSERT_EQ(idx.name(0), std::string("chr1:0:1000"));
}

void test_load_missing_column_throws() {
    std::string path = write_temp_file("bins_missing_col.tsv",
        "chrom\tstart\tend\n"
        "chr1\t0\t1000\n");
    bool threw = false;
    try {
        BinIndex::load(path);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);
}

void test_load_malformed_number_throws() {
    std::string path = write_temp_file("bins_malformed.tsv",
        "chrom\tstart\tend\tgc\n"
        "chr1\tNOTANUMBER\t1000\t0.4\n");
    bool threw = false;
    try {
        BinIndex::load(path);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);
}

void test_find_lookup() {
    std::string path = write_temp_file("bins_lookup.tsv",
        "chrom\tstart\tend\tgc\n"
        "chr1\t0\t1000\t0.4\n"
        "chr1\t1000\t2000\t0.5\n"
        "chr2\t0\t1000\t0.45\n");
    BinIndex idx = BinIndex::load(path);

    auto r1 = idx.find("chr1", 500);
    ASSERT_TRUE(r1.has_value());
    ASSERT_EQ(*r1, static_cast<size_t>(0));

    auto r2 = idx.find("chr1", 1500);
    ASSERT_TRUE(r2.has_value());
    ASSERT_EQ(*r2, static_cast<size_t>(1));

    auto r3 = idx.find("chr1", 999999);
    ASSERT_TRUE(!r3.has_value());

    auto r4 = idx.find("chrX", 10);
    ASSERT_TRUE(!r4.has_value());
}

void test_compute_bin_coordinate() {
    ASSERT_EQ(compute_bin_coordinate(100, 150, PositionMode::Start), static_cast<int64_t>(100));
    ASSERT_EQ(compute_bin_coordinate(100, 150, PositionMode::End), static_cast<int64_t>(149));
    ASSERT_EQ(compute_bin_coordinate(100, 150, PositionMode::Midpoint), static_cast<int64_t>(125));
}

void test_validate_chrom_naming_ok() {
    std::string err = validate_chrom_naming({"chr1", "chr2", "chr3"}, {"chr1", "chr2"});
    ASSERT_EQ(err, std::string(""));
}

void test_validate_chrom_naming_chr_prefix_mismatch() {
    std::string err = validate_chrom_naming({"chr1", "chr2", "chr3"}, {"1", "2"});
    ASSERT_TRUE(!err.empty());
}

void test_validate_chrom_naming_reverse_mismatch() {
    std::string err = validate_chrom_naming({"1", "2", "3"}, {"chr1", "chr2"});
    ASSERT_TRUE(!err.empty());
}

void test_validate_chrom_naming_unknown_chrom() {
    std::string err = validate_chrom_naming({"chr1", "chr2"}, {"chr1", "chrZ"});
    ASSERT_TRUE(!err.empty());
}

int main() {
    test_load_basic();
    test_load_columns_out_of_order();
    test_load_missing_column_throws();
    test_load_malformed_number_throws();
    test_find_lookup();
    test_compute_bin_coordinate();
    test_validate_chrom_naming_ok();
    test_validate_chrom_naming_chr_prefix_mismatch();
    test_validate_chrom_naming_reverse_mismatch();
    test_validate_chrom_naming_unknown_chrom();
    TEST_REPORT();
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++17 -Isrc -c tests/test_bins.cpp -o /tmp/test_bins.o`
Expected: FAIL — `bins.hpp` does not exist yet.

- [ ] **Step 3: Write `src/bins.hpp`**

```cpp
#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <optional>
#include "cli_args.hpp"

struct Bin {
    std::string chrom;
    int64_t start;
    int64_t end;
    double gc;
};

class BinIndex {
public:
    static BinIndex load(const std::string& path);

    size_t size() const { return bins_.size(); }
    std::string name(size_t i) const;
    std::optional<size_t> find(const std::string& chrom, int64_t coord) const;
    const std::vector<Bin>& bins() const { return bins_; }

private:
    std::vector<Bin> bins_;
    std::unordered_map<std::string, std::vector<std::pair<int64_t, size_t>>> by_chrom_;
};

int64_t compute_bin_coordinate(int64_t align_start, int64_t align_end, PositionMode mode);

std::string validate_chrom_naming(const std::vector<std::string>& bam_contigs,
                                   const std::vector<std::string>& bin_chroms);
```

- [ ] **Step 4: Write `src/bins.cpp`**

```cpp
#include "bins.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <unordered_set>

namespace {

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string field;
    while (std::getline(ss, field, '\t')) {
        fields.push_back(field);
    }
    return fields;
}

int64_t parse_i64(const std::string& s, const std::string& context) {
    try {
        size_t consumed = 0;
        int64_t v = std::stoll(s, &consumed);
        if (consumed != s.size()) throw std::invalid_argument("trailing characters");
        return v;
    } catch (const std::exception&) {
        throw std::runtime_error("bin file: expected an integer for " + context +
                                  ", got '" + s + "'");
    }
}

double parse_f64(const std::string& s, const std::string& context) {
    try {
        size_t consumed = 0;
        double v = std::stod(s, &consumed);
        if (consumed != s.size()) throw std::invalid_argument("trailing characters");
        return v;
    } catch (const std::exception&) {
        throw std::runtime_error("bin file: expected a number for " + context +
                                  ", got '" + s + "'");
    }
}

} // namespace

BinIndex BinIndex::load(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("cannot open bin file: " + path);
    }

    std::string header_line;
    if (!std::getline(in, header_line)) {
        throw std::runtime_error("bin file is empty: " + path);
    }
    std::vector<std::string> header = split_tab(header_line);

    auto col_index = [&](const std::string& name) -> int {
        for (size_t i = 0; i < header.size(); ++i) {
            if (header[i] == name) return static_cast<int>(i);
        }
        return -1;
    };

    int i_chrom = col_index("chrom");
    int i_start = col_index("start");
    int i_end = col_index("end");
    int i_gc = col_index("gc");

    std::vector<std::string> missing;
    if (i_chrom < 0) missing.push_back("chrom");
    if (i_start < 0) missing.push_back("start");
    if (i_end < 0) missing.push_back("end");
    if (i_gc < 0) missing.push_back("gc");
    if (!missing.empty()) {
        std::string msg = "bin file missing required column(s):";
        for (const auto& m : missing) msg += " " + m;
        throw std::runtime_error(msg + " (" + path + ")");
    }

    BinIndex idx;
    std::string line;
    size_t line_no = 1;
    while (std::getline(in, line)) {
        ++line_no;
        if (line.empty()) continue;
        std::vector<std::string> fields = split_tab(line);
        size_t max_needed = static_cast<size_t>(
            std::max({i_chrom, i_start, i_end, i_gc})) + 1;
        if (fields.size() < max_needed) {
            throw std::runtime_error("bin file: line " + std::to_string(line_no) +
                                      " has too few columns (" + path + ")");
        }
        Bin bin;
        bin.chrom = fields[i_chrom];
        std::string ctx = "line " + std::to_string(line_no);
        bin.start = parse_i64(fields[i_start], "start at " + ctx);
        bin.end = parse_i64(fields[i_end], "end at " + ctx);
        bin.gc = parse_f64(fields[i_gc], "gc at " + ctx);

        size_t original_index = idx.bins_.size();
        idx.bins_.push_back(bin);
        idx.by_chrom_[bin.chrom].emplace_back(bin.start, original_index);
    }

    for (auto& [chrom, vec] : idx.by_chrom_) {
        std::sort(vec.begin(), vec.end());
    }

    return idx;
}

std::string BinIndex::name(size_t i) const {
    const Bin& b = bins_[i];
    return b.chrom + ":" + std::to_string(b.start) + ":" + std::to_string(b.end);
}

std::optional<size_t> BinIndex::find(const std::string& chrom, int64_t coord) const {
    auto it = by_chrom_.find(chrom);
    if (it == by_chrom_.end()) return std::nullopt;

    const auto& vec = it->second;
    // Find the last entry with start <= coord.
    auto upper = std::upper_bound(vec.begin(), vec.end(),
        std::make_pair(coord, static_cast<size_t>(-1)));
    if (upper == vec.begin()) return std::nullopt;
    --upper;
    size_t bin_idx = upper->second;
    const Bin& b = bins_[bin_idx];
    if (coord >= b.start && coord < b.end) {
        return bin_idx;
    }
    return std::nullopt;
}

int64_t compute_bin_coordinate(int64_t align_start, int64_t align_end, PositionMode mode) {
    switch (mode) {
        case PositionMode::Start:
            return align_start;
        case PositionMode::End:
            return align_end > align_start ? align_end - 1 : align_start;
        case PositionMode::Midpoint:
            return align_start + (align_end - align_start) / 2;
    }
    return align_start;
}

std::string validate_chrom_naming(const std::vector<std::string>& bam_contigs,
                                   const std::vector<std::string>& bin_chroms) {
    std::unordered_set<std::string> bam_set(bam_contigs.begin(), bam_contigs.end());
    std::vector<std::string> missing;
    for (const auto& c : bin_chroms) {
        if (bam_set.find(c) == bam_set.end()) missing.push_back(c);
    }
    if (missing.empty()) return "";

    auto has_chr_prefix = [](const std::string& s) {
        return s.size() > 3 && s.compare(0, 3, "chr") == 0;
    };

    size_t bam_chr_count = 0;
    for (const auto& c : bam_contigs) if (has_chr_prefix(c)) bam_chr_count++;
    size_t missing_chr_count = 0;
    for (const auto& c : missing) if (has_chr_prefix(c)) missing_chr_count++;

    bool bam_mostly_chr = !bam_contigs.empty() && bam_chr_count * 2 > bam_contigs.size();
    bool missing_mostly_chr = !missing.empty() && missing_chr_count * 2 > missing.size();

    std::ostringstream oss;
    if (bam_mostly_chr && !missing_mostly_chr) {
        oss << "Chromosome-naming mismatch: BAM header uses 'chr'-prefixed contig "
               "names (e.g. 'chr1'), but " << missing.size()
            << " bin chromosome(s) do not (e.g. '" << missing.front()
            << "'). Fix the bin file or BAM to use a consistent naming convention.";
    } else if (!bam_mostly_chr && missing_mostly_chr) {
        oss << "Chromosome-naming mismatch: bin file uses 'chr'-prefixed chromosome "
               "names (e.g. '" << missing.front()
            << "'), but the BAM header's contig names do not. Fix the bin file or "
               "BAM to use a consistent naming convention.";
    } else {
        oss << "Bin file references " << missing.size()
            << " chromosome(s) not found in the BAM header, e.g. '"
            << missing.front() << "'.";
    }
    return oss.str();
}
```

- [ ] **Step 5: Run test to verify it passes**

Run: `g++ -std=c++17 -Isrc tests/test_bins.cpp src/bins.cpp src/cli_args.cpp -o /tmp/test_bins && /tmp/test_bins`
Expected: `0 failed`.

- [ ] **Step 6: Wire the Makefile**

Append before `clean:`:
```makefile
$(BUILD_DIR)/bins.o: src/bins.cpp src/bins.hpp src/cli_args.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_BUILD)/test_bins: tests/test_bins.cpp $(BUILD_DIR)/bins.o $(BUILD_DIR)/cli_args.o | $(TEST_BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)
```

Update `TESTS :=` line:
```makefile
TESTS := $(TEST_BUILD)/test_cli_args $(TEST_BUILD)/test_bins
```

- [ ] **Step 7: Run full build and test**

Run: `make clean && make && make test`
Expected: builds cleanly, both test binaries report `0 failed`.

- [ ] **Step 8: Commit**

```bash
git add src/bins.hpp src/bins.cpp tests/test_bins.cpp Makefile
git commit -m "Add bin file parsing, bin lookup index, and chromosome-naming validation"
```

---

### Task 4: Barcode index

**Files:**
- Create: `src/barcode_index.hpp`
- Create: `src/barcode_index.cpp`
- Test: `tests/test_barcode_index.cpp`
- Modify: `Makefile`

**Interfaces:**
- Produces: `class BarcodeIndex` with default constructor (no-allowlist mode), `static BarcodeIndex load_allowlist(const std::string& path)` (throws `std::runtime_error` on empty/unreadable file), `bool has_allowlist() const`, `std::optional<size_t> get_or_create(const std::string& barcode)`, `size_t size() const`, `const std::vector<std::string>& names() const`. Task 9 relies on all of these.

- [ ] **Step 1: Write the failing test**

`tests/test_barcode_index.cpp`:
```cpp
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++17 -Isrc -c tests/test_barcode_index.cpp -o /tmp/test_barcode_index.o`
Expected: FAIL — `barcode_index.hpp` does not exist yet.

- [ ] **Step 3: Write `src/barcode_index.hpp`**

```cpp
#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <optional>

class BarcodeIndex {
public:
    BarcodeIndex() = default;

    static BarcodeIndex load_allowlist(const std::string& path);

    bool has_allowlist() const { return has_allowlist_; }
    std::optional<size_t> get_or_create(const std::string& barcode);
    size_t size() const { return names_.size(); }
    const std::vector<std::string>& names() const { return names_; }

private:
    bool has_allowlist_ = false;
    std::vector<std::string> names_;
    std::unordered_map<std::string, size_t> index_;
};
```

- [ ] **Step 4: Write `src/barcode_index.cpp`**

```cpp
#include "barcode_index.hpp"
#include <fstream>
#include <stdexcept>

BarcodeIndex BarcodeIndex::load_allowlist(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("cannot read barcode allowlist file: " + path);
    }

    BarcodeIndex idx;
    idx.has_allowlist_ = true;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        idx.index_[line] = idx.names_.size();
        idx.names_.push_back(line);
    }

    if (idx.names_.empty()) {
        throw std::runtime_error("barcode allowlist file is empty: " + path);
    }

    return idx;
}

std::optional<size_t> BarcodeIndex::get_or_create(const std::string& barcode) {
    auto it = index_.find(barcode);
    if (it != index_.end()) return it->second;

    if (has_allowlist_) return std::nullopt;

    size_t new_index = names_.size();
    index_[barcode] = new_index;
    names_.push_back(barcode);
    return new_index;
}
```

- [ ] **Step 5: Run test to verify it passes**

Run: `g++ -std=c++17 -Isrc tests/test_barcode_index.cpp src/barcode_index.cpp -o /tmp/test_barcode_index && /tmp/test_barcode_index`
Expected: `0 failed`.

- [ ] **Step 6: Wire the Makefile**

Append before `clean:`:
```makefile
$(BUILD_DIR)/barcode_index.o: src/barcode_index.cpp src/barcode_index.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_BUILD)/test_barcode_index: tests/test_barcode_index.cpp $(BUILD_DIR)/barcode_index.o | $(TEST_BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)
```

Update `TESTS :=` line:
```makefile
TESTS := $(TEST_BUILD)/test_cli_args $(TEST_BUILD)/test_bins $(TEST_BUILD)/test_barcode_index
```

- [ ] **Step 7: Run full build and test**

Run: `make clean && make && make test`
Expected: builds cleanly, all three test binaries report `0 failed`.

- [ ] **Step 8: Commit**

```bash
git add src/barcode_index.hpp src/barcode_index.cpp tests/test_barcode_index.cpp Makefile
git commit -m "Add barcode index: allowlist order or first-seen order"
```

---

### Task 5: Counts matrix and gzip TSV writer

**Files:**
- Create: `src/counts_matrix.hpp`
- Create: `src/counts_matrix.cpp`
- Test: `tests/test_counts_matrix.cpp`
- Modify: `Makefile`

**Interfaces:**
- Produces: `class CountsMatrix` with `explicit CountsMatrix(size_t num_bins, size_t initial_barcodes = 0)`, `void grow_to(size_t num_barcodes)`, `void increment(size_t bin_idx, size_t barcode_idx)`, `void write_gz(const std::string& path, const std::vector<std::string>& bin_names, const std::vector<std::string>& barcode_names) const`, `int32_t at(size_t bin_idx, size_t barcode_idx) const`, `size_t num_bins() const`, `size_t num_barcodes() const`. Task 9 relies on all of these.

- [ ] **Step 1: Write the failing test**

`tests/test_counts_matrix.cpp`:
```cpp
#include "test_framework.hpp"
#include "counts_matrix.hpp"
#include <zlib.h>
#include <sstream>

void test_increment_and_read_back() {
    CountsMatrix m(2, 2);
    m.increment(0, 0);
    m.increment(0, 0);
    m.increment(1, 1);
    ASSERT_EQ(m.at(0, 0), 2);
    ASSERT_EQ(m.at(0, 1), 0);
    ASSERT_EQ(m.at(1, 1), 1);
}

void test_grow_preserves_existing_data() {
    CountsMatrix m(2, 1);
    m.increment(0, 0);
    m.increment(1, 0);
    m.grow_to(3);
    ASSERT_EQ(m.num_barcodes(), static_cast<size_t>(3));
    ASSERT_EQ(m.at(0, 0), 1);
    ASSERT_EQ(m.at(1, 0), 1);
    ASSERT_EQ(m.at(0, 2), 0);
    m.increment(0, 2);
    ASSERT_EQ(m.at(0, 2), 1);
}

void test_increment_grows_implicitly() {
    CountsMatrix m(1, 0);
    m.increment(0, 5);
    ASSERT_EQ(m.num_barcodes(), static_cast<size_t>(6));
    ASSERT_EQ(m.at(0, 5), 1);
}

std::string read_gz_file(const std::string& path) {
    gzFile f = gzopen(path.c_str(), "rb");
    std::string content;
    char buf[4096];
    int n;
    while ((n = gzread(f, buf, sizeof(buf))) > 0) {
        content.append(buf, n);
    }
    gzclose(f);
    return content;
}

void test_write_gz_format() {
    CountsMatrix m(2, 2);
    m.increment(0, 0);
    m.increment(0, 0);
    m.increment(1, 1);
    std::string path = "/tmp/test_counts_matrix_out.txt.gz";
    m.write_gz(path, {"chr1:0:1000", "chr1:1000:2000"}, {"AAAA-1", "CCCC-1"});

    std::string content = read_gz_file(path);
    std::string expected =
        "bin\tAAAA-1\tCCCC-1\n"
        "chr1:0:1000\t2\t0\n"
        "chr1:1000:2000\t0\t1\n";
    ASSERT_EQ(content, expected);
}

int main() {
    test_increment_and_read_back();
    test_grow_preserves_existing_data();
    test_increment_grows_implicitly();
    test_write_gz_format();
    TEST_REPORT();
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++17 -Isrc -c tests/test_counts_matrix.cpp -o /tmp/test_counts_matrix.o`
Expected: FAIL — `counts_matrix.hpp` does not exist yet.

- [ ] **Step 3: Write `src/counts_matrix.hpp`**

```cpp
#pragma once
#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>

class CountsMatrix {
public:
    explicit CountsMatrix(size_t num_bins, size_t initial_barcodes = 0);

    void grow_to(size_t num_barcodes);
    void increment(size_t bin_idx, size_t barcode_idx);

    void write_gz(const std::string& path,
                  const std::vector<std::string>& bin_names,
                  const std::vector<std::string>& barcode_names) const;

    int32_t at(size_t bin_idx, size_t barcode_idx) const;
    size_t num_bins() const { return num_bins_; }
    size_t num_barcodes() const { return num_barcodes_; }

private:
    size_t num_bins_;
    size_t num_barcodes_;
    size_t stride_;
    std::vector<int32_t> data_;
};
```

- [ ] **Step 4: Write `src/counts_matrix.cpp`**

```cpp
#include "counts_matrix.hpp"
#include <zlib.h>
#include <sstream>
#include <stdexcept>

CountsMatrix::CountsMatrix(size_t num_bins, size_t initial_barcodes)
    : num_bins_(num_bins),
      num_barcodes_(initial_barcodes),
      stride_(initial_barcodes == 0 ? 1 : initial_barcodes),
      data_(num_bins_ * stride_, 0) {}

void CountsMatrix::grow_to(size_t n) {
    if (n <= num_barcodes_) return;

    if (n > stride_) {
        size_t new_stride = stride_ == 0 ? 1 : stride_;
        while (new_stride < n) new_stride *= 2;

        std::vector<int32_t> new_data(num_bins_ * new_stride, 0);
        for (size_t bin = 0; bin < num_bins_; ++bin) {
            for (size_t bc = 0; bc < num_barcodes_; ++bc) {
                new_data[bin * new_stride + bc] = data_[bin * stride_ + bc];
            }
        }
        data_ = std::move(new_data);
        stride_ = new_stride;
    }
    num_barcodes_ = n;
}

void CountsMatrix::increment(size_t bin_idx, size_t barcode_idx) {
    if (barcode_idx >= num_barcodes_) grow_to(barcode_idx + 1);
    data_[bin_idx * stride_ + barcode_idx]++;
}

int32_t CountsMatrix::at(size_t bin_idx, size_t barcode_idx) const {
    return data_[bin_idx * stride_ + barcode_idx];
}

void CountsMatrix::write_gz(const std::string& path,
                             const std::vector<std::string>& bin_names,
                             const std::vector<std::string>& barcode_names) const {
    gzFile out = gzopen(path.c_str(), "wb");
    if (!out) {
        throw std::runtime_error("cannot open output file for writing: " + path);
    }

    std::ostringstream header;
    header << "bin";
    for (const auto& name : barcode_names) header << '\t' << name;
    header << '\n';
    std::string header_str = header.str();
    gzwrite(out, header_str.data(), static_cast<unsigned>(header_str.size()));

    for (size_t bin = 0; bin < num_bins_; ++bin) {
        std::ostringstream row;
        row << bin_names[bin];
        for (size_t bc = 0; bc < barcode_names.size(); ++bc) {
            row << '\t' << at(bin, bc);
        }
        row << '\n';
        std::string row_str = row.str();
        gzwrite(out, row_str.data(), static_cast<unsigned>(row_str.size()));
    }

    gzclose(out);
}
```

- [ ] **Step 5: Run test to verify it passes**

Run: `g++ -std=c++17 -Isrc tests/test_counts_matrix.cpp src/counts_matrix.cpp -o /tmp/test_counts_matrix $(pkg-config --cflags --libs htslib) -lz && /tmp/test_counts_matrix`
Expected: `0 failed`.

- [ ] **Step 6: Wire the Makefile**

Append before `clean:`:
```makefile
$(BUILD_DIR)/counts_matrix.o: src/counts_matrix.cpp src/counts_matrix.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_BUILD)/test_counts_matrix: tests/test_counts_matrix.cpp $(BUILD_DIR)/counts_matrix.o | $(TEST_BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)
```

Update `TESTS :=` line:
```makefile
TESTS := $(TEST_BUILD)/test_cli_args $(TEST_BUILD)/test_bins $(TEST_BUILD)/test_barcode_index $(TEST_BUILD)/test_counts_matrix
```

- [ ] **Step 7: Run full build and test**

Run: `make clean && make && make test`
Expected: builds cleanly, all four test binaries report `0 failed`.

- [ ] **Step 8: Commit**

```bash
git add src/counts_matrix.hpp src/counts_matrix.cpp tests/test_counts_matrix.cpp Makefile
git commit -m "Add dense growable counts matrix with gzip TSV writer"
```

---

### Task 6: Read filter classification (discordant detection + quality filters)

**Files:**
- Create: `src/read_filter.hpp`
- Create: `src/read_filter.cpp`
- Test: `tests/test_read_filter.cpp`
- Modify: `tests/test_framework.hpp` (add `make_test_record()` helper)
- Modify: `Makefile`

**Interfaces:**
- Consumes: `Args` from `cli_args.hpp`.
- Produces: `enum class DiscordantReason { None, SplitChimeric, CrossChromosome, ImproperPair }`, `bool is_unconditionally_dropped(const bam1_t* rec)`, `DiscordantReason classify_discordant(const bam1_t* rec, const Args& args)`, `bool passes_quality_filters(const bam1_t* rec, const Args& args)`. Task 9 relies on all three functions and the enum.

- [ ] **Step 1: Extend the test helper**

Append to `tests/test_framework.hpp` (after the existing macros, before nothing else needed):
```cpp
#include <htslib/sam.h>
#include <string>
#include <stdexcept>

// Builds a minimal synthetic alignment record for unit testing: a single
// `seq_len`M CIGAR op, no indels/clips. Caller must bam_destroy1() it.
inline bam1_t* make_test_record(const std::string& qname, uint16_t flag,
                                 int32_t tid, int64_t pos, uint8_t mapq,
                                 int32_t mtid, int64_t mpos, int64_t isize,
                                 int32_t seq_len) {
    bam1_t* rec = bam_init1();
    uint32_t cigar[1] = { bam_cigar_gen(seq_len, BAM_CMATCH) };
    std::string seq(seq_len, 'A');
    std::string qual(seq_len, static_cast<char>(30));
    int ret = bam_set1(rec, qname.size(), qname.c_str(), flag, tid, pos, mapq,
                        1, cigar, mtid, mpos, isize,
                        seq.size(), seq.c_str(), qual.c_str(), 32);
    if (ret < 0) {
        bam_destroy1(rec);
        throw std::runtime_error("bam_set1 failed for test record " + qname);
    }
    return rec;
}
```

- [ ] **Step 2: Write the failing test**

`tests/test_read_filter.cpp`:
```cpp
#include "test_framework.hpp"
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
```

- [ ] **Step 3: Run test to verify it fails**

Run: `g++ -std=c++17 -Isrc $(pkg-config --cflags htslib) -c tests/test_read_filter.cpp -o /tmp/test_read_filter.o`
Expected: FAIL — `read_filter.hpp` does not exist yet.

- [ ] **Step 4: Write `src/read_filter.hpp`**

```cpp
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
```

- [ ] **Step 5: Write `src/read_filter.cpp`**

```cpp
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
```

- [ ] **Step 6: Run test to verify it passes**

Run: `g++ -std=c++17 -Isrc $(pkg-config --cflags htslib) tests/test_read_filter.cpp src/read_filter.cpp src/cli_args.cpp -o /tmp/test_read_filter $(pkg-config --libs htslib) -lz && /tmp/test_read_filter`
Expected: `0 failed`.

- [ ] **Step 7: Wire the Makefile**

Append before `clean:`:
```makefile
$(BUILD_DIR)/read_filter.o: src/read_filter.cpp src/read_filter.hpp src/cli_args.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_BUILD)/test_read_filter: tests/test_read_filter.cpp $(BUILD_DIR)/read_filter.o $(BUILD_DIR)/cli_args.o | $(TEST_BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)
```

Update `TESTS :=` line:
```makefile
TESTS := $(TEST_BUILD)/test_cli_args $(TEST_BUILD)/test_bins $(TEST_BUILD)/test_barcode_index $(TEST_BUILD)/test_counts_matrix $(TEST_BUILD)/test_read_filter
```

- [ ] **Step 8: Run full build and test**

Run: `make clean && make && make test`
Expected: builds cleanly, all five test binaries report `0 failed`.

- [ ] **Step 9: Commit**

```bash
git add src/read_filter.hpp src/read_filter.cpp tests/test_read_filter.cpp tests/test_framework.hpp Makefile
git commit -m "Add discordant-read classification and quality filters"
```

---

### Task 7: Fragment-pairing buffer

**Files:**
- Create: `src/fragment_pairing.hpp`
- Create: `src/fragment_pairing.cpp`
- Test: `tests/test_fragment_pairing.cpp`
- Modify: `Makefile`

**Interfaces:**
- Consumes: `SingleEndCountingMode` from `cli_args.hpp`.
- Produces: `struct FragmentCandidate { std::string qname; int32_t tid; int64_t align_start; int64_t align_end; int64_t mate_pos; size_t barcode_idx; bool is_read1; }`, `class FragmentPairingBuffer` with `explicit FragmentPairingBuffer(SingleEndCountingMode mode)`, `std::optional<FragmentCandidate> process(const FragmentCandidate& incoming)`, `std::vector<FragmentCandidate> flush_up_to(int32_t tid, int64_t before_pos)`, `std::vector<FragmentCandidate> flush_all()`, `size_t pending_count() const`. Task 9 relies on all of these exact names and fields.

- [ ] **Step 1: Write the failing test**

`tests/test_fragment_pairing.cpp`:
```cpp
#include "test_framework.hpp"
#include "fragment_pairing.hpp"

FragmentCandidate make_candidate(std::string qname, int32_t tid, int64_t start,
                                  int64_t end, int64_t mate_pos, size_t barcode_idx,
                                  bool is_read1) {
    FragmentCandidate c;
    c.qname = qname;
    c.tid = tid;
    c.align_start = start;
    c.align_end = end;
    c.mate_pos = mate_pos;
    c.barcode_idx = barcode_idx;
    c.is_read1 = is_read1;
    return c;
}

void test_auto_prefers_longer_mate() {
    FragmentPairingBuffer buf(SingleEndCountingMode::Auto);
    auto r1 = make_candidate("q1", 0, 500, 527, 520, 0, true);   // 27bp
    auto r2 = make_candidate("q1", 0, 520, 571, 500, 0, false);  // 51bp

    auto result1 = buf.process(r1);
    ASSERT_TRUE(!result1.has_value());
    ASSERT_EQ(buf.pending_count(), static_cast<size_t>(1));

    auto result2 = buf.process(r2);
    ASSERT_TRUE(result2.has_value());
    ASSERT_EQ(result2->align_start, static_cast<int64_t>(520));
    ASSERT_TRUE(!result2->is_read1);
    ASSERT_EQ(buf.pending_count(), static_cast<size_t>(0));
}

void test_auto_tie_breaks_to_read1() {
    FragmentPairingBuffer buf(SingleEndCountingMode::Auto);
    auto r2 = make_candidate("q1", 0, 500, 550, 600, 0, false);
    auto r1 = make_candidate("q1", 0, 600, 650, 500, 0, true);

    buf.process(r2);
    auto result = buf.process(r1);
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(result->is_read1);
}

void test_r1_mode_picks_fixed_mate() {
    FragmentPairingBuffer buf(SingleEndCountingMode::R1);
    auto r1 = make_candidate("q1", 0, 500, 527, 520, 0, true);
    auto r2 = make_candidate("q1", 0, 520, 571, 500, 0, false);

    buf.process(r1);
    auto result = buf.process(r2);
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(result->is_read1);
    ASSERT_EQ(result->align_start, static_cast<int64_t>(500));
}

void test_r2_mode_picks_fixed_mate_when_both_present() {
    FragmentPairingBuffer buf(SingleEndCountingMode::R2);
    auto r1 = make_candidate("q1", 0, 500, 527, 520, 0, true);
    auto r2 = make_candidate("q1", 0, 520, 571, 500, 0, false);

    buf.process(r1);
    auto result = buf.process(r2);
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(!result->is_read1);
}

void test_auto_singleton_flush_is_counted() {
    FragmentPairingBuffer buf(SingleEndCountingMode::Auto);
    auto r1 = make_candidate("q1", 0, 500, 527, 520, 0, true);
    buf.process(r1);

    auto flushed = buf.flush_up_to(0, 300); // hasn't passed mate_pos=520 yet
    ASSERT_EQ(flushed.size(), static_cast<size_t>(0));
    ASSERT_EQ(buf.pending_count(), static_cast<size_t>(1));

    auto flushed2 = buf.flush_up_to(0, 600); // now past mate_pos=520
    ASSERT_EQ(flushed2.size(), static_cast<size_t>(1));
    ASSERT_EQ(flushed2[0].align_start, static_cast<int64_t>(500));
    ASSERT_EQ(buf.pending_count(), static_cast<size_t>(0));
}

void test_r1_mode_singleton_read2_only_is_dropped() {
    FragmentPairingBuffer buf(SingleEndCountingMode::R1);
    auto r2 = make_candidate("q1", 0, 520, 571, 500, 0, false);
    buf.process(r2);

    auto flushed = buf.flush_up_to(0, 600);
    ASSERT_EQ(flushed.size(), static_cast<size_t>(0));
}

void test_r1_mode_singleton_read1_only_is_kept() {
    FragmentPairingBuffer buf(SingleEndCountingMode::R1);
    auto r1 = make_candidate("q1", 0, 500, 527, 520, 0, true);
    buf.process(r1);

    auto flushed = buf.flush_up_to(0, 600);
    ASSERT_EQ(flushed.size(), static_cast<size_t>(1));
    ASSERT_TRUE(flushed[0].is_read1);
}

void test_flush_all_at_end_of_stream() {
    FragmentPairingBuffer buf(SingleEndCountingMode::Auto);
    buf.process(make_candidate("q1", 0, 500, 527, 520, 0, true));
    buf.process(make_candidate("q2", 1, 100, 150, 200, 1, true));

    auto flushed = buf.flush_all();
    ASSERT_EQ(flushed.size(), static_cast<size_t>(2));
    ASSERT_EQ(buf.pending_count(), static_cast<size_t>(0));
}

int main() {
    test_auto_prefers_longer_mate();
    test_auto_tie_breaks_to_read1();
    test_r1_mode_picks_fixed_mate();
    test_r2_mode_picks_fixed_mate_when_both_present();
    test_auto_singleton_flush_is_counted();
    test_r1_mode_singleton_read2_only_is_dropped();
    test_r1_mode_singleton_read1_only_is_kept();
    test_flush_all_at_end_of_stream();
    TEST_REPORT();
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++17 -Isrc -c tests/test_fragment_pairing.cpp -o /tmp/test_fragment_pairing.o`
Expected: FAIL — `fragment_pairing.hpp` does not exist yet.

- [ ] **Step 3: Write `src/fragment_pairing.hpp`**

```cpp
#pragma once
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>
#include <optional>
#include "cli_args.hpp"

struct FragmentCandidate {
    std::string qname;
    int32_t tid = -1;
    int64_t align_start = 0;
    int64_t align_end = 0;
    int64_t mate_pos = -1;
    size_t barcode_idx = 0;
    bool is_read1 = true;
};

class FragmentPairingBuffer {
public:
    explicit FragmentPairingBuffer(SingleEndCountingMode mode);

    std::optional<FragmentCandidate> process(const FragmentCandidate& incoming);
    std::vector<FragmentCandidate> flush_up_to(int32_t tid, int64_t before_pos);
    std::vector<FragmentCandidate> flush_all();

    size_t pending_count() const { return pending_.size(); }

private:
    SingleEndCountingMode mode_;
    std::unordered_map<std::string, FragmentCandidate> pending_;

    std::optional<FragmentCandidate> resolve_pair(const FragmentCandidate& a,
                                                   const FragmentCandidate& b) const;
    std::optional<FragmentCandidate> resolve_singleton(const FragmentCandidate& c) const;
};
```

- [ ] **Step 4: Write `src/fragment_pairing.cpp`**

```cpp
#include "fragment_pairing.hpp"

FragmentPairingBuffer::FragmentPairingBuffer(SingleEndCountingMode mode) : mode_(mode) {}

std::optional<FragmentCandidate> FragmentPairingBuffer::process(
        const FragmentCandidate& incoming) {
    auto it = pending_.find(incoming.qname);
    if (it == pending_.end()) {
        pending_.emplace(incoming.qname, incoming);
        return std::nullopt;
    }
    FragmentCandidate other = it->second;
    pending_.erase(it);
    return resolve_pair(incoming, other);
}

std::optional<FragmentCandidate> FragmentPairingBuffer::resolve_pair(
        const FragmentCandidate& a, const FragmentCandidate& b) const {
    switch (mode_) {
        case SingleEndCountingMode::Auto: {
            int64_t len_a = a.align_end - a.align_start;
            int64_t len_b = b.align_end - b.align_start;
            if (len_a != len_b) return (len_a > len_b) ? a : b;
            return a.is_read1 ? a : b;
        }
        case SingleEndCountingMode::R1:
            if (a.is_read1) return a;
            if (b.is_read1) return b;
            return std::nullopt;
        case SingleEndCountingMode::R2:
            if (!a.is_read1) return a;
            if (!b.is_read1) return b;
            return std::nullopt;
    }
    return std::nullopt;
}

std::optional<FragmentCandidate> FragmentPairingBuffer::resolve_singleton(
        const FragmentCandidate& c) const {
    switch (mode_) {
        case SingleEndCountingMode::Auto:
            return c;
        case SingleEndCountingMode::R1:
            return c.is_read1 ? std::optional<FragmentCandidate>(c) : std::nullopt;
        case SingleEndCountingMode::R2:
            return !c.is_read1 ? std::optional<FragmentCandidate>(c) : std::nullopt;
    }
    return std::nullopt;
}

std::vector<FragmentCandidate> FragmentPairingBuffer::flush_up_to(
        int32_t tid, int64_t before_pos) {
    std::vector<FragmentCandidate> results;
    for (auto it = pending_.begin(); it != pending_.end(); ) {
        const FragmentCandidate& c = it->second;
        if (c.tid == tid && c.mate_pos < before_pos) {
            auto resolved = resolve_singleton(c);
            if (resolved) results.push_back(*resolved);
            it = pending_.erase(it);
        } else {
            ++it;
        }
    }
    return results;
}

std::vector<FragmentCandidate> FragmentPairingBuffer::flush_all() {
    std::vector<FragmentCandidate> results;
    for (auto& entry : pending_) {
        auto resolved = resolve_singleton(entry.second);
        if (resolved) results.push_back(*resolved);
    }
    pending_.clear();
    return results;
}
```

- [ ] **Step 5: Run test to verify it passes**

Run: `g++ -std=c++17 -Isrc tests/test_fragment_pairing.cpp src/fragment_pairing.cpp src/cli_args.cpp -o /tmp/test_fragment_pairing && /tmp/test_fragment_pairing`
Expected: `0 failed`.

- [ ] **Step 6: Wire the Makefile**

Append before `clean:`:
```makefile
$(BUILD_DIR)/fragment_pairing.o: src/fragment_pairing.cpp src/fragment_pairing.hpp src/cli_args.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_BUILD)/test_fragment_pairing: tests/test_fragment_pairing.cpp $(BUILD_DIR)/fragment_pairing.o $(BUILD_DIR)/cli_args.o | $(TEST_BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)
```

Update `TESTS :=` line:
```makefile
TESTS := $(TEST_BUILD)/test_cli_args $(TEST_BUILD)/test_bins $(TEST_BUILD)/test_barcode_index $(TEST_BUILD)/test_counts_matrix $(TEST_BUILD)/test_read_filter $(TEST_BUILD)/test_fragment_pairing
```

- [ ] **Step 7: Run full build and test**

Run: `make clean && make && make test`
Expected: builds cleanly, all six test binaries report `0 failed`.

- [ ] **Step 8: Commit**

```bash
git add src/fragment_pairing.hpp src/fragment_pairing.cpp tests/test_fragment_pairing.cpp Makefile
git commit -m "Add fragment-pairing buffer: auto/r1/r2 mate resolution and singleton flushing"
```

---

### Task 8: Discordant BAM writer

**Files:**
- Create: `src/discordant_writer.hpp`
- Create: `src/discordant_writer.cpp`
- Test: `tests/test_discordant_writer.cpp`
- Modify: `Makefile`

**Interfaces:**
- Produces: `class DiscordantWriter` with `DiscordantWriter(const std::string& path, sam_hdr_t* hdr)` (throws `std::runtime_error` on failure), `void write(const bam1_t* rec)`, `void close()`, destructor. Task 9 relies on this exact interface.

- [ ] **Step 1: Write the failing test**

`tests/test_discordant_writer.cpp`:
```cpp
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++17 -Isrc $(pkg-config --cflags htslib) -c tests/test_discordant_writer.cpp -o /tmp/test_discordant_writer.o`
Expected: FAIL — `discordant_writer.hpp` does not exist yet.

- [ ] **Step 3: Write `src/discordant_writer.hpp`**

```cpp
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
```

- [ ] **Step 4: Write `src/discordant_writer.cpp`**

```cpp
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
```

- [ ] **Step 5: Run test to verify it passes**

Run: `g++ -std=c++17 -Isrc $(pkg-config --cflags htslib) tests/test_discordant_writer.cpp src/discordant_writer.cpp -o /tmp/test_discordant_writer $(pkg-config --libs htslib) -lz && /tmp/test_discordant_writer`
Expected: `0 failed`.

- [ ] **Step 6: Wire the Makefile**

Append before `clean:`:
```makefile
$(BUILD_DIR)/discordant_writer.o: src/discordant_writer.cpp src/discordant_writer.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_BUILD)/test_discordant_writer: tests/test_discordant_writer.cpp $(BUILD_DIR)/discordant_writer.o | $(TEST_BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)
```

Update `TESTS :=` line:
```makefile
TESTS := $(TEST_BUILD)/test_cli_args $(TEST_BUILD)/test_bins $(TEST_BUILD)/test_barcode_index $(TEST_BUILD)/test_counts_matrix $(TEST_BUILD)/test_read_filter $(TEST_BUILD)/test_fragment_pairing $(TEST_BUILD)/test_discordant_writer
```

- [ ] **Step 7: Run full build and test**

Run: `make clean && make && make test`
Expected: builds cleanly, all seven test binaries report `0 failed`.

- [ ] **Step 8: Commit**

```bash
git add src/discordant_writer.hpp src/discordant_writer.cpp tests/test_discordant_writer.cpp Makefile
git commit -m "Add discordant BAM writer with automatic indexing"
```

---

### Task 9: Main orchestration + minimal end-to-end smoke test

**Files:**
- Modify: `src/main.cpp` (full rewrite: wires all modules into the streaming pipeline)
- Modify: `Makefile` (link all modules into `$(BIN)`)
- Create: `tests/fixtures/gen_smoke_fixture.py`
- Create: `tests/integration/test_smoke.py`

**Interfaces:**
- Consumes: every module from Tasks 2–8 by their exact names.
- Produces: a working `dfn_bin` binary implementing the full pipeline described in the spec's Algorithm section, tested end-to-end for the "clean reads" path (discordant/pairing-edge-case coverage is Task 10's job).

- [ ] **Step 1: Set up the Python virtualenv**

Run:
```bash
python3 -m venv .venv
.venv/bin/pip install --upgrade pip
.venv/bin/pip install pysam pytest
```
Expected: installs succeed; `.venv/bin/python -c "import pysam; import pytest"` runs with no output/errors.

- [ ] **Step 2: Write the smoke-test fixture generator**

`tests/fixtures/gen_smoke_fixture.py`:
```python
#!/usr/bin/env python3
"""Generates a tiny BAM + bins.tsv covering only the 'clean read' path,
for Task 9's end-to-end smoke test. Comprehensive edge-case fixtures
(discordant reads, fragment-pairing modes, error paths) are in
gen_fixtures.py, added in Task 10."""
import pysam
import sys

def main(out_dir):
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 10000}, {"SN": "chr2", "LN": 10000}],
    }

    bam_path = f"{out_dir}/smoke.bam"
    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        def add_read(qname, tid, pos, cb, mapq=60, seq_len=50):
            a = pysam.AlignedSegment()
            a.query_name = qname
            a.query_sequence = "A" * seq_len
            a.flag = 0
            a.reference_id = tid
            a.reference_start = pos
            a.mapping_quality = mapq
            a.cigarstring = f"{seq_len}M"
            a.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
            a.tags = [("CB", cb)]
            out.write(a)

        add_read("clean1", 0, 150, "AAAA-1")
        add_read("clean2", 0, 1200, "AAAA-1")
        add_read("clean3", 0, 1300, "CCCC-1")
        add_read("clean4", 1, 150, "CCCC-1")

    pysam.index(bam_path)

    bins_path = f"{out_dir}/smoke_bins.tsv"
    with open(bins_path, "w") as f:
        f.write("chrom\tstart\tend\tgc\n")
        f.write("chr1\t0\t1000\t0.40\n")
        f.write("chr1\t1000\t2000\t0.50\n")
        f.write("chr2\t0\t1000\t0.45\n")

    return bam_path, bins_path

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
```

- [ ] **Step 3: Run it and confirm fixtures are generated**

Run: `mkdir -p /tmp/dfn_smoke && .venv/bin/python tests/fixtures/gen_smoke_fixture.py /tmp/dfn_smoke && ls /tmp/dfn_smoke`
Expected: lists `smoke.bam`, `smoke.bam.bai`, `smoke_bins.tsv`.

- [ ] **Step 4: Write the failing integration test**

`tests/integration/test_smoke.py`:
```python
import gzip
import subprocess
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "fixtures"))
import gen_smoke_fixture

REPO_ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
DFN_BIN = os.path.join(REPO_ROOT, "dfn_bin")

def read_matrix(path):
    with gzip.open(path, "rt") as f:
        lines = [line.rstrip("\n").split("\t") for line in f]
    header, rows = lines[0], lines[1:]
    return header, {row[0]: row[1:] for row in rows}

def test_clean_reads_produce_expected_counts(tmp_path):
    bam_path, bins_path = gen_smoke_fixture.main(str(tmp_path))
    out_prefix = str(tmp_path / "out")

    result = subprocess.run(
        [DFN_BIN, "--bam", bam_path, "--bins", bins_path, "--out-prefix", out_prefix],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr

    header, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    barcodes = header[1:]
    assert set(barcodes) == {"AAAA-1", "CCCC-1"}

    aaaa_idx = barcodes.index("AAAA-1")
    cccc_idx = barcodes.index("CCCC-1")

    assert int(rows["chr1:0:1000"][aaaa_idx]) == 1   # clean1
    assert int(rows["chr1:1000:2000"][aaaa_idx]) == 1  # clean2
    assert int(rows["chr1:1000:2000"][cccc_idx]) == 1  # clean3
    assert int(rows["chr2:0:1000"][cccc_idx]) == 1     # clean4
```

- [ ] **Step 5: Run test to verify it fails**

Run: `.venv/bin/pytest tests/integration/test_smoke.py -v`
Expected: FAIL — `dfn_bin` binary doesn't implement the pipeline yet (either missing outputs or all-zero counts).

- [ ] **Step 6: Rewrite `src/main.cpp` with the full pipeline**

```cpp
#include "cli_args.hpp"
#include "bins.hpp"
#include "barcode_index.hpp"
#include "counts_matrix.hpp"
#include "read_filter.hpp"
#include "fragment_pairing.hpp"
#include "discordant_writer.hpp"

#include <htslib/sam.h>
#include <cstdio>
#include <vector>
#include <string>

namespace {

void bin_candidate(const FragmentCandidate& c, const BinIndex& bin_index,
                    const std::vector<std::string>& tid_to_chrom,
                    const Args& args, CountsMatrix& matrix) {
    int64_t coord = compute_bin_coordinate(c.align_start, c.align_end, args.position);
    auto bin_idx = bin_index.find(tid_to_chrom[c.tid], coord);
    if (bin_idx) {
        matrix.increment(*bin_idx, c.barcode_idx);
    }
}

int run(const Args& args) {
    // --- Validate inputs up front ---
    samFile* bam_fp = sam_open(args.bam_path.c_str(), "r");
    if (!bam_fp) {
        std::fprintf(stderr, "error: cannot open BAM file: %s\n", args.bam_path.c_str());
        return 1;
    }
    sam_hdr_t* hdr = sam_hdr_read(bam_fp);
    if (!hdr) {
        std::fprintf(stderr, "error: cannot read BAM header: %s\n", args.bam_path.c_str());
        return 1;
    }

    hts_idx_t* idx = sam_index_load(bam_fp, args.bam_path.c_str());
    if (!idx) {
        std::fprintf(stderr,
            "error: no BAM index found for %s (run `samtools index %s`)\n",
            args.bam_path.c_str(), args.bam_path.c_str());
        return 1;
    }
    hts_idx_destroy(idx);

    BinIndex bin_index;
    try {
        bin_index = BinIndex::load(args.bins_path);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }

    std::vector<std::string> bam_contigs;
    std::vector<std::string> tid_to_chrom(sam_hdr_nref(hdr));
    for (int i = 0; i < sam_hdr_nref(hdr); ++i) {
        std::string name = sam_hdr_tid2name(hdr, i);
        bam_contigs.push_back(name);
        tid_to_chrom[i] = name;
    }
    std::vector<std::string> bin_chroms;
    for (const auto& b : bin_index.bins()) bin_chroms.push_back(b.chrom);
    std::string naming_error = validate_chrom_naming(bam_contigs, bin_chroms);
    if (!naming_error.empty()) {
        std::fprintf(stderr, "error: %s\n", naming_error.c_str());
        return 1;
    }

    BarcodeIndex barcode_index;
    if (!args.barcodes_path.empty()) {
        try {
            barcode_index = BarcodeIndex::load_allowlist(args.barcodes_path);
        } catch (const std::exception& e) {
            std::fprintf(stderr, "error: %s\n", e.what());
            return 1;
        }
    }

    // --- Set up outputs ---
    std::string discordant_path = args.out_prefix + ".discordant.bam";
    DiscordantWriter discordant_writer(discordant_path, hdr);

    CountsMatrix matrix(bin_index.size(), barcode_index.size());

    if (args.threads > 1) {
        hts_set_threads(bam_fp, args.threads);
    }

    // --- Stream the BAM once ---
    FragmentPairingBuffer pairing_buffer(args.single_end_counting);
    bam1_t* rec = bam_init1();

    while (sam_read1(bam_fp, hdr, rec) >= 0) {
        if (is_unconditionally_dropped(rec)) continue;

        DiscordantReason reason = classify_discordant(rec, args);
        if (reason != DiscordantReason::None) {
            discordant_writer.write(rec);
            for (const auto& flushed : pairing_buffer.flush_up_to(rec->core.tid, rec->core.pos)) {
                bin_candidate(flushed, bin_index, tid_to_chrom, args, matrix);
            }
            continue;
        }

        for (const auto& flushed : pairing_buffer.flush_up_to(rec->core.tid, rec->core.pos)) {
            bin_candidate(flushed, bin_index, tid_to_chrom, args, matrix);
        }

        if (!passes_quality_filters(rec, args)) continue;

        uint8_t* bc_aux = bam_aux_get(rec, args.barcode_tag.c_str());
        if (!bc_aux) continue;
        std::string barcode = bam_aux2Z(bc_aux);

        auto barcode_idx = barcode_index.get_or_create(barcode);
        if (!barcode_idx) continue;

        int64_t align_start = rec->core.pos;
        int64_t align_end = bam_endpos(rec);

        bool treat_as_unpaired =
            !(rec->core.flag & BAM_FPAIRED) || (rec->core.flag & BAM_FMUNMAP);

        if (treat_as_unpaired) {
            int64_t coord = compute_bin_coordinate(align_start, align_end, args.position);
            auto bin_idx = bin_index.find(tid_to_chrom[rec->core.tid], coord);
            if (bin_idx) matrix.increment(*bin_idx, *barcode_idx);
            continue;
        }

        FragmentCandidate candidate;
        candidate.qname = bam_get_qname(rec);
        candidate.tid = rec->core.tid;
        candidate.align_start = align_start;
        candidate.align_end = align_end;
        candidate.mate_pos = rec->core.mpos;
        candidate.barcode_idx = *barcode_idx;
        candidate.is_read1 = (rec->core.flag & BAM_FREAD1) != 0;

        auto resolved = pairing_buffer.process(candidate);
        if (resolved) {
            bin_candidate(*resolved, bin_index, tid_to_chrom, args, matrix);
        }
    }

    for (const auto& flushed : pairing_buffer.flush_all()) {
        bin_candidate(flushed, bin_index, tid_to_chrom, args, matrix);
    }

    bam_destroy1(rec);
    discordant_writer.close();
    sam_hdr_destroy(hdr);
    sam_close(bam_fp);

    // --- Write output matrix ---
    std::vector<std::string> bin_names;
    for (size_t i = 0; i < bin_index.size(); ++i) bin_names.push_back(bin_index.name(i));

    matrix.write_gz(args.out_prefix + ".raw_counts.txt.gz", bin_names, barcode_index.names());

    return 0;
}

} // namespace

int main(int argc, char** argv) {
    Args args;
    bool help_requested = false;

    if (!parse_args(argc, argv, args, help_requested)) {
        return 1;
    }
    if (help_requested) {
        std::fputs(usage_text().c_str(), stdout);
        return 0;
    }

    return run(args);
}
```

- [ ] **Step 7: Wire the Makefile to link all modules**

Replace the `$(BIN):` rule:
```makefile
$(BIN): $(BUILD_DIR)/main.o $(BUILD_DIR)/cli_args.o $(BUILD_DIR)/bins.o \
        $(BUILD_DIR)/barcode_index.o $(BUILD_DIR)/counts_matrix.o \
        $(BUILD_DIR)/read_filter.o $(BUILD_DIR)/fragment_pairing.o \
        $(BUILD_DIR)/discordant_writer.o
	$(CXX) $^ -o $@ $(LDLIBS)
```

- [ ] **Step 8: Build and run the smoke test**

Run: `make clean && make`
Expected: builds cleanly, produces `./dfn_bin`.

Run: `.venv/bin/pytest tests/integration/test_smoke.py -v`
Expected: `1 passed`.

- [ ] **Step 9: Run the full C++ unit test suite too, to make sure nothing regressed**

Run: `make test`
Expected: all seven test binaries report `0 failed`.

- [ ] **Step 10: Commit**

```bash
git add src/main.cpp Makefile .venv .gitignore tests/fixtures/gen_smoke_fixture.py tests/integration/test_smoke.py
git status  # confirm .venv is excluded by .gitignore, not staged
git add -A -- ':!.venv'
git commit -m "Wire full streaming pipeline in main.cpp; add end-to-end smoke test"
```

---

### Task 10: Comprehensive fixtures and full integration test suite

**Files:**
- Create: `tests/fixtures/gen_fixtures.py`
- Create: `tests/integration/test_dfn_bin.py`

**Interfaces:**
- Consumes: the `dfn_bin` binary built in Task 9, via subprocess.
- Produces: full coverage of every scenario in the spec's Testing section: discordant detection (all three criteria + `--max-insert-size` override), fragment-pairing modes (`auto`/`r1`/`r2`, unequal-length mates, singleton handling), barcode allowlist filtering, `--exclude-dups` toggle, and all four error-handling paths.

**Fixture design note:** every record is written in true ascending `(tid, pos)` order (unmapped reads last), matching a real coordinate-sorted BAM — the tool's fragment-pairing buffer depends on that ordering, so a fixture that isn't actually sorted would produce misleading test results even though `pysam.index()` wouldn't itself complain. Each fragment-pairing/discordant scenario below gets its own dedicated, well-separated bin so a test can assert an exact count without needing to account for interactions between scenarios sharing a bin.

- [ ] **Step 1: Write the comprehensive fixture generator**

`tests/fixtures/gen_fixtures.py`:
```python
#!/usr/bin/env python3
"""Builds the full fixture set for Task 10's integration tests: one BAM
covering every filtering/discordant/pairing scenario, plus supporting
bins/allowlist/malformed files for the error-handling tests.

Records are listed below in the exact order they're written, which must be
true ascending (tid, pos) order (unmapped last) to match a real
coordinate-sorted BAM -- the tool's fragment-pairing buffer assumes this.
"""
import pysam

CHR1_LEN = 400000
CHR2_LEN = 400000

# SAM FLAG bits used below (see design spec's CLI contract / discordant section):
PAIRED = 1
PROPER_PAIR = 2
UNMAP = 4
MUNMAP = 8
REVERSE = 16
MREVERSE = 32
READ1 = 64
READ2 = 128
SECONDARY = 256
QCFAIL = 512
DUP = 1024
SUPPLEMENTARY = 2048


def _add(out, qname, flag, tid, pos, mapq, mtid, mpos, isize, cb=None,
         seq_len=50, sa=None):
    a = pysam.AlignedSegment()
    a.query_name = qname
    a.query_sequence = "A" * seq_len
    a.flag = flag
    a.reference_id = tid
    a.reference_start = pos
    a.mapping_quality = mapq
    a.cigarstring = f"{seq_len}M"
    a.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
    a.next_reference_id = mtid
    a.next_reference_start = mpos
    a.template_length = isize
    tags = []
    if cb is not None:
        tags.append(("CB", cb))
    if sa is not None:
        tags.append(("SA", sa))
    a.tags = tags
    out.write(a)


def build_main_bam(out_dir):
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": CHR1_LEN}, {"SN": "chr2", "LN": CHR2_LEN}],
    }
    bam_path = f"{out_dir}/fixtures.bam"

    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        # ===== chr1 (tid=0), strictly ascending pos =====

        # --- bin chr1:0-1000: clean1 plus decoys that must NOT inflate it ---
        _add(out, "clean1", 0, 0, 150, 60, -1, -1, 0, cb="AAAA-1")
        _add(out, "qcfail1", QCFAIL, 0, 160, 60, -1, -1, 0, cb="AAAA-1")
        _add(out, "lowmapq1", 0, 0, 170, 5, -1, -1, 0, cb="AAAA-1")
        _add(out, "secondary1", SECONDARY, 0, 180, 60, -1, -1, 0, cb="AAAA-1")
        _add(out, "dup1", DUP, 0, 190, 60, -1, -1, 0, cb="AAAA-1")
        _add(out, "missingcb1", 0, 0, 200, 60, -1, -1, 0, cb=None)

        # --- bin chr1:1000-2000: two different barcodes ---
        _add(out, "clean2", 0, 0, 1200, 60, -1, -1, 0, cb="AAAA-1")
        _add(out, "clean3", 0, 0, 1300, 60, -1, -1, 0, cb="BBBB-1")

        # --- bin chr1:2000-3000: proper pair, small insert -> counted ONCE ---
        _add(out, "properpair_normal", PAIRED | PROPER_PAIR | MREVERSE | READ1,
             0, 2400, 60, 0, 2600, 250, cb="AAAA-1")
        _add(out, "properpair_normal", PAIRED | PROPER_PAIR | REVERSE | READ2,
             0, 2600, 60, 0, 2400, -250, cb="AAAA-1")

        # --- autolen1: 27bp read1 in bin chr1:3000-4000, 51bp read2 in bin
        # chr1:4000-5000 -- separate bins so "which mate won" is directly
        # observable by which bin gets the count. ---
        _add(out, "autolen1", PAIRED | PROPER_PAIR | MREVERSE | READ1,
             0, 3100, 60, 0, 4100, 1050, cb="AAAA-1", seq_len=27)

        # --- improper1: interleaved here, diverted regardless of position ---
        _add(out, "improper1", PAIRED | READ1, 0, 3500, 60, 0, 3600, 150, cb="AAAA-1")
        _add(out, "improper1", PAIRED | READ2, 0, 3600, 60, 0, 3500, -150, cb="AAAA-1")

        _add(out, "autolen1", PAIRED | PROPER_PAIR | REVERSE | READ2,
             0, 4100, 60, 0, 3100, -1050, cb="AAAA-1", seq_len=51)

        # --- bin chr1:5000-6000: singleton1 (read2 filtered by low MAPQ
        # before ever reaching the pairing buffer) ---
        _add(out, "singleton1", PAIRED | PROPER_PAIR | MREVERSE | READ1,
             0, 5100, 60, 0, 5200, 150, cb="AAAA-1")
        _add(out, "singleton1", PAIRED | PROPER_PAIR | REVERSE | READ2,
             0, 5200, 5, 0, 5100, -150, cb="AAAA-1")

        # --- bin chr1:6000-7000: unpaired read whose arrival forces the
        # pairing buffer to flush singleton1 (mate_pos=5200 < 6100) ---
        _add(out, "flush_trigger1", 0, 0, 6100, 60, -1, -1, 0, cb="AAAA-1")

        # --- properpair_bigisize: proper-pair flag set, huge isize (only
        # the isize tag is huge -- actual positions stay in separate,
        # observable bins chr1:7000-8000 / chr1:8000-9000). Default run
        # trusts the proper-pair flag; --max-insert-size overrides it. ---
        _add(out, "properpair_bigisize", PAIRED | PROPER_PAIR | MREVERSE | READ1,
             0, 7100, 60, 0, 8100, 290000, cb="AAAA-1")

        # --- crosschrom1's chr1-side mate, diverted regardless of position ---
        _add(out, "crosschrom1", PAIRED | MREVERSE | READ1, 0, 7500, 60, 1, 500, 0, cb="AAAA-1")

        _add(out, "properpair_bigisize", PAIRED | PROPER_PAIR | REVERSE | READ2,
             0, 8100, 60, 0, 7100, -290000, cb="AAAA-1")

        # --- bin chr1:9000-10000: barcode not on the allowlist (when one is used) ---
        _add(out, "notallowed1", 0, 0, 9100, 60, -1, -1, 0, cb="ZZZZ-1")

        # --- split/chimeric, diverted regardless of position ---
        _add(out, "supplonoSA1", SUPPLEMENTARY, 0, 9500, 60, -1, -1, 0, cb="AAAA-1")
        _add(out, "split1", 0, 0, 9600, 60, -1, -1, 0, cb="AAAA-1",
             sa="chr2,601,+,50M,60,0;")

        # --- off-bin: far past every defined bin ---
        _add(out, "offbin1", 0, 0, 350000, 60, -1, -1, 0, cb="AAAA-1")

        # ===== chr2 (tid=1), strictly ascending pos =====
        _add(out, "clean4", 0, 1, 150, 60, -1, -1, 0, cb="BBBB-1")
        _add(out, "crosschrom1", PAIRED | REVERSE | READ2, 1, 500, 60, 0, 7500, 0, cb="AAAA-1")
        _add(out, "split1", SUPPLEMENTARY, 1, 600, 60, -1, -1, 0, cb="AAAA-1")

        # ===== unmapped, sorted last regardless of chromosome =====
        _add(out, "unmapped1", UNMAP, -1, -1, 0, -1, -1, 0, cb="AAAA-1")

    pysam.index(bam_path)
    return bam_path


def build_bins_file(out_dir):
    path = f"{out_dir}/bins.tsv"
    with open(path, "w") as f:
        f.write("chrom\tstart\tend\tgc\n")
        f.write("chr1\t0\t1000\t0.40\n")
        f.write("chr1\t1000\t2000\t0.50\n")
        f.write("chr1\t2000\t3000\t0.45\n")
        f.write("chr1\t3000\t4000\t0.50\n")
        f.write("chr1\t4000\t5000\t0.55\n")
        f.write("chr1\t5000\t6000\t0.45\n")
        f.write("chr1\t6000\t7000\t0.50\n")
        f.write("chr1\t7000\t8000\t0.45\n")
        f.write("chr1\t8000\t9000\t0.50\n")
        f.write("chr1\t9000\t10000\t0.55\n")
        f.write("chr2\t0\t1000\t0.45\n")
    return path


def build_mismatched_bins_file(out_dir):
    # Uses "1" instead of "chr1" to trigger the naming-mismatch error.
    path = f"{out_dir}/bins_mismatched.tsv"
    with open(path, "w") as f:
        f.write("chrom\tstart\tend\tgc\n")
        f.write("1\t0\t1000\t0.40\n")
    return path


def build_malformed_bins_file(out_dir):
    path = f"{out_dir}/bins_malformed.tsv"
    with open(path, "w") as f:
        f.write("chrom\tstart\tend\n")  # missing 'gc' column
        f.write("chr1\t0\t1000\n")
    return path


def build_allowlist_file(out_dir):
    path = f"{out_dir}/allowlist.txt"
    with open(path, "w") as f:
        f.write("AAAA-1\n")
        f.write("BBBB-1\n")
    return path


def build_empty_allowlist_file(out_dir):
    path = f"{out_dir}/allowlist_empty.txt"
    open(path, "w").close()
    return path


def build_all(out_dir):
    return {
        "bam": build_main_bam(out_dir),
        "bins": build_bins_file(out_dir),
        "bins_mismatched": build_mismatched_bins_file(out_dir),
        "bins_malformed": build_malformed_bins_file(out_dir),
        "allowlist": build_allowlist_file(out_dir),
        "allowlist_empty": build_empty_allowlist_file(out_dir),
    }
```

- [ ] **Step 2: Generate and sanity-check the fixture with samtools**

Run:
```bash
mkdir -p /tmp/dfn_fixtures
.venv/bin/python -c "import sys; sys.path.insert(0,'tests/fixtures'); import gen_fixtures; print(gen_fixtures.build_all('/tmp/dfn_fixtures'))"
samtools view /tmp/dfn_fixtures/fixtures.bam | wc -l
samtools view /tmp/dfn_fixtures/fixtures.bam | cut -f4 | uniq -c   # eyeball ascending order
```
Expected: prints the dict of file paths; `wc -l` reports 28 records (6 in bin0 + 2 clean + 2 properpair_normal + 1 autolen_r1 + 2 improper + 1 autolen_r2 + 2 singleton1 + 1 flush_trigger + 1 properpair_bigisize_r1 + 1 crosschrom1_r1 + 1 properpair_bigisize_r2 + 1 notallowed1 + 1 supplonoSA1 + 1 split1_primary + 1 offbin1 = 24 on chr1, plus clean4 + crosschrom1_r2 + split1_supplementary = 3 on chr2, plus unmapped1 = 1; 24+3+1=28); the `cut -f4` column should be non-decreasing within each chromosome block.

- [ ] **Step 3: Write the failing integration test suite**

`tests/integration/test_dfn_bin.py`:
```python
import gzip
import subprocess
import sys
import os
import pysam
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "fixtures"))
import gen_fixtures

REPO_ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
DFN_BIN = os.path.join(REPO_ROOT, "dfn_bin")


@pytest.fixture
def fixtures(tmp_path):
    return gen_fixtures.build_all(str(tmp_path))


def run_dfn_bin(fixtures, tmp_path, extra_args=None, out_name="out"):
    out_prefix = str(tmp_path / out_name)
    cmd = [DFN_BIN, "--bam", fixtures["bam"], "--bins", fixtures["bins"],
           "--out-prefix", out_prefix]
    if extra_args:
        cmd += extra_args
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result, out_prefix


def read_matrix(path):
    with gzip.open(path, "rt") as f:
        lines = [line.rstrip("\n").split("\t") for line in f]
    header, rows = lines[0], lines[1:]
    barcodes = header[1:]
    return barcodes, {row[0]: {bc: int(v) for bc, v in zip(barcodes, row[1:])} for row in rows}


def discordant_qnames(bam_path):
    with pysam.AlignmentFile(bam_path, "rb") as f:
        return [r.query_name for r in f]


def test_clean_reads_counted_correctly(fixtures, tmp_path):
    result, out_prefix = run_dfn_bin(fixtures, tmp_path)
    assert result.returncode == 0, result.stderr

    barcodes, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    assert rows["chr1:1000:2000"]["AAAA-1"] == 1  # clean2
    assert rows["chr1:1000:2000"]["BBBB-1"] == 1  # clean3
    assert rows["chr2:0:1000"]["BBBB-1"] == 1     # clean4


def test_decoys_do_not_inflate_bin_zero(fixtures, tmp_path):
    # qcfail1/lowmapq1/secondary1/dup1/missingcb1 all sit in chr1:0:1000
    # alongside clean1, each meant to be dropped by a different filter.
    result, out_prefix = run_dfn_bin(fixtures, tmp_path)
    assert result.returncode == 0, result.stderr
    barcodes, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    assert rows["chr1:0:1000"]["AAAA-1"] == 1  # clean1 only


def test_offbin_read_not_counted(fixtures, tmp_path):
    result, out_prefix = run_dfn_bin(fixtures, tmp_path)
    barcodes, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    total = sum(v for row in rows.values() for v in row.values())
    # offbin1 (chr1:350000) has no containing bin, so it must never
    # appear anywhere -- cross-check via the discordant BAM being the
    # only other place a read could "disappear" to (it isn't there either).
    assert "offbin1" not in discordant_qnames(out_prefix + ".discordant.bam")
    assert total > 0  # sanity: matrix isn't trivially empty


def test_properpair_normal_counted_once_not_twice(fixtures, tmp_path):
    result, out_prefix = run_dfn_bin(fixtures, tmp_path)
    assert result.returncode == 0, result.stderr
    barcodes, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    assert rows["chr1:2000:3000"]["AAAA-1"] == 1  # not 2


def test_discordant_bam_contains_split_and_cross_chrom(fixtures, tmp_path):
    result, out_prefix = run_dfn_bin(fixtures, tmp_path)
    assert result.returncode == 0, result.stderr

    qnames = discordant_qnames(out_prefix + ".discordant.bam")
    assert qnames.count("split1") == 2
    assert qnames.count("supplonoSA1") == 1
    assert qnames.count("crosschrom1") == 2
    assert qnames.count("improper1") == 2
    assert "properpair_normal" not in qnames
    assert "clean1" not in qnames


def test_max_insert_size_default_trusts_proper_pair_flag(fixtures, tmp_path):
    result, out_prefix = run_dfn_bin(fixtures, tmp_path)
    assert result.returncode == 0, result.stderr
    qnames = discordant_qnames(out_prefix + ".discordant.bam")
    assert "properpair_bigisize" not in qnames
    # Trusted as a normal pair -> resolves to read1 (tie, equal 50bp
    # lengths) -> counted once in its bin.
    _, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    assert rows["chr1:7000:8000"]["AAAA-1"] == 1
    assert rows["chr1:8000:9000"]["AAAA-1"] == 0


def test_max_insert_size_override_diverts_big_isize_pair(fixtures, tmp_path):
    result, out_prefix = run_dfn_bin(fixtures, tmp_path, extra_args=["--max-insert-size", "10000"])
    assert result.returncode == 0, result.stderr
    qnames = discordant_qnames(out_prefix + ".discordant.bam")
    assert qnames.count("properpair_bigisize") == 2
    _, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    assert rows["chr1:7000:8000"]["AAAA-1"] == 0
    assert rows["chr1:8000:9000"]["AAAA-1"] == 0


def test_auto_mode_keeps_longer_mate(fixtures, tmp_path):
    result, out_prefix = run_dfn_bin(fixtures, tmp_path)
    assert result.returncode == 0, result.stderr
    _, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    # autolen1: 27bp read1 -> chr1:3000:4000, 51bp read2 -> chr1:4000:5000.
    # auto must keep the 51bp one.
    assert rows["chr1:3000:4000"]["AAAA-1"] == 0
    assert rows["chr1:4000:5000"]["AAAA-1"] == 1


def test_r1_mode_keeps_shorter_mate(fixtures, tmp_path):
    result, out_prefix = run_dfn_bin(fixtures, tmp_path, extra_args=["--single-end-counting", "r1"])
    assert result.returncode == 0, result.stderr
    _, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    # r1 mode always keeps read1 (the 27bp one), regardless of length.
    assert rows["chr1:3000:4000"]["AAAA-1"] == 1
    assert rows["chr1:4000:5000"]["AAAA-1"] == 0


def test_auto_mode_counts_singleton_when_mate_filtered(fixtures, tmp_path):
    result, out_prefix = run_dfn_bin(fixtures, tmp_path)
    assert result.returncode == 0, result.stderr
    _, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    assert rows["chr1:5000:6000"]["AAAA-1"] == 1  # singleton1's read1


def test_r1_mode_keeps_singleton_that_is_read1(fixtures, tmp_path):
    result, out_prefix = run_dfn_bin(fixtures, tmp_path, extra_args=["--single-end-counting", "r1"])
    assert result.returncode == 0, result.stderr
    _, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    assert rows["chr1:5000:6000"]["AAAA-1"] == 1


def test_r2_mode_drops_singleton_that_is_read1(fixtures, tmp_path):
    result, out_prefix = run_dfn_bin(fixtures, tmp_path, extra_args=["--single-end-counting", "r2"])
    assert result.returncode == 0, result.stderr
    _, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    # singleton1's only survivor is read1; r2 mode requires read2 -> dropped.
    assert rows["chr1:5000:6000"]["AAAA-1"] == 0


def test_flush_trigger_read_counted_independently(fixtures, tmp_path):
    result, out_prefix = run_dfn_bin(fixtures, tmp_path)
    assert result.returncode == 0, result.stderr
    _, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    assert rows["chr1:6000:7000"]["AAAA-1"] == 1


def test_barcode_allowlist_filters_unlisted_barcodes(fixtures, tmp_path):
    result, out_prefix = run_dfn_bin(fixtures, tmp_path, extra_args=["--barcodes", fixtures["allowlist"]])
    assert result.returncode == 0, result.stderr
    barcodes, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    assert barcodes == ["AAAA-1", "BBBB-1"]  # allowlist order
    assert "ZZZZ-1" not in barcodes


def test_notallowed_barcode_counted_when_no_allowlist_given(fixtures, tmp_path):
    result, out_prefix = run_dfn_bin(fixtures, tmp_path)
    assert result.returncode == 0, result.stderr
    barcodes, rows = read_matrix(out_prefix + ".raw_counts.txt.gz")
    assert "ZZZZ-1" in barcodes
    assert rows["chr1:9000:10000"]["ZZZZ-1"] == 1


def test_exclude_dups_toggle(fixtures, tmp_path):
    result_on, prefix_on = run_dfn_bin(fixtures, tmp_path, out_name="on")
    result_off, prefix_off = run_dfn_bin(fixtures, tmp_path, extra_args=["--exclude-dups", "false"],
                                          out_name="off")
    assert result_on.returncode == 0, result_on.stderr
    assert result_off.returncode == 0, result_off.stderr
    _, rows_on = read_matrix(prefix_on + ".raw_counts.txt.gz")
    _, rows_off = read_matrix(prefix_off + ".raw_counts.txt.gz")
    assert rows_on["chr1:0:1000"]["AAAA-1"] == 1    # dup1 excluded
    assert rows_off["chr1:0:1000"]["AAAA-1"] == 2   # dup1 now counted too


def test_missing_bam_index_errors_clearly(fixtures, tmp_path):
    os.remove(fixtures["bam"] + ".bai")
    result, _ = run_dfn_bin(fixtures, tmp_path)
    assert result.returncode != 0
    assert "index" in result.stderr.lower()


def test_malformed_bins_file_errors_clearly(fixtures, tmp_path):
    out_prefix = str(tmp_path / "out")
    result = subprocess.run(
        [DFN_BIN, "--bam", fixtures["bam"], "--bins", fixtures["bins_malformed"],
         "--out-prefix", out_prefix],
        capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "gc" in result.stderr.lower() or "column" in result.stderr.lower()


def test_chrom_naming_mismatch_errors_clearly(fixtures, tmp_path):
    out_prefix = str(tmp_path / "out")
    result = subprocess.run(
        [DFN_BIN, "--bam", fixtures["bam"], "--bins", fixtures["bins_mismatched"],
         "--out-prefix", out_prefix],
        capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "chromosome" in result.stderr.lower()


def test_empty_allowlist_errors_clearly(fixtures, tmp_path):
    result, _ = run_dfn_bin(fixtures, tmp_path, extra_args=["--barcodes", fixtures["allowlist_empty"]])
    assert result.returncode != 0
    assert "allowlist" in result.stderr.lower() or "empty" in result.stderr.lower()
```

- [ ] **Step 4: Run the test suite and fix any discrepancies**

Run: `.venv/bin/pytest tests/integration/test_dfn_bin.py -v`
Expected: all tests pass. If a count assertion is off, use `samtools view` on the fixture BAM and re-derive the exact expected count by hand from the filter/pairing rules in the design spec — don't adjust the assertion to match unexplained output. In particular, re-verify: (a) the fixture is genuinely coordinate-sorted (Step 2's `cut -f4 | uniq -c` check), since the pairing buffer's flush logic depends on it; (b) which single bin each fragment's *winning* representative lands in for each `--single-end-counting` mode, by re-tracing `resolve_pair`/`resolve_singleton` from `fragment_pairing.cpp` against that fragment's actual flag/length values in the fixture.

- [ ] **Step 5: Run everything together**

Run: `make test && .venv/bin/pytest tests/ -v`
Expected: all C++ unit tests and all Python integration tests pass.

- [ ] **Step 6: Commit**

```bash
git add tests/fixtures/gen_fixtures.py tests/integration/test_dfn_bin.py
git commit -m "Add comprehensive fixtures and integration tests for discordant detection, fragment-pairing modes, and error handling"
```

---

### Task 11: `correct_gc.py` companion script

**Files:**
- Create: `correct_gc.py`
- Create: `tests_py/test_correct_gc.py`

**Interfaces:**
- Consumes: `<prefix>.raw_counts.txt.gz` (as written by `CountsMatrix::write_gz`) and the bins TSV's `gc` column.
- Produces: `<prefix>.gc_corrected.txt.gz`, same shape/labels as the raw matrix.

- [ ] **Step 1: Install the additional Python dependencies**

Run: `.venv/bin/pip install numpy statsmodels`
Expected: installs succeed.

- [ ] **Step 2: Write the failing test**

`tests_py/test_correct_gc.py`:
```python
import gzip
import subprocess
import sys
import os
import numpy as np

REPO_ROOT = os.path.join(os.path.dirname(__file__), "..")
SCRIPT = os.path.join(REPO_ROOT, "correct_gc.py")


def write_counts(path, bins, barcodes, matrix):
    with gzip.open(path, "wt") as f:
        f.write("bin\t" + "\t".join(barcodes) + "\n")
        for bin_name, row in zip(bins, matrix):
            f.write(bin_name + "\t" + "\t".join(str(v) for v in row) + "\n")


def write_bins(path, bins, gcs):
    with open(path, "w") as f:
        f.write("chrom\tstart\tend\tgc\n")
        for bin_name, gc in zip(bins, gcs):
            chrom, start, end = bin_name.split(":")
            f.write(f"{chrom}\t{start}\t{end}\t{gc}\n")


def test_gc_trend_is_flattened(tmp_path):
    rng = np.random.default_rng(0)
    n_bins = 200
    gcs = rng.uniform(0.3, 0.6, size=n_bins)
    baseline = 100
    # Inject a strong linear GC trend plus small noise.
    counts = (baseline * (1 + 3.0 * (gcs - 0.45)) + rng.normal(0, 2, size=n_bins))
    counts = np.clip(counts, 0, None).astype(int)

    bins = [f"chr1:{i*1000}:{(i+1)*1000}" for i in range(n_bins)]
    barcodes = ["AAAA-1"]

    counts_path = str(tmp_path / "raw.txt.gz")
    bins_path = str(tmp_path / "bins.tsv")
    out_prefix = str(tmp_path / "out")

    write_counts(counts_path, bins, barcodes, counts.reshape(-1, 1))
    write_bins(bins_path, bins, gcs)

    result = subprocess.run(
        [sys.executable, SCRIPT, "--counts", counts_path, "--bins", bins_path,
         "--out-prefix", out_prefix],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr

    with gzip.open(out_prefix + ".gc_corrected.txt.gz", "rt") as f:
        lines = [line.rstrip("\n").split("\t") for line in f]
    header, rows = lines[0], lines[1:]
    assert header == ["bin", "AAAA-1"]
    corrected = np.array([int(float(row[1])) for row in rows])

    raw_corr = abs(np.corrcoef(counts, gcs)[0, 1])
    corrected_corr = abs(np.corrcoef(corrected, gcs)[0, 1])
    assert raw_corr > 0.5          # confirms the injected trend is real
    assert corrected_corr < 0.15   # confirms correction flattened it


def test_shape_and_labels_match_raw_matrix(tmp_path):
    bins = ["chr1:0:1000", "chr1:1000:2000"]
    barcodes = ["AAAA-1", "BBBB-1"]
    counts = np.array([[10, 20], [15, 5]])
    gcs = [0.4, 0.5]

    counts_path = str(tmp_path / "raw.txt.gz")
    bins_path = str(tmp_path / "bins.tsv")
    out_prefix = str(tmp_path / "out")

    write_counts(counts_path, bins, barcodes, counts)
    write_bins(bins_path, bins, gcs)

    result = subprocess.run(
        [sys.executable, SCRIPT, "--counts", counts_path, "--bins", bins_path,
         "--out-prefix", out_prefix],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr

    with gzip.open(out_prefix + ".gc_corrected.txt.gz", "rt") as f:
        lines = [line.rstrip("\n").split("\t") for line in f]
    assert lines[0] == ["bin", "AAAA-1", "BBBB-1"]
    assert [row[0] for row in lines[1:]] == bins
```

- [ ] **Step 3: Run test to verify it fails**

Run: `.venv/bin/python -m pytest tests_py/test_correct_gc.py -v`
Expected: FAIL — `correct_gc.py` does not exist yet.

- [ ] **Step 4: Write `correct_gc.py`**

```python
#!/usr/bin/env python3
"""GC-bias correction for dfn_bin's raw counts matrix.

Method: for each cell (column), fit a LOWESS curve of raw counts against
each bin's GC fraction (statsmodels.nonparametric.lowess), then rescale
the fitted values so the corrected column's median equals that cell's own
genome-wide median raw count. This is the QDNAseq/Ginkgo-style approach
(correct each bin toward the genome-wide median), applied per cell -- the
single-cell CNV norm, at the cost of one LOWESS fit per cell rather than
one fit pooled across all cells (see design spec, Component 2).
"""
import argparse
import gzip
import sys

import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess


def read_counts_matrix(path):
    with gzip.open(path, "rt") as f:
        lines = [line.rstrip("\n").split("\t") for line in f]
    header = lines[0]
    barcodes = header[1:]
    bins = []
    matrix = []
    for row in lines[1:]:
        bins.append(row[0])
        matrix.append([float(v) for v in row[1:]])
    return bins, barcodes, np.array(matrix)


def read_bin_gc(path, bins):
    gc_by_bin = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        i_chrom = header.index("chrom")
        i_start = header.index("start")
        i_end = header.index("end")
        i_gc = header.index("gc")
        for line in f:
            fields = line.rstrip("\n").split("\t")
            name = f"{fields[i_chrom]}:{fields[i_start]}:{fields[i_end]}"
            gc_by_bin[name] = float(fields[i_gc])
    return np.array([gc_by_bin[b] for b in bins])


def correct_gc(matrix, gc):
    """matrix: (n_bins, n_cells) raw counts. gc: (n_bins,) GC fraction.
    Returns a matrix of the same shape, GC-corrected per cell."""
    n_bins, n_cells = matrix.shape
    corrected = np.zeros_like(matrix)
    for cell in range(n_cells):
        counts = matrix[:, cell]
        median = np.median(counts)
        if median == 0:
            corrected[:, cell] = counts
            continue
        fitted = lowess(counts, gc, return_sorted=False)
        fitted = np.where(fitted <= 0, np.finfo(float).eps, fitted)
        corrected[:, cell] = counts * (median / fitted)
    return corrected


def main():
    parser = argparse.ArgumentParser(description="GC-correct a dfn_bin raw counts matrix")
    parser.add_argument("--counts", required=True, help="Path to <prefix>.raw_counts.txt.gz")
    parser.add_argument("--bins", required=True, help="Path to the bins.tsv used by dfn_bin")
    parser.add_argument("--out-prefix", required=True, help="Prefix for the output file")
    args = parser.parse_args()

    bins, barcodes, matrix = read_counts_matrix(args.counts)
    gc = read_bin_gc(args.bins, bins)
    corrected = correct_gc(matrix, gc)

    out_path = args.out_prefix + ".gc_corrected.txt.gz"
    with gzip.open(out_path, "wt") as f:
        f.write("bin\t" + "\t".join(barcodes) + "\n")
        for bin_name, row in zip(bins, corrected):
            f.write(bin_name + "\t" + "\t".join(f"{v:.4f}" for v in row) + "\n")


if __name__ == "__main__":
    main()
```

- [ ] **Step 5: Run test to verify it passes**

Run: `.venv/bin/python -m pytest tests_py/test_correct_gc.py -v`
Expected: `2 passed`.

- [ ] **Step 6: Commit**

```bash
git add correct_gc.py tests_py/test_correct_gc.py
git commit -m "Add correct_gc.py: per-cell LOWESS GC-bias correction"
```

---

### Task 12: README and usage documentation

**Files:**
- Create: `README.md`

**Interfaces:**
- None (documentation only).

- [ ] **Step 1: Write `README.md`**

```markdown
# dfn_bin

Turns a coordinate-sorted, indexed DEFND-seq gDNA BAM into a bin×cell
fragment-count matrix for single-cell copy-number analysis, plus a
companion BAM of structural-rearrangement evidence. Full design rationale:
`docs/superpowers/specs/2026-07-21-defnd-seq-binning-design.md`.

## Build

Requires htslib (`brew install htslib` / `apt install libhtslib-dev`) and a
C++17 compiler.

```bash
make
```

Produces `./dfn_bin`. Run `make test` to run the C++ unit test suite.

## Usage

```bash
./dfn_bin --bam sample.bam --bins bins.tsv --out-prefix sample
python3 correct_gc.py --counts sample.raw_counts.txt.gz --bins bins.tsv --out-prefix sample
```

`bins.tsv` is a tab-separated file with a header row containing (in any
column order) `chrom`, `start`, `end`, `gc`.

### Outputs

- `sample.raw_counts.txt.gz` — bins (rows) × cells (columns), raw integer
  fragment counts. Bin rows are named `<chrom>:<start>:<end>`; columns are
  cell barcodes.
- `sample.gc_corrected.txt.gz` — same shape, GC-bias-corrected (via
  `correct_gc.py`).
- `sample.discordant.bam` (+ `.bai`) — split/chimeric reads and discordant
  read pairs (cross-chromosome, or same-chromosome with the aligner's
  proper-pair flag unset), for manual/downstream structural-variant review.
  Excluded from both counts matrices.

See `./dfn_bin --help` for the full flag list and defaults.

## Validation

After a run, sanity-check the output before trusting it:

- **Assignment rate**: total counts per cell should be a large majority of
  that cell's passing-filter read count. A low fraction suggests a bins
  file that doesn't cover most of the genome, or a chromosome-naming issue
  that silently passed validation.
- **GC-correction check**: plot raw count vs. bin GC fraction — it should
  show a visible trend. Plot the same thing post-correction (using
  `sample.gc_corrected.txt.gz`) — the trend should be flat.
- **Discordant rate**: reads landing in `sample.discordant.bam` should be a
  small minority of total reads. A large fraction suggests a library or
  alignment QC problem rather than genuine genome-wide rearrangement.
  Spot-check a handful of entries in IGV.

## Testing

```bash
make test                                    # C++ unit tests
.venv/bin/pytest tests/integration -v        # end-to-end integration tests
.venv/bin/pytest tests_py -v                 # correct_gc.py tests
```

The Python virtualenv (`.venv`) and its dependencies (`pysam`, `pytest`,
`numpy`, `statsmodels`) are set up in the course of the implementation plan
(`docs/superpowers/plans/2026-07-21-defnd-seq-binning-tool.md`, Tasks 9 and 11).
```

- [ ] **Step 2: Verify the usage example actually works end to end**

Run:
```bash
.venv/bin/python -c "import sys; sys.path.insert(0,'tests/fixtures'); import gen_fixtures; print(gen_fixtures.build_all('/tmp/dfn_readme_check'))"
./dfn_bin --bam /tmp/dfn_readme_check/fixtures.bam --bins /tmp/dfn_readme_check/bins.tsv --out-prefix /tmp/dfn_readme_check/sample
.venv/bin/python correct_gc.py --counts /tmp/dfn_readme_check/sample.raw_counts.txt.gz --bins /tmp/dfn_readme_check/bins.tsv --out-prefix /tmp/dfn_readme_check/sample
ls /tmp/dfn_readme_check/sample*
```
Expected: lists `sample.raw_counts.txt.gz`, `sample.gc_corrected.txt.gz`, `sample.discordant.bam`, `sample.discordant.bam.bai`.

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "Add README with build, usage, and validation instructions"
```

---

## Self-Review

**Spec coverage:** Purpose/Background (README + Task 9 wiring) ✓. Related work sections (informational, no code needed) ✓. Scope/non-goals (enforced by design: no config files, no extra formats, no Tn5 logic) ✓. Architecture's two artifacts (`dfn_bin` Tasks 1–10, `correct_gc.py` Task 11) ✓. Assumed scale → dense matrix (Task 5) ✓. CLI contract — every flag (Task 2) ✓. `--single-end-counting` semantics (Task 7) ✓. Discordant-read detection, all 3 criteria + `--max-insert-size` (Task 6, wired in Task 9, tested in Task 10) ✓. Algorithm steps 1–6 (Task 9's `run()`) ✓. Fragment-pairing buffer (Task 7) ✓. Chromosome-naming validation (Task 3) ✓. Storage/dense growable matrix (Task 5) ✓. Threading (`hts_set_threads` in Task 9) ✓. `correct_gc.py`'s method (Task 11) ✓. Outputs, all three files (Tasks 5, 8, 11) ✓. Error handling, all four paths (Task 10's four error tests) ✓. Build/Makefile (Task 1, extended every task) ✓. Testing section's every fixture case (Task 10) ✓. Validation manual checks (README, Task 12) ✓.

**Placeholder scan:** no TBD/TODO; every step has complete, runnable code; no "similar to Task N" shortcuts — each test file is fully written out even where structurally similar to a sibling.

**Type consistency:** `Args`, `PositionMode`, `SingleEndCountingMode` (Task 2) used identically in `bins.hpp` (Task 3), `read_filter.hpp` (Task 6), `fragment_pairing.hpp` (Task 7), and `main.cpp` (Task 9). `FragmentCandidate`'s fields (`qname`, `tid`, `align_start`, `align_end`, `mate_pos`, `barcode_idx`, `is_read1`) are the same across Task 7's header, its tests, and Task 9's construction of candidates. `DiscordantReason` (Task 6) matches its three non-`None` values used in Task 9's dispatch. `BinIndex`, `BarcodeIndex`, `CountsMatrix`, `DiscordantWriter` method names match between their defining tasks and Task 9's usage.

**Arithmetic verification:** Task 10's fixture went through a full manual hand-trace against `classify_discordant`/`passes_quality_filters`/`resolve_pair`/`resolve_singleton`'s actual logic (not just intent) — the first draft had three real errors (a barcode counted under the wrong allowlist assumption, a proper-pair fixture that should have counted but was assumed diverted, and a mate-selection test whose two mates coincidentally landed in the same bin so the assertion couldn't actually distinguish which mate won). The fixture was redesigned so every fragment-pairing/discordant scenario gets its own dedicated bin, making each assertion a direct, single-cause check rather than a brittle aggregate total. Also caught: the original fixture's record order wasn't true ascending `(tid, pos)`, which the pairing buffer's flush logic depends on — fixed by writing records in verified sorted order per chromosome block.

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-07-21-defnd-seq-binning-tool.md`. Two execution options:

1. **Subagent-Driven (recommended)** — I dispatch a fresh subagent per task, review between tasks, fast iteration.
2. **Inline Execution** — Execute tasks in this session using executing-plans, batch execution with checkpoints.

Which approach?
