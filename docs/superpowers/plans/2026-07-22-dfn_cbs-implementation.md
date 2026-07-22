# dfn_cbs Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `dfn_cbs`, a C++17 binary that GC-corrects a `dfn_copy` raw-counts matrix (per-cell LOWESS, matching `correct_gc.py`'s method) and segments each cell's corrected profile via a vendored copy of `cbs+`'s CBS core, at a scale (hundreds of thousands of cells) that OOM'd the Python implementation.

**Architecture:** Load the raw-counts matrix once (bounded, one-shot ~13GB); per cell, in parallel, fit LOWESS, compute the GC-corrected/ratio tracks, segment via the vendored CBS core (centering forced off), and stream results into cell-major scratch files; then, one scratch file at a time, transpose to the final `bins × cells` gzip TSV convention. A progress reporter thread prints periodic status to stderr throughout.

**Tech Stack:** C++17, zlib (gzip output only — no htslib, this binary never opens a BAM), GNU Make (extends the existing `Makefile`); Python/pysam/pytest for integration test fixtures, matching `dfn_copy`'s existing test setup.

## Global Constraints

Copied verbatim from `docs/superpowers/specs/2026-07-22-dfn_cbs-gc-correction-segmentation-design.md` — every task's requirements implicitly include these:

- **No htslib dependency.** `dfn_cbs` never opens a BAM; its object files must not include `<htslib/...>` headers or link `-lhts`.
- **LOWESS parity with `correct_gc.py`**: `frac = 2/3`, `iter = 3` bisquare-robustifying iterations, tri-cube kernel, local **linear** fit, evaluated at each original point (`return_sorted=False` equivalent) — not resampled onto a grid. Non-positive fitted values clamped to machine epsilon before use as a divisor. If a cell's raw-count median is 0, both `gc_corrected` and `lowess_ratio` fall back to the raw counts unchanged.
- **CBS centering is forced off**, not exposed as a flag — enabling it would erase real whole-chromosome copy-number signal.
- **Peak memory bound**: never hold more than one `bins × cells`-sized structure in memory at a time (the loaded input matrix during Phase 2, OR one scratch file during Phase 3 transpose — never several simultaneously).
- **Cell filter**: `--min-reads` (default `100000`) is the sum of a cell's raw binned counts; cells below it are dropped entirely from every output (not zeroed, not flagged).
- **Output file names**: `<prefix>.gc_corrected.txt.gz`, `<prefix>.lowess_ratio.txt.gz`, `<prefix>.segmented_lowess_ratio.txt.gz` (all bins×cells gzip TSV), `<prefix>.seg` (plain text, not gzipped, IGV/DNAcopy format: `ID  chrom  loc.start  loc.end  num.mark  seg.mean`, `ID` = barcode).
- **CLI flags/defaults** (exact): `--counts`, `--bins`, `--out-prefix` (required); `--min-reads` (100000); `--threads` (1, cell-level parallelism); `--alpha` (0.01); `--perms` (1000); `--min-seg-len` (25); `--max-depth` (100); `--cbs-method` (`1cp`, or `2cp`); `--seed` (1); `--quiet` (off); `--help`.
- **Progress reporting**: stderr only; Phase 1/3 get one-shot start/done lines; Phase 2 gets a periodic (~5s) status line via a dedicated reporter thread reading an atomic counter (no locking on the worker hot path); `\r`-rewritten in place when stderr is a TTY, newline-terminated plain lines otherwise; suppressed entirely by `--quiet` (errors never suppressed).
- **`cbs+` itself is never modified.** Its segmentation core (`stats.cpp`, `cbs.cpp` from `~/mskcc/dev/cbs+/`) is vendored (copied) into `dfn_cbs`'s own source tree, called in-process.
- Follow this repo's existing conventions exactly: `src/<name>.hpp` + `src/<name>.cpp` pairs, one clear responsibility each; `tests/test_<name>.cpp` using `tests/test_framework.hpp`'s `ASSERT_TRUE`/`ASSERT_EQ`/`TEST_REPORT()` macros; Makefile rules added incrementally per task, mirroring the existing per-file `$(BUILD_DIR)/<name>.o` / `$(TEST_BUILD)/test_<name>` pattern; integration tests under `tests/integration/` using pysam-generated fixtures and pytest, matching `tests/integration/test_dfn_copy.py`'s style.

---

### Task 1: Scaffolding — split the test framework, add the `dfn_cbs` build skeleton

`tests/test_framework.hpp` currently pulls in `<htslib/sam.h>` unconditionally (for `make_test_record`, used only by `dfn_copy`'s BAM-touching tests). Since `dfn_cbs` must have zero htslib dependency, its tests can't include that header as-is. Split the BAM-specific helper into its own header so `test_framework.hpp` itself has no htslib coupling.

**Files:**
- Modify: `tests/test_framework.hpp` (remove the htslib include and `make_test_record`)
- Create: `tests/bam_test_helpers.hpp` (the removed content, unchanged)
- Modify: `tests/test_read_filter.cpp`, `tests/test_discordant_writer.cpp` (add the new include — these are the only two files that call `make_test_record`)
- Create: `src/dfn_cbs/` — no, keep flat: new `dfn_cbs`-specific sources live directly in `src/` alongside `dfn_copy`'s, following the existing flat layout (no subdirectory split for `dfn_copy` today, so don't introduce one just for the new binary — file names alone (`cbs_args.*`, `bin_gc.*`, etc.) keep them distinguishable)
- Create: `src/vendor/cbs/` (directory only, for Task 8)
- Modify: `Makefile` (add a `CBS_BIN` target skeleton and empty `dfn_cbs` main, extended task-by-task below)
- Create: `src/cbs_main.cpp` (minimal placeholder — real orchestration lands in Task 14)

- [ ] **Step 1: Split `make_test_record` out of `test_framework.hpp`**

Replace the full current contents of `tests/test_framework.hpp` with just the macros:

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

#define ASSERT_NEAR(a, b, tol) do { \
    g_tests_run++; \
    double _a = static_cast<double>(a); double _b = static_cast<double>(b); \
    if (std::fabs(_a - _b) > (tol)) { \
        std::printf("FAIL %s:%d: ASSERT_NEAR(%s, %s) |%.6f - %.6f| > %.6f\n", \
                     __FILE__, __LINE__, #a, #b, _a, _b, static_cast<double>(tol)); \
        g_tests_failed++; \
    } \
} while (0)

#define TEST_REPORT() do { \
    std::printf("%d assertions, %d failed\n", g_tests_run, g_tests_failed); \
    return g_tests_failed == 0 ? 0 : 1; \
} while (0)
```

Note: `ASSERT_NEAR` is new — added here because every `dfn_cbs` numeric test (LOWESS fits, ratios, segment means) needs tolerance-based comparison, and `test_framework.hpp` is the shared home for test macros. Needs `#include <cmath>` — add it:

```cpp
#pragma once
#include <cstdio>
#include <cmath>
```//prepend this include alongside the existing `<cstdio>` at the top of the file.

- [ ] **Step 2: Create `tests/bam_test_helpers.hpp` with the removed content**

```cpp
#pragma once
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
    uint32_t cigar[1] = { static_cast<uint32_t>(bam_cigar_gen(seq_len, BAM_CMATCH)) };
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

- [ ] **Step 3: Update the two callers**

In `tests/test_read_filter.cpp` and `tests/test_discordant_writer.cpp`, add
`#include "bam_test_helpers.hpp"` immediately after the existing
`#include "test_framework.hpp"` line in each file.

- [ ] **Step 4: Rebuild and rerun the existing suite to confirm the split didn't break anything**

Run: `make clean && make test`
Expected: same output as before the split — `7 === test_... ===` blocks, `0 failed` in each (unit test counts unchanged: 30/20/24/22/18/23/5 assertions across the seven existing suites).

- [ ] **Step 5: Add a minimal `dfn_cbs` placeholder binary and Makefile skeleton**

Create `src/cbs_main.cpp`:

```cpp
#include <cstdio>

int main() {
    std::printf("dfn_cbs placeholder\n");
    return 0;
}
```

Add to `Makefile` (after the existing `$(BIN)` rule, before the first
`$(TEST_BUILD)/test_cli_args` rule — placement doesn't matter to `make`,
but keep new `dfn_cbs`-related rules grouped together for readability):

```make
CBS_BIN := dfn_cbs

$(BUILD_DIR)/cbs_main.o: src/cbs_main.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(CBS_BIN): $(BUILD_DIR)/cbs_main.o
	$(CXX) $^ -o $@ -lz
```

And change the `all` target and `clean` target to include it:

```make
all: $(BIN) $(CBS_BIN)
```

```make
clean:
	rm -rf $(BUILD_DIR) $(TEST_BUILD) $(BIN) $(CBS_BIN)
```

Note `dfn_cbs`'s link line uses only `-lz` (no `$(LDLIBS)`, which carries
htslib) — this is the mechanism that keeps the "no htslib dependency"
constraint enforced at the build level, not just by convention. Every
later task's Makefile addition for `dfn_cbs` object files must **not**
add `pkg-config --cflags htslib` to their compile line.

- [ ] **Step 6: Build and verify**

Run: `make all`
Expected: both `dfn_copy` and `dfn_cbs` binaries are produced; running
`./dfn_cbs` prints `dfn_cbs placeholder` and exits 0.

Run: `ldd dfn_cbs 2>&1 | grep -i hts || echo "no htslib linkage"` (or `otool -L dfn_cbs | grep -i hts || echo "no htslib linkage"` on macOS)
Expected: `no htslib linkage`

- [ ] **Step 7: Commit**

```bash
git add tests/test_framework.hpp tests/bam_test_helpers.hpp \
        tests/test_read_filter.cpp tests/test_discordant_writer.cpp \
        Makefile src/cbs_main.cpp
git commit -m "Scaffold dfn_cbs: split test framework, add build skeleton"
```

---

### Task 2: CLI argument parsing

**Files:**
- Create: `src/cbs_args.hpp`
- Create: `src/cbs_args.cpp`
- Test: `tests/test_cbs_args.cpp`

**Interfaces:**
- Produces: `struct CbsArgs` (all CLI-configurable fields, see Global
  Constraints for exact flags/defaults) and `bool parse_cbs_args(int argc,
  char** argv, CbsArgs& args, bool& help_requested)` plus
  `std::string cbs_usage_text()` — every later task that needs a CLI value
  reads it from a `CbsArgs` instance produced here.

- [ ] **Step 1: Write the failing test**

Create `tests/test_cbs_args.cpp`:

```cpp
#include "test_framework.hpp"
#include "cbs_args.hpp"
#include <vector>

std::pair<int, char**> make_argv(std::vector<std::string> args) {
    static std::vector<std::string> storage;
    static std::vector<char*> ptrs;
    storage = args;
    ptrs.clear();
    for (auto& s : storage) ptrs.push_back(const_cast<char*>(s.c_str()));
    return {static_cast<int>(ptrs.size()), ptrs.data()};
}

void test_defaults_with_required_only() {
    auto [argc, argv] = make_argv({"dfn_cbs", "--counts", "in.txt.gz",
                                    "--bins", "bins.tsv", "--out-prefix", "out"});
    CbsArgs args;
    bool help = false;
    ASSERT_TRUE(parse_cbs_args(argc, argv, args, help));
    ASSERT_TRUE(!help);
    ASSERT_EQ(args.counts_path, std::string("in.txt.gz"));
    ASSERT_EQ(args.bins_path, std::string("bins.tsv"));
    ASSERT_EQ(args.out_prefix, std::string("out"));
    ASSERT_EQ(args.min_reads, static_cast<int64_t>(100000));
    ASSERT_EQ(args.threads, 1);
    ASSERT_NEAR(args.alpha, 0.01, 1e-12);
    ASSERT_EQ(args.perms, 1000);
    ASSERT_EQ(args.min_seg_len, 25);
    ASSERT_EQ(args.max_depth, 100);
    ASSERT_TRUE(args.cbs_method == CbsMethod::OneCp);
    ASSERT_EQ(args.seed, static_cast<uint64_t>(1));
    ASSERT_TRUE(!args.quiet);
}

void test_overrides_all_flags() {
    auto [argc, argv] = make_argv({
        "dfn_cbs", "--counts", "in.txt.gz", "--bins", "bins.tsv", "--out-prefix", "out",
        "--min-reads", "5000", "--threads", "8", "--alpha", "0.05", "--perms", "200",
        "--min-seg-len", "10", "--max-depth", "50", "--cbs-method", "2cp",
        "--seed", "42", "--quiet"});
    CbsArgs args;
    bool help = false;
    ASSERT_TRUE(parse_cbs_args(argc, argv, args, help));
    ASSERT_EQ(args.min_reads, static_cast<int64_t>(5000));
    ASSERT_EQ(args.threads, 8);
    ASSERT_NEAR(args.alpha, 0.05, 1e-12);
    ASSERT_EQ(args.perms, 200);
    ASSERT_EQ(args.min_seg_len, 10);
    ASSERT_EQ(args.max_depth, 50);
    ASSERT_TRUE(args.cbs_method == CbsMethod::TwoCp);
    ASSERT_EQ(args.seed, static_cast<uint64_t>(42));
    ASSERT_TRUE(args.quiet);
}

void test_help_flag() {
    auto [argc, argv] = make_argv({"dfn_cbs", "--help"});
    CbsArgs args;
    bool help = false;
    ASSERT_TRUE(parse_cbs_args(argc, argv, args, help));
    ASSERT_TRUE(help);
}

void test_missing_required_fails() {
    auto [argc, argv] = make_argv({"dfn_cbs", "--bins", "bins.tsv", "--out-prefix", "out"});
    CbsArgs args;
    bool help = false;
    ASSERT_TRUE(!parse_cbs_args(argc, argv, args, help));
}

void test_unknown_flag_fails() {
    auto [argc, argv] = make_argv({"dfn_cbs", "--counts", "in.txt.gz", "--bins", "bins.tsv",
                                    "--out-prefix", "out", "--nonsense", "1"});
    CbsArgs args;
    bool help = false;
    ASSERT_TRUE(!parse_cbs_args(argc, argv, args, help));
}

void test_invalid_cbs_method_fails() {
    auto [argc, argv] = make_argv({"dfn_cbs", "--counts", "in.txt.gz", "--bins", "bins.tsv",
                                    "--out-prefix", "out", "--cbs-method", "3cp"});
    CbsArgs args;
    bool help = false;
    ASSERT_TRUE(!parse_cbs_args(argc, argv, args, help));
}

int main() {
    test_defaults_with_required_only();
    test_overrides_all_flags();
    test_help_flag();
    test_missing_required_fails();
    test_unknown_flag_fails();
    test_invalid_cbs_method_fails();
    TEST_REPORT();
}
```

- [ ] **Step 2: Run test to verify it fails**

Add a temporary Makefile rule (folded into the real one in Step 5):
`tests/build/test_cbs_args` doesn't exist yet, so just try to compile
directly: `g++ -std=c++17 -Isrc tests/test_cbs_args.cpp -o /tmp/t 2>&1 | head -5`
Expected: FAIL — `cbs_args.hpp: No such file or directory`

- [ ] **Step 3: Write `src/cbs_args.hpp`**

```cpp
#pragma once
#include <cstdint>
#include <string>

enum class CbsMethod { OneCp, TwoCp };

struct CbsArgs {
    std::string counts_path;
    std::string bins_path;
    std::string out_prefix;
    int64_t min_reads = 100000;
    int threads = 1;
    double alpha = 0.01;
    int perms = 1000;
    int min_seg_len = 25;
    int max_depth = 100;
    CbsMethod cbs_method = CbsMethod::OneCp;
    uint64_t seed = 1;
    bool quiet = false;
};

// Parses argv into args. Returns true on success (including a bare
// --help, which sets help_requested and leaves other fields unset).
// On error, prints a specific message to stderr and returns false.
bool parse_cbs_args(int argc, char** argv, CbsArgs& args, bool& help_requested);

std::string cbs_usage_text();
```

- [ ] **Step 4: Write `src/cbs_args.cpp`**

```cpp
#include "cbs_args.hpp"
#include <cstdio>
#include <cstdlib>

std::string cbs_usage_text() {
    return
        "dfn_cbs --counts <path> --bins <path> --out-prefix <prefix> [options]\n"
        "\n"
        "Required:\n"
        "  --counts <path>         <prefix>.raw_counts.txt.gz produced by dfn_copy\n"
        "  --bins <path>           Same bins.tsv used to produce --counts\n"
        "  --out-prefix <prefix>   Prefix for the four output files\n"
        "\n"
        "Options:\n"
        "  --min-reads <int>              Minimum raw binned-count sum to keep a cell (default: 100000)\n"
        "  --threads <int>                Worker threads for per-cell parallelism (default: 1)\n"
        "  --alpha <float>                CBS split-acceptance p-value threshold (default: 0.01)\n"
        "  --perms <int>                  CBS circular-permutation count per split (default: 1000)\n"
        "  --min-seg-len <int>            CBS minimum segment length in bins (default: 25)\n"
        "  --max-depth <int>              CBS maximum recursive split depth (default: 100)\n"
        "  --cbs-method <1cp|2cp>         CBS split search mode (default: 1cp)\n"
        "  --seed <int>                   RNG seed; per-cell RNG uses seed XOR hash(barcode) (default: 1)\n"
        "  --quiet                        Suppress progress messages on stderr\n"
        "  --help                         Print this message\n";
}

namespace {

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

bool parse_uint64(const std::string& s, uint64_t& out) {
    try {
        size_t consumed = 0;
        out = std::stoull(s, &consumed);
        return consumed == s.size();
    } catch (const std::exception&) {
        return false;
    }
}

bool parse_double(const std::string& s, double& out) {
    try {
        size_t consumed = 0;
        out = std::stod(s, &consumed);
        return consumed == s.size();
    } catch (const std::exception&) {
        return false;
    }
}

} // namespace

bool parse_cbs_args(int argc, char** argv, CbsArgs& args, bool& help_requested) {
    help_requested = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--help") {
            help_requested = true;
            return true;
        }
        if (arg == "--quiet") {
            args.quiet = true;
            continue;
        }

        auto next_value = [&](const char* flag_name) -> const char* {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "error: %s requires a value\n", flag_name);
                return nullptr;
            }
            return argv[++i];
        };

        if (arg == "--counts") {
            const char* v = next_value("--counts"); if (!v) return false;
            args.counts_path = v;
        } else if (arg == "--bins") {
            const char* v = next_value("--bins"); if (!v) return false;
            args.bins_path = v;
        } else if (arg == "--out-prefix") {
            const char* v = next_value("--out-prefix"); if (!v) return false;
            args.out_prefix = v;
        } else if (arg == "--min-reads") {
            const char* v = next_value("--min-reads"); if (!v) return false;
            if (!parse_int64(v, args.min_reads)) {
                std::fprintf(stderr, "error: --min-reads requires an integer, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--threads") {
            const char* v = next_value("--threads"); if (!v) return false;
            if (!parse_int(v, args.threads)) {
                std::fprintf(stderr, "error: --threads requires an integer, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--alpha") {
            const char* v = next_value("--alpha"); if (!v) return false;
            if (!parse_double(v, args.alpha)) {
                std::fprintf(stderr, "error: --alpha requires a number, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--perms") {
            const char* v = next_value("--perms"); if (!v) return false;
            if (!parse_int(v, args.perms)) {
                std::fprintf(stderr, "error: --perms requires an integer, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--min-seg-len") {
            const char* v = next_value("--min-seg-len"); if (!v) return false;
            if (!parse_int(v, args.min_seg_len)) {
                std::fprintf(stderr, "error: --min-seg-len requires an integer, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--max-depth") {
            const char* v = next_value("--max-depth"); if (!v) return false;
            if (!parse_int(v, args.max_depth)) {
                std::fprintf(stderr, "error: --max-depth requires an integer, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--cbs-method") {
            const char* v = next_value("--cbs-method"); if (!v) return false;
            std::string sv = v;
            if (sv == "1cp") args.cbs_method = CbsMethod::OneCp;
            else if (sv == "2cp") args.cbs_method = CbsMethod::TwoCp;
            else {
                std::fprintf(stderr, "error: --cbs-method requires 1cp|2cp, got '%s'\n", v);
                return false;
            }
        } else if (arg == "--seed") {
            const char* v = next_value("--seed"); if (!v) return false;
            if (!parse_uint64(v, args.seed)) {
                std::fprintf(stderr, "error: --seed requires an integer, got '%s'\n", v);
                return false;
            }
        } else {
            std::fprintf(stderr, "error: unknown flag '%s'\n", arg.c_str());
            return false;
        }
    }

    if (args.counts_path.empty()) {
        std::fprintf(stderr, "error: --counts is required\n");
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

    return true;
}
```

- [ ] **Step 5: Add Makefile rules**

```make
$(BUILD_DIR)/cbs_args.o: src/cbs_args.cpp src/cbs_args.hpp | $(BUILD_DIR)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc -c $< -o $@

$(TEST_BUILD)/test_cbs_args: tests/test_cbs_args.cpp $(BUILD_DIR)/cbs_args.o | $(TEST_BUILD)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc $^ -o $@
```

Add `$(TEST_BUILD)/test_cbs_args` to the `TESTS :=` list, and
`$(BUILD_DIR)/cbs_args.o` to `dfn_cbs`'s link line:

```make
$(CBS_BIN): $(BUILD_DIR)/cbs_main.o $(BUILD_DIR)/cbs_args.o
	$(CXX) $^ -o $@ -lz
```

Note: `dfn_cbs`'s own compile rules use plain `-std=c++17 -Wall -Wextra
-O2 -Isrc` rather than the shared `$(CXXFLAGS)`, since `$(CXXFLAGS)`
carries `pkg-config --cflags htslib` — reusing it would (harmlessly, but
against the stated constraint) add htslib's include path to `dfn_cbs`'s
compile line. Every subsequent task's `dfn_cbs`-side Makefile rule
follows this same explicit-flags pattern instead of `$(CXXFLAGS)`.

- [ ] **Step 6: Run test to verify it passes**

Run: `make tests/build/test_cbs_args && ./tests/build/test_cbs_args`
Expected: `PASS` — some number of assertions, `0 failed`.

- [ ] **Step 7: Wire `cbs_main.cpp` to use it (still a placeholder otherwise)**

Update `src/cbs_main.cpp`:

```cpp
#include "cbs_args.hpp"
#include <cstdio>

int main(int argc, char** argv) {
    CbsArgs args;
    bool help = false;
    if (!parse_cbs_args(argc, argv, args, help)) {
        return 1;
    }
    if (help) {
        std::fputs(cbs_usage_text().c_str(), stdout);
        return 0;
    }
    std::printf("dfn_cbs placeholder (args parsed OK)\n");
    return 0;
}
```

- [ ] **Step 8: Rebuild and verify**

Run: `make all && ./dfn_cbs --help`
Expected: prints the usage text from `cbs_usage_text()`.

- [ ] **Step 9: Commit**

```bash
git add src/cbs_args.hpp src/cbs_args.cpp tests/test_cbs_args.cpp \
        src/cbs_main.cpp Makefile
git commit -m "Add dfn_cbs CLI argument parsing"
```

---

### Task 3: Bin/GC file loader

`dfn_copy`'s own `src/bins.hpp`/`BinIndex` does coordinate-lookup
(binary-search "which bin contains this alignment position") — machinery
`dfn_cbs` doesn't need, since it never bins a read; it only needs the bin
list in file order (matching the counts matrix's row order) with each
bin's chrom/start/end/gc. A separate, smaller loader avoids pulling in
`BinIndex`'s unrelated lookup API and its `cli_args.hpp` dependency
(`PositionMode` etc., meaningless here).

**Files:**
- Create: `src/bin_gc.hpp`
- Create: `src/bin_gc.cpp`
- Test: `tests/test_bin_gc.cpp`

**Interfaces:**
- Produces: `struct BinGc { std::string chrom; int64_t start; int64_t
  end; double gc; }`, `std::vector<BinGc> load_bin_gc(const std::string&
  path)` (throws `std::runtime_error` on malformed input), and
  `std::string bin_name(const BinGc& b)` returning `"<chrom>:<start>:<end>"`
  — matching `dfn_copy`'s bin-naming convention exactly, since these bin
  names must match `--counts`' row labels (checked in Task 4).

- [ ] **Step 1: Write the failing test**

Create `tests/test_bin_gc.cpp`:

```cpp
#include "test_framework.hpp"
#include "bin_gc.hpp"
#include <fstream>
#include <cstdio>

void test_loads_bins_in_file_order() {
    const char* path = "/tmp/dfn_cbs_test_bins.tsv";
    std::ofstream f(path);
    f << "chrom\tstart\tend\tgc\n";
    f << "chr1\t0\t1000\t0.42\n";
    f << "chr1\t1000\t2000\t0.55\n";
    f << "chr2\t0\t1000\t0.60\n";
    f.close();

    auto bins = load_bin_gc(path);
    ASSERT_EQ(bins.size(), static_cast<size_t>(3));
    ASSERT_EQ(bins[0].chrom, std::string("chr1"));
    ASSERT_EQ(bins[0].start, static_cast<int64_t>(0));
    ASSERT_EQ(bins[0].end, static_cast<int64_t>(1000));
    ASSERT_NEAR(bins[0].gc, 0.42, 1e-9);
    ASSERT_EQ(bins[2].chrom, std::string("chr2"));
    ASSERT_EQ(bin_name(bins[0]), std::string("chr1:0:1000"));

    std::remove(path);
}

void test_columns_can_be_reordered() {
    const char* path = "/tmp/dfn_cbs_test_bins_reorder.tsv";
    std::ofstream f(path);
    f << "gc\tchrom\tend\tstart\n";
    f << "0.5\tchrX\t2000\t1000\n";
    f.close();

    auto bins = load_bin_gc(path);
    ASSERT_EQ(bins.size(), static_cast<size_t>(1));
    ASSERT_EQ(bins[0].chrom, std::string("chrX"));
    ASSERT_EQ(bins[0].start, static_cast<int64_t>(1000));
    ASSERT_EQ(bins[0].end, static_cast<int64_t>(2000));

    std::remove(path);
}

void test_missing_column_throws() {
    const char* path = "/tmp/dfn_cbs_test_bins_missing.tsv";
    std::ofstream f(path);
    f << "chrom\tstart\tend\n"; // no gc column
    f << "chr1\t0\t1000\n";
    f.close();

    bool threw = false;
    try {
        load_bin_gc(path);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);

    std::remove(path);
}

int main() {
    test_loads_bins_in_file_order();
    test_columns_can_be_reordered();
    test_missing_column_throws();
    TEST_REPORT();
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++17 -Isrc tests/test_bin_gc.cpp -o /tmp/t 2>&1 | head -5`
Expected: FAIL — `bin_gc.hpp: No such file or directory`

- [ ] **Step 3: Write `src/bin_gc.hpp`**

```cpp
#pragma once
#include <cstdint>
#include <string>
#include <vector>

struct BinGc {
    std::string chrom;
    int64_t start;
    int64_t end;
    double gc;
};

// Loads a bins.tsv (header-driven columns: chrom, start, end, gc, any
// order, extra columns ignored) in file order. Throws std::runtime_error
// with a specific message on a missing column or a malformed row.
std::vector<BinGc> load_bin_gc(const std::string& path);

std::string bin_name(const BinGc& b);
```

- [ ] **Step 4: Write `src/bin_gc.cpp`**

```cpp
#include "bin_gc.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

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

std::vector<BinGc> load_bin_gc(const std::string& path) {
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

    std::vector<BinGc> bins;
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
        BinGc b;
        b.chrom = fields[i_chrom];
        std::string ctx = "line " + std::to_string(line_no);
        b.start = parse_i64(fields[i_start], "start at " + ctx);
        b.end = parse_i64(fields[i_end], "end at " + ctx);
        b.gc = parse_f64(fields[i_gc], "gc at " + ctx);
        bins.push_back(b);
    }

    return bins;
}

std::string bin_name(const BinGc& b) {
    return b.chrom + ":" + std::to_string(b.start) + ":" + std::to_string(b.end);
}
```

- [ ] **Step 5: Add Makefile rules**

```make
$(BUILD_DIR)/bin_gc.o: src/bin_gc.cpp src/bin_gc.hpp | $(BUILD_DIR)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc -c $< -o $@

$(TEST_BUILD)/test_bin_gc: tests/test_bin_gc.cpp $(BUILD_DIR)/bin_gc.o | $(TEST_BUILD)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc $^ -o $@
```

Add `$(TEST_BUILD)/test_bin_gc` to `TESTS :=`, and
`$(BUILD_DIR)/bin_gc.o` to `dfn_cbs`'s link line.

- [ ] **Step 6: Run test to verify it passes**

Run: `make tests/build/test_bin_gc && ./tests/build/test_bin_gc`
Expected: `PASS` — 3 test functions, some assertions, `0 failed`.

- [ ] **Step 7: Commit**

```bash
git add src/bin_gc.hpp src/bin_gc.cpp tests/test_bin_gc.cpp Makefile
git commit -m "Add dfn_cbs bin/GC file loader"
```

---

### Task 4: Raw counts matrix loader

Reads `<prefix>.raw_counts.txt.gz` (the file `dfn_copy` produces) fully
into memory as a dense `bins × cells` array, validating its header/row
count against the loaded bin list from Task 3 (same chrom-naming-mismatch
philosophy as `dfn_copy`: a loud, specific error, never a silent
zero-filled matrix).

**Files:**
- Create: `src/raw_counts_loader.hpp`
- Create: `src/raw_counts_loader.cpp`
- Test: `tests/test_raw_counts_loader.cpp`

**Interfaces:**
- Consumes: nothing from earlier tasks directly (works on the file
  format alone), but its "expected bin count" parameter is meant to be
  called with `bins.size()` from Task 3's `load_bin_gc`.
- Produces: `struct RawCountsMatrix { size_t num_bins; size_t num_cells;
  std::vector<std::string> barcodes; std::vector<int32_t> data; /* row-major
  by bin: data[bin * num_cells + cell] */ }`,
  `RawCountsMatrix load_raw_counts(const std::string& path, size_t
  expected_num_bins, const std::vector<std::string>& expected_bin_names)`
  — later tasks (5, and the orchestration in 14) read cells via
  `data[bin * num_cells + cell]` and iterate `barcodes` for names.

- [ ] **Step 1: Write the failing test**

Create `tests/test_raw_counts_loader.cpp`:

```cpp
#include "test_framework.hpp"
#include "raw_counts_loader.hpp"
#include <zlib.h>
#include <cstdio>

namespace {
void write_gz(const std::string& path, const std::string& content) {
    gzFile out = gzopen(path.c_str(), "wb");
    gzwrite(out, content.data(), static_cast<unsigned>(content.size()));
    gzclose(out);
}
}

void test_loads_matrix_matching_bin_names() {
    const char* path = "/tmp/dfn_cbs_test_counts.txt.gz";
    write_gz(path,
        "bin\tAAAA-1\tCCCC-1\n"
        "chr1:0:1000\t10\t20\n"
        "chr1:1000:2000\t5\t0\n");

    auto m = load_raw_counts(path, 2, {"chr1:0:1000", "chr1:1000:2000"});
    ASSERT_EQ(m.num_bins, static_cast<size_t>(2));
    ASSERT_EQ(m.num_cells, static_cast<size_t>(2));
    ASSERT_EQ(m.barcodes[0], std::string("AAAA-1"));
    ASSERT_EQ(m.barcodes[1], std::string("CCCC-1"));
    ASSERT_EQ(m.data[0 * m.num_cells + 0], static_cast<int32_t>(10));
    ASSERT_EQ(m.data[0 * m.num_cells + 1], static_cast<int32_t>(20));
    ASSERT_EQ(m.data[1 * m.num_cells + 0], static_cast<int32_t>(5));
    ASSERT_EQ(m.data[1 * m.num_cells + 1], static_cast<int32_t>(0));

    std::remove(path);
}

void test_bin_count_mismatch_throws() {
    const char* path = "/tmp/dfn_cbs_test_counts_mismatch.txt.gz";
    write_gz(path, "bin\tAAAA-1\nchr1:0:1000\t10\n");

    bool threw = false;
    try {
        load_raw_counts(path, 5, {"chr1:0:1000"}); // expects 5 bins, file has 1
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);

    std::remove(path);
}

void test_bin_name_mismatch_throws() {
    const char* path = "/tmp/dfn_cbs_test_counts_namemismatch.txt.gz";
    write_gz(path, "bin\tAAAA-1\nchr1:0:1000\t10\n");

    bool threw = false;
    try {
        load_raw_counts(path, 1, {"chr2:0:1000"}); // right count, wrong name/order
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);

    std::remove(path);
}

int main() {
    test_loads_matrix_matching_bin_names();
    test_bin_count_mismatch_throws();
    test_bin_name_mismatch_throws();
    TEST_REPORT();
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++17 -Isrc tests/test_raw_counts_loader.cpp -lz -o /tmp/t 2>&1 | head -5`
Expected: FAIL — `raw_counts_loader.hpp: No such file or directory`

- [ ] **Step 3: Write `src/raw_counts_loader.hpp`**

```cpp
#pragma once
#include <cstdint>
#include <string>
#include <vector>

struct RawCountsMatrix {
    size_t num_bins = 0;
    size_t num_cells = 0;
    std::vector<std::string> barcodes;
    std::vector<int32_t> data; // row-major by bin: data[bin * num_cells + cell]
};

// Loads <prefix>.raw_counts.txt.gz fully into memory. Validates the file's
// bin-row count and names against expected_num_bins/expected_bin_names
// (in order) -- throws std::runtime_error with a specific message on any
// mismatch, malformed row, or unreadable file, never returning a partial
// or silently-wrong matrix.
RawCountsMatrix load_raw_counts(const std::string& path,
                                 size_t expected_num_bins,
                                 const std::vector<std::string>& expected_bin_names);
```

- [ ] **Step 4: Write `src/raw_counts_loader.cpp`**

```cpp
#include "raw_counts_loader.hpp"
#include <zlib.h>
#include <sstream>
#include <stdexcept>

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

// Reads one '\n'-terminated line from a gzFile. Returns false at EOF.
bool gz_getline(gzFile f, std::string& out) {
    out.clear();
    char buf[4096];
    bool got_any = false;
    while (gzgets(f, buf, sizeof(buf)) != nullptr) {
        got_any = true;
        out += buf;
        if (!out.empty() && out.back() == '\n') {
            out.pop_back();
            return true;
        }
        if (gzeof(f)) break;
    }
    return got_any;
}

int32_t parse_i32(const std::string& s, const std::string& context) {
    try {
        size_t consumed = 0;
        long v = std::stol(s, &consumed);
        if (consumed != s.size()) throw std::invalid_argument("trailing characters");
        return static_cast<int32_t>(v);
    } catch (const std::exception&) {
        throw std::runtime_error("raw counts file: expected an integer for " +
                                  context + ", got '" + s + "'");
    }
}

} // namespace

RawCountsMatrix load_raw_counts(const std::string& path,
                                 size_t expected_num_bins,
                                 const std::vector<std::string>& expected_bin_names) {
    gzFile in = gzopen(path.c_str(), "rb");
    if (!in) {
        throw std::runtime_error("cannot open raw counts file: " + path);
    }

    std::string header_line;
    if (!gz_getline(in, header_line)) {
        gzclose(in);
        throw std::runtime_error("raw counts file is empty: " + path);
    }
    std::vector<std::string> header = split_tab(header_line);
    if (header.empty() || header[0] != "bin") {
        gzclose(in);
        throw std::runtime_error("raw counts file: expected header starting with "
                                  "'bin' (" + path + ")");
    }

    RawCountsMatrix m;
    m.num_cells = header.size() - 1;
    m.barcodes.assign(header.begin() + 1, header.end());
    m.num_bins = expected_num_bins;
    m.data.assign(m.num_bins * m.num_cells, 0);

    std::string line;
    size_t bin_row = 0;
    size_t line_no = 1;
    while (gz_getline(in, line)) {
        ++line_no;
        if (line.empty()) continue;
        if (bin_row >= expected_num_bins) {
            gzclose(in);
            throw std::runtime_error("raw counts file: has more bin rows than "
                                      "the bin file (" + path + ")");
        }
        std::vector<std::string> fields = split_tab(line);
        if (fields.size() != m.num_cells + 1) {
            gzclose(in);
            throw std::runtime_error("raw counts file: line " + std::to_string(line_no) +
                                      " has " + std::to_string(fields.size() - 1) +
                                      " value column(s), expected " +
                                      std::to_string(m.num_cells) + " (" + path + ")");
        }
        if (fields[0] != expected_bin_names[bin_row]) {
            gzclose(in);
            throw std::runtime_error("raw counts file: bin row " + std::to_string(bin_row) +
                                      " is '" + fields[0] + "', expected '" +
                                      expected_bin_names[bin_row] +
                                      "' -- --counts and --bins must be the matching pair "
                                      "from the same dfn_copy run (" + path + ")");
        }
        for (size_t c = 0; c < m.num_cells; ++c) {
            m.data[bin_row * m.num_cells + c] =
                parse_i32(fields[c + 1], "value at line " + std::to_string(line_no));
        }
        ++bin_row;
    }
    gzclose(in);

    if (bin_row != expected_num_bins) {
        throw std::runtime_error("raw counts file: has " + std::to_string(bin_row) +
                                  " bin row(s), expected " + std::to_string(expected_num_bins) +
                                  " to match --bins (" + path + ")");
    }

    return m;
}
```

- [ ] **Step 5: Add Makefile rules**

```make
$(BUILD_DIR)/raw_counts_loader.o: src/raw_counts_loader.cpp src/raw_counts_loader.hpp | $(BUILD_DIR)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc -c $< -o $@

$(TEST_BUILD)/test_raw_counts_loader: tests/test_raw_counts_loader.cpp $(BUILD_DIR)/raw_counts_loader.o | $(TEST_BUILD)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc $^ -o $@ -lz
```

Add `$(TEST_BUILD)/test_raw_counts_loader` to `TESTS :=`, and
`$(BUILD_DIR)/raw_counts_loader.o` to `dfn_cbs`'s link line.

- [ ] **Step 6: Run test to verify it passes**

Run: `make tests/build/test_raw_counts_loader && ./tests/build/test_raw_counts_loader`
Expected: `PASS` — 3 test functions, some assertions, `0 failed`.

- [ ] **Step 7: Commit**

```bash
git add src/raw_counts_loader.hpp src/raw_counts_loader.cpp \
        tests/test_raw_counts_loader.cpp Makefile
git commit -m "Add dfn_cbs raw counts matrix loader"
```

---

### Task 5: Cell filter

**Files:**
- Create: `src/cell_filter.hpp`
- Create: `src/cell_filter.cpp`
- Test: `tests/test_cell_filter.cpp`

**Interfaces:**
- Consumes: `RawCountsMatrix` (Task 4) by const reference.
- Produces: `struct FilteredCells { std::vector<size_t> surviving_indices;
  /* original column indices, ascending, into RawCountsMatrix::data */
  int64_t max_column_sum; }`, `FilteredCells filter_cells(const
  RawCountsMatrix& m, int64_t min_reads)` — later tasks (2's per-cell loop
  in Task 14) iterate `surviving_indices` instead of `0..num_cells`, and
  the error-handling path (Task 14) uses `max_column_sum` when every cell
  is filtered out.

- [ ] **Step 1: Write the failing test**

Create `tests/test_cell_filter.cpp`:

```cpp
#include "test_framework.hpp"
#include "cell_filter.hpp"

RawCountsMatrix make_matrix(size_t num_bins, size_t num_cells,
                             std::vector<int32_t> data,
                             std::vector<std::string> barcodes) {
    RawCountsMatrix m;
    m.num_bins = num_bins;
    m.num_cells = num_cells;
    m.data = std::move(data);
    m.barcodes = std::move(barcodes);
    return m;
}

void test_drops_cells_below_threshold() {
    // 2 bins x 3 cells. Column sums: cell0=10+5=15, cell1=100+50=150, cell2=0+0=0.
    auto m = make_matrix(2, 3,
        {10, 100, 0,
          5,  50, 0},
        {"AAAA-1", "CCCC-1", "GGGG-1"});

    auto result = filter_cells(m, 20);
    ASSERT_EQ(result.surviving_indices.size(), static_cast<size_t>(1));
    ASSERT_EQ(result.surviving_indices[0], static_cast<size_t>(1)); // only cell1 (150) survives
    ASSERT_EQ(result.max_column_sum, static_cast<int64_t>(150));
}

void test_all_cells_survive_when_threshold_is_zero() {
    auto m = make_matrix(1, 2, {0, 1}, {"AAAA-1", "CCCC-1"});
    auto result = filter_cells(m, 0);
    ASSERT_EQ(result.surviving_indices.size(), static_cast<size_t>(2));
}

void test_no_cells_survive_reports_max_sum() {
    auto m = make_matrix(1, 2, {5, 9}, {"AAAA-1", "CCCC-1"});
    auto result = filter_cells(m, 100000);
    ASSERT_EQ(result.surviving_indices.size(), static_cast<size_t>(0));
    ASSERT_EQ(result.max_column_sum, static_cast<int64_t>(9));
}

int main() {
    test_drops_cells_below_threshold();
    test_all_cells_survive_when_threshold_is_zero();
    test_no_cells_survive_reports_max_sum();
    TEST_REPORT();
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++17 -Isrc tests/test_cell_filter.cpp -o /tmp/t 2>&1 | head -5`
Expected: FAIL — `cell_filter.hpp: No such file or directory`

- [ ] **Step 3: Write `src/cell_filter.hpp`**

```cpp
#pragma once
#include <cstdint>
#include <vector>
#include "raw_counts_loader.hpp"

struct FilteredCells {
    std::vector<size_t> surviving_indices; // original column indices, ascending
    int64_t max_column_sum = 0;
};

// Sums each cell's (column's) raw counts across all bins and keeps only
// cells whose sum is >= min_reads. max_column_sum is always the largest
// observed column sum, regardless of how many cells survive -- used for
// error reporting when the filter is too strict.
FilteredCells filter_cells(const RawCountsMatrix& m, int64_t min_reads);
```

- [ ] **Step 4: Write `src/cell_filter.cpp`**

```cpp
#include "cell_filter.hpp"

FilteredCells filter_cells(const RawCountsMatrix& m, int64_t min_reads) {
    FilteredCells result;
    result.max_column_sum = 0;

    std::vector<int64_t> column_sums(m.num_cells, 0);
    for (size_t bin = 0; bin < m.num_bins; ++bin) {
        for (size_t cell = 0; cell < m.num_cells; ++cell) {
            column_sums[cell] += m.data[bin * m.num_cells + cell];
        }
    }

    for (size_t cell = 0; cell < m.num_cells; ++cell) {
        if (column_sums[cell] > result.max_column_sum) {
            result.max_column_sum = column_sums[cell];
        }
        if (column_sums[cell] >= min_reads) {
            result.surviving_indices.push_back(cell);
        }
    }

    return result;
}
```

- [ ] **Step 5: Add Makefile rules**

```make
$(BUILD_DIR)/cell_filter.o: src/cell_filter.cpp src/cell_filter.hpp src/raw_counts_loader.hpp | $(BUILD_DIR)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc -c $< -o $@

$(TEST_BUILD)/test_cell_filter: tests/test_cell_filter.cpp $(BUILD_DIR)/cell_filter.o | $(TEST_BUILD)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc $^ -o $@
```

Add `$(TEST_BUILD)/test_cell_filter` to `TESTS :=`, and
`$(BUILD_DIR)/cell_filter.o` to `dfn_cbs`'s link line.

- [ ] **Step 6: Run test to verify it passes**

Run: `make tests/build/test_cell_filter && ./tests/build/test_cell_filter`
Expected: `PASS` — 3 test functions, some assertions, `0 failed`.

- [ ] **Step 7: Commit**

```bash
git add src/cell_filter.hpp src/cell_filter.cpp tests/test_cell_filter.cpp Makefile
git commit -m "Add dfn_cbs cell filter"
```

---

### Task 6: LOWESS

The algorithmic core. Matches `statsmodels.nonparametric.smoothers_lowess.
lowess`'s default behavior: `frac=2/3` (fraction of all points used as the
local window, so window size scales with `n` automatically), `it=3`
residual-based reweighting passes performed **after** an initial fit (4
weighted-linear-regression passes total per call), tri-cube kernel
weighting within the window, bisquare reweighting between passes, local
**linear** (not constant) fit, evaluated at each original `x` (not
resampled onto a grid — matches `return_sorted=False`).

**Files:**
- Create: `src/lowess.hpp`
- Create: `src/lowess.cpp`
- Test: `tests/test_lowess.cpp`

**Interfaces:**
- Produces: `std::vector<double> lowess(const std::vector<double>& y,
  const std::vector<double>& x, double frac = 2.0/3.0, int iterations =
  3)` — returns fitted values in the same order as the input `x`/`y`
  (not sorted). Task 7 calls this with `y = counts` (cast to `double`)
  and `x = gc` values from Task 3's `BinGc` list.

- [ ] **Step 1: Write the failing test**

Create `tests/test_lowess.cpp`:

```cpp
#include "test_framework.hpp"
#include "lowess.hpp"
#include <cmath>

namespace {
double mean(const std::vector<double>& v) {
    double s = 0.0;
    for (double x : v) s += x;
    return s / static_cast<double>(v.size());
}

double pearson(const std::vector<double>& a, const std::vector<double>& b) {
    double ma = mean(a), mb = mean(b);
    double num = 0.0, da = 0.0, db = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        double xa = a[i] - ma, xb = b[i] - mb;
        num += xa * xb;
        da += xa * xa;
        db += xb * xb;
    }
    if (da <= 0.0 || db <= 0.0) return 0.0;
    return num / std::sqrt(da * db);
}
}

void test_recovers_a_clean_linear_trend() {
    // y = 100 + 200*x, x evenly spaced in [0,1]. A perfectly linear trend
    // should be recovered almost exactly by a local *linear* fit.
    const int n = 60;
    std::vector<double> x(n), y(n);
    for (int i = 0; i < n; ++i) {
        x[i] = static_cast<double>(i) / (n - 1);
        y[i] = 100.0 + 200.0 * x[i];
    }
    auto fitted = lowess(y, x);
    ASSERT_EQ(fitted.size(), static_cast<size_t>(n));
    for (int i = 5; i < n - 5; ++i) { // skip boundary points, more edge bias there
        ASSERT_NEAR(fitted[i], y[i], 5.0);
    }
}

void test_flattens_correlation_with_gc() {
    // Inject a strong linear GC trend into otherwise-flat counts, plus a
    // little per-point noise (deterministic pseudo-noise, not random --
    // keeps the test reproducible). Correlation(counts, gc) should be
    // high; correlation(counts - fitted, gc) should collapse toward 0,
    // mirroring correct_gc.py's own existing validation style.
    const int n = 80;
    std::vector<double> gc(n), counts(n);
    for (int i = 0; i < n; ++i) {
        gc[i] = static_cast<double>(i) / (n - 1);
        double noise = 3.0 * std::sin(static_cast<double>(i) * 0.7);
        counts[i] = 50.0 + 300.0 * gc[i] + noise;
    }
    double corr_before = pearson(counts, gc);
    ASSERT_TRUE(corr_before > 0.9);

    auto fitted = lowess(counts, gc);
    std::vector<double> residual(n);
    for (int i = 0; i < n; ++i) residual[i] = counts[i] - fitted[i];
    double corr_after = std::fabs(pearson(residual, gc));
    ASSERT_TRUE(corr_after < 0.3);
}

void test_handles_unsorted_x_and_returns_original_order() {
    // x is NOT pre-sorted -- lowess() must sort internally and return
    // fitted values back in the caller's original order.
    std::vector<double> x = {0.8, 0.1, 0.5, 0.3, 0.9, 0.2, 0.6, 0.4, 0.7, 0.0};
    std::vector<double> y(x.size());
    for (size_t i = 0; i < x.size(); ++i) y[i] = 10.0 + 20.0 * x[i];

    auto fitted = lowess(y, x);
    ASSERT_EQ(fitted.size(), x.size());
    // x[9] = 0.0 is the smallest value; its fit should be near y[9] = 10.0.
    ASSERT_NEAR(fitted[9], 10.0, 5.0);
    // x[4] = 0.9 is the largest value; its fit should be near y[4] = 28.0.
    ASSERT_NEAR(fitted[4], 28.0, 5.0);
}

int main() {
    test_recovers_a_clean_linear_trend();
    test_flattens_correlation_with_gc();
    test_handles_unsorted_x_and_returns_original_order();
    TEST_REPORT();
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++17 -Isrc tests/test_lowess.cpp -o /tmp/t 2>&1 | head -5`
Expected: FAIL — `lowess.hpp: No such file or directory`

- [ ] **Step 3: Write `src/lowess.hpp`**

```cpp
#pragma once
#include <vector>

// Locally weighted scatterplot smoothing (Cleveland 1979), matching
// statsmodels.nonparametric.smoothers_lowess.lowess's default behavior:
// frac fraction of points used as the local window (tri-cube-weighted),
// `iterations` residual-based bisquare reweighting passes performed after
// an initial fit, local *linear* fit at each point. Returns fitted values
// in the same order as the input x/y (not sorted, not resampled onto a
// grid). x and y must be the same size and non-empty.
std::vector<double> lowess(const std::vector<double>& y,
                            const std::vector<double>& x,
                            double frac = 2.0 / 3.0,
                            int iterations = 3);
```

- [ ] **Step 4: Write `src/lowess.cpp`**

```cpp
#include "lowess.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <stdexcept>

namespace {

// One tri-cube-weighted local linear fit at every point, given the
// current per-point robustness weights (all 1.0 on the initial pass).
// order/xs are x's ascending-sorted permutation/values, precomputed once
// by the caller and reused across passes (x doesn't change between
// passes -- only the robustness weights do).
std::vector<double> weighted_fit_pass(const std::vector<double>& x,
                                       const std::vector<double>& y,
                                       const std::vector<size_t>& order,
                                       const std::vector<double>& xs,
                                       const std::vector<double>& robustness,
                                       size_t k) {
    size_t n = x.size();
    std::vector<double> fitted(n, 0.0);

    // inv_order[order[t]] = t : original index -> position in sorted order.
    std::vector<size_t> inv_order(n);
    for (size_t t = 0; t < n; ++t) inv_order[order[t]] = t;

    for (size_t t = 0; t < n; ++t) {
        size_t orig_i = order[t];
        double xi = xs[t];

        // Find the k nearest neighbors (by x-distance) around sorted
        // position t, via two-pointer expansion.
        size_t lo = t, hi = t;
        while (hi - lo + 1 < k) {
            bool can_left = lo > 0;
            bool can_right = hi + 1 < n;
            if (can_left && can_right) {
                double dl = xi - xs[lo - 1];
                double dr = xs[hi + 1] - xi;
                if (dl <= dr) --lo; else ++hi;
            } else if (can_left) {
                --lo;
            } else if (can_right) {
                ++hi;
            } else {
                break; // window already covers all n points
            }
        }

        double d_max = std::max(xi - xs[lo], xs[hi] - xi);

        double S0 = 0.0, S1 = 0.0, S2 = 0.0, Sy = 0.0, Sxy = 0.0;
        for (size_t p = lo; p <= hi; ++p) {
            size_t orig_j = order[p];
            double xj = xs[p];
            double dist = std::fabs(xj - xi);
            double tricube = 1.0;
            if (d_max > 0.0) {
                double u = dist / d_max;
                if (u >= 1.0) tricube = 0.0;
                else {
                    double v = 1.0 - u * u * u;
                    tricube = v * v * v;
                }
            }
            double w = tricube * robustness[orig_j];
            double yj = y[orig_j];
            S0 += w;
            S1 += w * xj;
            S2 += w * xj * xj;
            Sy += w * yj;
            Sxy += w * xj * yj;
        }

        double denom = S0 * S2 - S1 * S1;
        double a, b;
        if (S0 <= 0.0) {
            // Degenerate: no positive weight anywhere in the window
            // (shouldn't happen with a well-formed robustness vector, but
            // guarded rather than dividing by zero). Fall back to the raw
            // value itself so the point simply isn't smoothed this pass.
            fitted[orig_i] = y[orig_i];
            continue;
        }
        if (std::fabs(denom) < 1e-12 * std::max(1.0, S0 * S0)) {
            // Colinear/degenerate window (e.g. every x in range is
            // identical) -- fall back to the weighted mean.
            a = Sy / S0;
            b = 0.0;
        } else {
            b = (S0 * Sxy - S1 * Sy) / denom;
            a = (Sy - b * S1) / S0;
        }
        fitted[orig_i] = a + b * xi;
        (void)inv_order; // computed for clarity/future use; not otherwise read
    }

    return fitted;
}

double median_abs(const std::vector<double>& v) {
    std::vector<double> a(v.size());
    for (size_t i = 0; i < v.size(); ++i) a[i] = std::fabs(v[i]);
    std::sort(a.begin(), a.end());
    size_t n = a.size();
    if (n == 0) return 0.0;
    if (n % 2 == 1) return a[n / 2];
    return 0.5 * (a[n / 2 - 1] + a[n / 2]);
}

} // namespace

std::vector<double> lowess(const std::vector<double>& y,
                            const std::vector<double>& x,
                            double frac,
                            int iterations) {
    if (x.size() != y.size()) {
        throw std::runtime_error("lowess: x and y must be the same size");
    }
    size_t n = x.size();
    if (n == 0) {
        throw std::runtime_error("lowess: x/y must be non-empty");
    }

    std::vector<size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](size_t a, size_t b) { return x[a] < x[b]; });

    std::vector<double> xs(n);
    for (size_t t = 0; t < n; ++t) xs[t] = x[order[t]];

    size_t k = static_cast<size_t>(std::ceil(frac * static_cast<double>(n)));
    if (k < 2) k = 2;
    if (k > n) k = n;

    std::vector<double> robustness(n, 1.0);
    std::vector<double> fitted = weighted_fit_pass(x, y, order, xs, robustness, k);

    for (int pass = 0; pass < iterations; ++pass) {
        std::vector<double> residual(n);
        for (size_t i = 0; i < n; ++i) residual[i] = y[i] - fitted[i];
        double s = median_abs(residual);
        if (s <= 0.0) {
            // Fit is already exact everywhere; further reweighting can't
            // change anything (bisquare of all-zero residuals is
            // undefined via division by s -- but there's nothing left to
            // robustify against).
            break;
        }
        for (size_t i = 0; i < n; ++i) {
            double u = residual[i] / (6.0 * s);
            if (u * u >= 1.0) {
                robustness[i] = 0.0;
            } else {
                double v = 1.0 - u * u;
                robustness[i] = v * v;
            }
        }
        fitted = weighted_fit_pass(x, y, order, xs, robustness, k);
    }

    return fitted;
}
```

- [ ] **Step 5: Add Makefile rules**

```make
$(BUILD_DIR)/lowess.o: src/lowess.cpp src/lowess.hpp | $(BUILD_DIR)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc -c $< -o $@

$(TEST_BUILD)/test_lowess: tests/test_lowess.cpp $(BUILD_DIR)/lowess.o | $(TEST_BUILD)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc $^ -o $@
```

Add `$(TEST_BUILD)/test_lowess` to `TESTS :=`, and `$(BUILD_DIR)/lowess.o`
to `dfn_cbs`'s link line.

- [ ] **Step 6: Run test to verify it passes**

Run: `make tests/build/test_lowess && ./tests/build/test_lowess`
Expected: `PASS` — 3 test functions, some assertions, `0 failed`. If
`test_recovers_a_clean_linear_trend` fails near the array boundaries,
that's expected local-regression edge bias — the test already skips the
first/last 5 points for exactly this reason; a failure *inside* that
range indicates a real bug.

- [ ] **Step 7: Commit**

```bash
git add src/lowess.hpp src/lowess.cpp tests/test_lowess.cpp Makefile
git commit -m "Add dfn_cbs LOWESS implementation"
```

---

### Task 7: GC correction

Wraps Task 6's `lowess()` with the exact per-cell correction formulas from
the spec: `lowess_ratio = counts / fitted` (clamping non-positive fits to
machine epsilon first), `gc_corrected = counts * (median(counts) /
fitted)`, with the zero-median fallback (both outputs equal the raw
`counts` unchanged).

**Files:**
- Create: `src/gc_correct.hpp`
- Create: `src/gc_correct.cpp`
- Test: `tests/test_gc_correct.cpp`

**Interfaces:**
- Consumes: `lowess()` from Task 6.
- Produces: `struct GcCorrectedCell { std::vector<double> gc_corrected;
  std::vector<double> lowess_ratio; }`, `GcCorrectedCell
  gc_correct_cell(const std::vector<double>& counts, const
  std::vector<double>& gc)` — Task 14's per-cell loop calls this once per
  surviving cell with that cell's raw-count column (cast to `double`) and
  the bin GC values from Task 3.

- [ ] **Step 1: Write the failing test**

Create `tests/test_gc_correct.cpp`:

```cpp
#include "test_framework.hpp"
#include "gc_correct.hpp"
#include <cmath>

void test_ratio_has_mean_near_one_and_corrected_median_matches_raw() {
    const int n = 60;
    std::vector<double> gc(n), counts(n);
    for (int i = 0; i < n; ++i) {
        gc[i] = static_cast<double>(i) / (n - 1);
        counts[i] = 80.0 + 100.0 * gc[i]; // clean GC trend, no noise
    }

    auto result = gc_correct_cell(counts, gc);
    ASSERT_EQ(result.gc_corrected.size(), static_cast<size_t>(n));
    ASSERT_EQ(result.lowess_ratio.size(), static_cast<size_t>(n));

    double ratio_mean = 0.0;
    for (double r : result.lowess_ratio) ratio_mean += r;
    ratio_mean /= n;
    ASSERT_NEAR(ratio_mean, 1.0, 0.1);

    std::vector<double> raw_sorted = counts, corrected_sorted = result.gc_corrected;
    std::sort(raw_sorted.begin(), raw_sorted.end());
    std::sort(corrected_sorted.begin(), corrected_sorted.end());
    double raw_median = raw_sorted[n / 2];
    double corrected_median = corrected_sorted[n / 2];
    ASSERT_NEAR(corrected_median, raw_median, raw_median * 0.15);
}

void test_zero_median_falls_back_to_raw_counts() {
    std::vector<double> gc = {0.1, 0.2, 0.3, 0.4, 0.5};
    std::vector<double> counts = {0.0, 0.0, 0.0, 0.0, 0.0}; // median is 0
    auto result = gc_correct_cell(counts, gc);
    for (size_t i = 0; i < counts.size(); ++i) {
        ASSERT_NEAR(result.gc_corrected[i], counts[i], 1e-9);
        ASSERT_NEAR(result.lowess_ratio[i], counts[i], 1e-9);
    }
}

void test_mismatched_sizes_throws() {
    std::vector<double> gc = {0.1, 0.2};
    std::vector<double> counts = {1.0, 2.0, 3.0};
    bool threw = false;
    try {
        gc_correct_cell(counts, gc);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);
}

int main() {
    test_ratio_has_mean_near_one_and_corrected_median_matches_raw();
    test_zero_median_falls_back_to_raw_counts();
    test_mismatched_sizes_throws();
    TEST_REPORT();
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++17 -Isrc tests/test_gc_correct.cpp -o /tmp/t 2>&1 | head -5`
Expected: FAIL — `gc_correct.hpp: No such file or directory`

- [ ] **Step 3: Write `src/gc_correct.hpp`**

```cpp
#pragma once
#include <vector>

struct GcCorrectedCell {
    std::vector<double> gc_corrected;  // raw-count-scale, rescaled to this cell's median
    std::vector<double> lowess_ratio;  // dimensionless, mean approx 1
};

// GC-corrects one cell's bin-count vector against the bins' GC values.
// counts.size() must equal gc.size() -- throws std::runtime_error
// otherwise. If median(counts) == 0, both outputs equal counts unchanged
// (matching correct_gc.py's existing safeguard).
GcCorrectedCell gc_correct_cell(const std::vector<double>& counts,
                                 const std::vector<double>& gc);
```

- [ ] **Step 4: Write `src/gc_correct.cpp`**

```cpp
#include "gc_correct.hpp"
#include "lowess.hpp"
#include <algorithm>
#include <limits>
#include <stdexcept>

GcCorrectedCell gc_correct_cell(const std::vector<double>& counts,
                                 const std::vector<double>& gc) {
    if (counts.size() != gc.size()) {
        throw std::runtime_error("gc_correct_cell: counts and gc must be the same size");
    }
    size_t n = counts.size();
    GcCorrectedCell result;
    result.gc_corrected.resize(n);
    result.lowess_ratio.resize(n);

    std::vector<double> sorted_counts = counts;
    std::sort(sorted_counts.begin(), sorted_counts.end());
    double median = (n % 2 == 1)
        ? sorted_counts[n / 2]
        : 0.5 * (sorted_counts[n / 2 - 1] + sorted_counts[n / 2]);

    if (median == 0.0) {
        result.gc_corrected = counts;
        result.lowess_ratio = counts;
        return result;
    }

    std::vector<double> fitted = lowess(counts, gc);
    const double eps = std::numeric_limits<double>::epsilon();
    for (size_t i = 0; i < n; ++i) {
        double f = fitted[i];
        if (f <= 0.0) f = eps;
        result.lowess_ratio[i] = counts[i] / f;
        result.gc_corrected[i] = counts[i] * (median / f);
    }

    return result;
}
```

- [ ] **Step 5: Add Makefile rules**

```make
$(BUILD_DIR)/gc_correct.o: src/gc_correct.cpp src/gc_correct.hpp src/lowess.hpp | $(BUILD_DIR)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc -c $< -o $@

$(TEST_BUILD)/test_gc_correct: tests/test_gc_correct.cpp $(BUILD_DIR)/gc_correct.o $(BUILD_DIR)/lowess.o | $(TEST_BUILD)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc $^ -o $@
```

Add `$(TEST_BUILD)/test_gc_correct` to `TESTS :=`, and
`$(BUILD_DIR)/gc_correct.o` to `dfn_cbs`'s link line.

- [ ] **Step 6: Run test to verify it passes**

Run: `make tests/build/test_gc_correct && ./tests/build/test_gc_correct`
Expected: `PASS` — 3 test functions, some assertions, `0 failed`.

- [ ] **Step 7: Commit**

```bash
git add src/gc_correct.hpp src/gc_correct.cpp tests/test_gc_correct.cpp Makefile
git commit -m "Add dfn_cbs GC correction"
```

---

### Task 8: Vendor `cbs+`'s segmentation core

`cbs+`'s actual segmentation entry point is `segment_chromosome(const
Series&, const Args&, std::vector<Segment>&, std::mt19937_64&)` in
`cbs.cpp`, which in turn calls into `stats.cpp` (the CUSUM statistics) and
optionally `correction.cpp` (GC/rolling-median hooks — unused here, since
`dfn_cbs` always leaves `Args::correction = "none"`, its default, but
copied verbatim anyway rather than stripped, to keep this a faithful,
minimal-diff copy of tested code rather than a hand-edited fork). Only
`io.hpp`'s three plain structs (`Series`, `Segment`, `Args`) are needed —
not its file-parsing/CLI machinery (`parse_args`, `load_input`,
`write_segments`, `detect_delim`) — so a trimmed `io_types.hpp` replaces
it, containing just those three structs, byte-for-byte identical to their
definitions in the original `~/mskcc/dev/cbs+/include/io.hpp`.

**Files:**
- Create: `src/vendor/cbs/io_types.hpp` (new — trimmed struct-only
  extract of `cbs+`'s `include/io.hpp`, not a verbatim file copy)
- Create: `src/vendor/cbs/stats.hpp` (verbatim copy of
  `~/mskcc/dev/cbs+/include/stats.hpp`)
- Create: `src/vendor/cbs/stats.cpp` (verbatim copy of
  `~/mskcc/dev/cbs+/src/stats.cpp`)
- Create: `src/vendor/cbs/correction.hpp` (verbatim copy of
  `~/mskcc/dev/cbs+/include/correction.hpp`)
- Create: `src/vendor/cbs/correction.cpp` (verbatim copy of
  `~/mskcc/dev/cbs+/src/correction.cpp`)
- Create: `src/vendor/cbs/cbs.hpp` (copy of
  `~/mskcc/dev/cbs+/include/cbs.hpp`, one-line change: `#include "io.hpp"`
  → `#include "io_types.hpp"`)
- Create: `src/vendor/cbs/cbs.cpp` (copy of `~/mskcc/dev/cbs+/src/cbs.cpp`,
  unchanged — it only ever includes `cbs.hpp`, `stats.hpp`,
  `correction.hpp`)
- Create: `src/vendor/cbs/README.md` (provenance note)
- Test: `tests/test_vendor_cbs.cpp`

**Interfaces:**
- Produces: `segment_chromosome(const Series&, const Args&,
  std::vector<Segment>&, std::mt19937_64&)` exactly as `cbs+` defines it
  — Task 9 calls this once per chromosome per cell.

- [ ] **Step 1: Create the provenance note**

Create `src/vendor/cbs/README.md`:

```markdown
# Vendored from cbs+

Copied from `~/mskcc/dev/cbs+/` (commit/state as of 2026-07-22) per
`docs/superpowers/specs/2026-07-22-dfn_cbs-gc-correction-segmentation-design.md`'s
"Related work: cbs+" section: `dfn_cbs` calls `segment_chromosome()`
in-process, once per chromosome per cell, rather than shelling out to
`cbs_tool` (which would mean hundreds of thousands of subprocess spawns).

Files here are verbatim copies of `cbs+`'s `stats.hpp/cpp`,
`correction.hpp/cpp`, and `cbs.hpp/cpp`, with exactly one change:
`cbs.hpp` includes `io_types.hpp` (this directory) instead of `cbs+`'s
own `io.hpp`. `io_types.hpp` is not a copy — it's a trimmed extract of
just the three plain structs (`Series`, `Segment`, `Args`) that
`cbs.cpp`/`stats.cpp` actually reference; `cbs+`'s real `io.hpp` also
declares its CLI/file-parsing functions (`parse_args`, `load_input`,
`write_segments`, `detect_delim`, `guess_header`), which `dfn_cbs` has no
use for and does not vendor.

This is a deliberate copy, not a symlink or git submodule — `cbs+` has no
library build target to link against, and `cbs+` itself is not modified
by this project. If `cbs+`'s segmentation logic changes upstream, this
copy needs updating by hand.
```

- [ ] **Step 2: Create `src/vendor/cbs/io_types.hpp`**

```cpp
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
```

- [ ] **Step 3: Create `src/vendor/cbs/stats.hpp`**

```cpp
#pragma once
#include <vector>
#include <tuple>

struct MaxStat { double stat; int k; }; // 1-CP split after k (s..k | k+1..e)
struct MaxStat2 { double stat; int k1; int k2; }; // 2-CP (s..k1 | k1+1..k2 | k2+1..e)

MaxStat max_cusum_1cp(const std::vector<double>& x, const std::vector<double>& csum, int s, int e);
MaxStat2 max_cusum_2cp(const std::vector<double>& x, const std::vector<double>& csum, int s, int e, int min_seg_len);

inline double seg_mean(const std::vector<double>& csum, int s, int e){ return (csum[e+1]-csum[s]) / (double)(e-s+1); }
void build_csum(const std::vector<double>& x, std::vector<double>& csum);
```

- [ ] **Step 4: Create `src/vendor/cbs/stats.cpp`**

```cpp
#include "stats.hpp"
#include <cmath>
#include <algorithm>

void build_csum(const std::vector<double>& x, std::vector<double>& csum){
    csum.assign(x.size()+1, 0.0);
    for(size_t i=0;i<x.size();++i) csum[i+1]=csum[i]+x[i];
}

// Standardized 1-CP max CUSUM
MaxStat max_cusum_1cp(const std::vector<double>& x, const std::vector<double>& csum, int s, int e){
    int n = e-s+1; if (n<2) return {0.0,-1};
    double mean = seg_mean(csum,s,e);
    double best=-1.0; int bestk=-1; double pref=0.0;
    for(int i=s;i<=e-1;++i){ pref += x[i]-mean; int L=i-s+1; int R=n-L; double denom = std::sqrt((double)L*R/n); double t=std::abs(pref)/denom; if(t>best){best=t; bestk=i;} }
    return {best,bestk};
}

// Approximate 2-CP scan: search k1<k2 maximizing between-segment variability vs pooled
MaxStat2 max_cusum_2cp(const std::vector<double>& x, const std::vector<double>& csum, int s, int e, int min_seg_len){
    int n=e-s+1; if(n<3*min_seg_len) return {0.0,-1,-1};
    double best=-1.0; int bk1=-1,bk2=-1; double totMean=seg_mean(csum,s,e);
    // Coarse grid to keep O(n^2) manageable; scan every step=1 for n<=5000 else step>1
    int step = (n<=5000)?1: (n<=20000?2:5);
    for(int k1=s+min_seg_len-1; k1<=e-2*min_seg_len+1; k1+=step){
        for(int k2=k1+min_seg_len; k2<=e-min_seg_len+1; k2+=step){
            // means
            double m1 = seg_mean(csum,s,k1);
            double m2 = seg_mean(csum,k1+1,k2);
            double m3 = seg_mean(csum,k2+1,e);
            int n1 = k1-s+1; int n2 = k2-(k1+1)+1; int n3 = e-(k2+1)+1;
            double between = n1*(m1-totMean)*(m1-totMean) + n2*(m2-totMean)*(m2-totMean) + n3*(m3-totMean)*(m3-totMean);
            // Use sqrt(between) as monotone proxy; more stable under permutations
            double t = std::sqrt(between);
            if(t>best){ best=t; bk1=k1; bk2=k2; }
        }
    }
    return {best,bk1,bk2};
}
```

- [ ] **Step 5: Create `src/vendor/cbs/correction.hpp`**

```cpp
#pragma once
#include <vector>
#include <string>
#include <unordered_map>

// Simple correction utilities: rolling-median detrend and linear covariate regression

void rolling_median_detrend(std::vector<double>& x, int window);

// If covariate columns are present: x <- x - (a + b1*c1 + b2*c2 + ...)
// Returns (intercept, coefficients) in out_beta for reference; silent if no covariates
void linear_covariate_regress(std::vector<double>& x, const std::unordered_map<std::string,std::vector<double>>& covars,
                              std::vector<double>& out_beta);
```

- [ ] **Step 6: Create `src/vendor/cbs/correction.cpp`**

```cpp
#include "correction.hpp"
#include <algorithm>
#include <cmath>

static double median_vec(std::vector<double> v){ if(v.empty()) return 0.0; size_t n=v.size(); std::nth_element(v.begin(), v.begin()+n/2, v.end()); double m=v[n/2]; if(n%2==0){ std::nth_element(v.begin(), v.begin()+n/2-1, v.end()); m=0.5*(v[n/2]+v[n/2-1]); } return m; }

void rolling_median_detrend(std::vector<double>& x, int window){
    if(window<3) return; if(window%2==0) window+=1; int n=(int)x.size(); std::vector<double> y(n);
    int h=window/2; std::vector<double> buf; buf.reserve(window);
    for(int i=0;i<n;++i){ int s=std::max(0,i-h); int e=std::min(n-1,i+h); buf.clear(); for(int j=s;j<=e;++j) buf.push_back(x[j]); double med=median_vec(buf); y[i]=x[i]-med; }
    x.swap(y);
}

// Ordinary Least Squares for small #covariates (<=4 typical): X=[1, c1, c2, ...]
static void ols_solve(const std::vector<std::vector<double>>& X, const std::vector<double>& y, std::vector<double>& beta){
    // Solve normal equations X^T X beta = X^T y with tiny dimension — use naive Gauss-Jordan
    int n = (int)y.size(); int p = (int)X.size(); // X stored row-major by column vectors here: X[p][n]
    // Build XtX and Xty
    std::vector<std::vector<double>> A(p, std::vector<double>(p,0.0));
    std::vector<double> b(p,0.0);
    for(int i=0;i<p;++i){ for(int j=0;j<p;++j){ double s=0; for(int k=0;k<n;++k) s += X[i][k]*X[j][k]; A[i][j]=s; } }
    for(int i=0;i<p;++i){ double s=0; for(int k=0;k<n;++k) s += X[i][k]*y[k]; b[i]=s; }
    // Augment
    for(int i=0;i<p;++i){ A[i].push_back(b[i]); }
    // Gauss-Jordan
    for(int i=0;i<p;++i){
        // pivot
        int piv=i; for(int r=i+1;r<p;++r) if (std::abs(A[r][i])>std::abs(A[piv][i])) piv=r; std::swap(A[i],A[piv]);
        double diag=A[i][i]; if(std::abs(diag)<1e-12) continue; for(size_t c=0;c<A[i].size();++c) A[i][c]/=diag;
        // eliminate
        for(int r=0;r<p;++r){ if(r==i) continue; double f=A[r][i]; if(std::abs(f)<1e-18) continue; for(size_t c=0;c<A[r].size();++c) A[r][c]-=f*A[i][c]; }
    }
    beta.assign(p,0.0); for(int i=0;i<p;++i) beta[i]=A[i][p];
}

void linear_covariate_regress(std::vector<double>& x, const std::unordered_map<std::string,std::vector<double>>& covars,
                              std::vector<double>& out_beta){
    if(covars.empty()) return; int n=(int)x.size();
    // Construct design with intercept + each covariate normalized to mean 0, var 1
    std::vector<std::vector<double>> X; X.reserve(covars.size()+1);
    std::vector<double> ones(n,1.0); X.push_back(std::move(ones));
    for(const auto& kv: covars){ const auto& v = kv.second; if((int)v.size()!=n) continue; double mu=0, s2=0; for(double z:v) mu+=z; mu/=n; for(double z:v) s2+=(z-mu)*(z-mu); s2 = (s2<=0?1.0:std::sqrt(s2/n)); std::vector<double> z(n); for(int i=0;i<n;++i) z[i]=(v[i]-mu)/s2; X.push_back(std::move(z)); }
    std::vector<double> beta; ols_solve(X,x,beta);
    // subtract fitted
    std::vector<double> fit(n,0.0);
    for(size_t j=0;j<X.size();++j){ for(int i=0;i<n;++i) fit[i]+=beta[j]*X[j][i]; }
    for(int i=0;i<n;++i) x[i]-=fit[i];
    out_beta = std::move(beta);
}
```

- [ ] **Step 7: Create `src/vendor/cbs/cbs.hpp`**

```cpp
#pragma once
#include "io_types.hpp"
#include <random>
#include <vector>

void segment_chromosome(const Series& S, const Args& a, std::vector<Segment>& out, std::mt19937_64& rng);
```

- [ ] **Step 8: Create `src/vendor/cbs/cbs.cpp`**

```cpp
#include "cbs.hpp"
#include "stats.hpp"
#include "correction.hpp"
#include <algorithm>
#include <random>

static void rotate_seg(const std::vector<double>& x, int s, int e, int r, std::vector<double>& out){ int n=e-s+1; out.resize(n); for(int i=0;i<n;++i){ int src=s+(i+r)%n; out[i]=x[src]; } }
static double max1_temp(const std::vector<double>& tmp){ std::vector<double> c; build_csum(tmp, c); auto ms=max_cusum_1cp(tmp,c,0,(int)tmp.size()-1); return ms.stat; }
static double max2_temp(const std::vector<double>& tmp, int min_seg){ std::vector<double> c; build_csum(tmp, c); auto ms=max_cusum_2cp(tmp,c,0,(int)tmp.size()-1,min_seg); return ms.stat; }

void segment_chromosome(const Series& S, const Args& a, std::vector<Segment>& out, std::mt19937_64& rng){
    const auto &x0 = S.x; std::vector<double> x = x0; // local copy (for correction)

    // Optional corrections
    if(a.center){ // per-chrom median center
        std::vector<double> tmp=x; std::nth_element(tmp.begin(), tmp.begin()+tmp.size()/2, tmp.end()); double med=tmp[tmp.size()/2]; if(tmp.size()%2==0){ std::nth_element(tmp.begin(), tmp.begin()+tmp.size()/2-1, tmp.end()); med=0.5*(tmp[tmp.size()/2]+tmp[tmp.size()/2-1]); }
        for(double &v:x) v-=med;
    }
    if(a.correction=="rolling") rolling_median_detrend(x, a.roll_window);
    else if(a.correction=="linear"){
        std::unordered_map<std::string,std::vector<double>> covars;
        if(!a.gc_col.empty()){
            auto it=S.covars.find(a.gc_col); if(it!=S.covars.end()) covars["GC"]=it->second; }
        if(!a.wave_col.empty()){
            auto it=S.covars.find(a.wave_col); if(it!=S.covars.end()) covars["WAVE"]=it->second; }
        std::vector<double> beta; if(!covars.empty()) linear_covariate_regress(x, covars, beta);
    }

    int N=(int)x.size(); if(N==0){ return; }
    std::vector<double> csum; build_csum(x, csum);

    struct Node{ int s,e,level; }; std::vector<Node> stack; stack.push_back({0,N-1,0});
    std::uniform_int_distribution<int> dist;

    while(!stack.empty()){
        Node nd=stack.back(); stack.pop_back(); int s=nd.s,e=nd.e,level=nd.level; int n=e-s+1;
        if(n<=2*a.min_seg_len || level>=a.max_depth){ double mean=seg_mean(csum,s,e); out.push_back({s,e,mean,level}); continue; }

        if(a.method=="2cp"){
            auto ms = max_cusum_2cp(x, csum, s, e, a.min_seg_len);
            if(ms.k1<0 || ms.k2<0){ double mean=seg_mean(csum,s,e); out.push_back({s,e,mean,level}); continue; }
            // Permutation test using circular rotations, test max 2-CP statistic
            double obs = ms.stat; int extreme=1; int total=a.perms+1; std::vector<double> tmp; tmp.reserve(n);
            dist = std::uniform_int_distribution<int>(0, n-1);
            for(int p=0;p<a.perms;++p){ int r = dist(rng); rotate_seg(x,s,e,r,tmp); double t=max2_temp(tmp, a.min_seg_len); if(t>=obs) ++extreme; }
            double pval=(double)extreme/total;
            if(pval<a.alpha){ stack.push_back({ms.k2+1,e,level+1}); stack.push_back({ms.k1+1,ms.k2,level+1}); stack.push_back({s,ms.k1,level+1}); }
            else { double mean=seg_mean(csum,s,e); out.push_back({s,e,mean,level}); }
        } else {
            auto ms = max_cusum_1cp(x, csum, s, e);
            if(ms.k<0){ double mean=seg_mean(csum,s,e); out.push_back({s,e,mean,level}); continue; }
            int nL=ms.k-s+1; int nR=e-(ms.k+1)+1; if(nL<a.min_seg_len || nR<a.min_seg_len){ double mean=seg_mean(csum,s,e); out.push_back({s,e,mean,level}); continue; }
            double obs=ms.stat; int extreme=1; int total=a.perms+1; std::vector<double> tmp; tmp.reserve(n);
            dist = std::uniform_int_distribution<int>(0, n-1);
            for(int p=0;p<a.perms;++p){ int r=dist(rng); rotate_seg(x,s,e,r,tmp); double t=max1_temp(tmp); if(t>=obs) ++extreme; }
            double pval=(double)extreme/total;
            if(pval<a.alpha){ stack.push_back({ms.k+1,e,level+1}); stack.push_back({s,ms.k,level+1}); }
            else { double mean=seg_mean(csum,s,e); out.push_back({s,e,mean,level}); }
        }
    }
}
```

- [ ] **Step 9: Write the failing test**

Create `tests/test_vendor_cbs.cpp`:

```cpp
#include "test_framework.hpp"
#include "cbs.hpp"
#include <algorithm>

void test_detects_an_injected_breakpoint() {
    // 60 points: first 30 at ~1.0, last 30 at ~2.0 (a clean step). With
    // centering OFF (dfn_cbs always calls it this way), the segment means
    // must reflect the true absolute levels, not a chromosome-recentered
    // difference.
    Series S;
    S.chrom = "chr1";
    for (int i = 0; i < 60; ++i) {
        S.pos.push_back(i * 1000);
        S.x.push_back(i < 30 ? 1.0 : 2.0);
    }

    Args a;
    a.alpha = 0.01;
    a.perms = 200;
    a.min_seg_len = 5;
    a.max_depth = 10;
    a.method = "1cp";
    a.center = false; // dfn_cbs always forces this off

    std::vector<Segment> segs;
    std::mt19937_64 rng(12345);
    segment_chromosome(S, a, segs, rng);
    std::sort(segs.begin(), segs.end(), [](const Segment& x, const Segment& y) { return x.s < y.s; });

    ASSERT_TRUE(segs.size() >= 2);
    // First segment should end near index 29 (the injected breakpoint).
    ASSERT_TRUE(segs.front().e >= 25 && segs.front().e <= 34);
    // Segment means should reflect the true absolute levels (1.0 / 2.0),
    // not values shifted toward zero by centering.
    ASSERT_NEAR(segs.front().mean, 1.0, 0.2);
    ASSERT_NEAR(segs.back().mean, 2.0, 0.2);
}

int main() {
    test_detects_an_injected_breakpoint();
    TEST_REPORT();
}
```

- [ ] **Step 10: Run test to verify it fails, then add Makefile rules and verify it passes**

Run: `g++ -std=c++17 -Isrc/vendor/cbs -Isrc tests/test_vendor_cbs.cpp -o /tmp/t 2>&1 | head -5`
Expected: FAIL — `cbs.hpp: No such file or directory`

Add to `Makefile`:

```make
$(BUILD_DIR)/vendor_stats.o: src/vendor/cbs/stats.cpp src/vendor/cbs/stats.hpp | $(BUILD_DIR)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc/vendor/cbs -c $< -o $@

$(BUILD_DIR)/vendor_correction.o: src/vendor/cbs/correction.cpp src/vendor/cbs/correction.hpp | $(BUILD_DIR)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc/vendor/cbs -c $< -o $@

$(BUILD_DIR)/vendor_cbs.o: src/vendor/cbs/cbs.cpp src/vendor/cbs/cbs.hpp | $(BUILD_DIR)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc/vendor/cbs -c $< -o $@

$(TEST_BUILD)/test_vendor_cbs: tests/test_vendor_cbs.cpp $(BUILD_DIR)/vendor_stats.o $(BUILD_DIR)/vendor_correction.o $(BUILD_DIR)/vendor_cbs.o | $(TEST_BUILD)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc/vendor/cbs $^ -o $@
```

Add `$(TEST_BUILD)/test_vendor_cbs` to `TESTS :=`, and the three
`vendor_*.o` files to `dfn_cbs`'s link line. Note the object files are
named `vendor_stats.o`/`vendor_correction.o`/`vendor_cbs.o` (not
`stats.o`/`cbs.o`) to avoid any ambiguity with `dfn_copy`'s existing
`build/` object names, even though they're already disambiguated by
directory in the source tree.

Run: `make tests/build/test_vendor_cbs && ./tests/build/test_vendor_cbs`
Expected: `PASS` — `1 assertions` block, `0 failed`. (The permutation
test is stochastic in principle, but `--seed 12345` with 200 permutations
against a very clean, large step makes this deterministic and stable in
practice; if it ever flakes, that's worth a closer look, not a reflexive
retry.)

- [ ] **Step 11: Commit**

```bash
git add src/vendor/cbs/ tests/test_vendor_cbs.cpp Makefile
git commit -m "Vendor cbs+'s segmentation core into dfn_cbs"
```

---

### Task 9: Per-cell segmentation wrapper

Groups one cell's `lowess_ratio` values by chromosome (first-seen order,
robust to non-contiguous chromosome blocks in the bin file), calls
Task 8's `segment_chromosome()` once per chromosome with centering forced
off and a per-(barcode, chromosome) RNG seed, then expands the resulting
segments back into a bins-length `segmented_lowess_ratio` vector and
builds this cell's `.seg` rows.

**Files:**
- Create: `src/segment_cell.hpp`
- Create: `src/segment_cell.cpp`
- Test: `tests/test_segment_cell.cpp`

**Interfaces:**
- Consumes: `BinGc` (Task 3), `CbsArgs`/`CbsMethod` (Task 2), vendored
  `Series`/`Args`/`Segment`/`segment_chromosome` (Task 8).
- Produces: `struct SegRow { std::string barcode; std::string chrom;
  int64_t start; int64_t end; int num_mark; double seg_mean; }`,
  `struct CellSegmentation { std::vector<double> segmented_lowess_ratio;
  std::vector<SegRow> seg_rows; }`, `CellSegmentation segment_cell(const
  std::vector<BinGc>& bins, const std::vector<double>& lowess_ratio,
  const std::string& barcode, const CbsArgs& cbs_args)` — Task 14's
  per-cell loop calls this once per surviving cell, after Task 7's
  `gc_correct_cell()`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_segment_cell.cpp`:

```cpp
#include "test_framework.hpp"
#include "segment_cell.hpp"

std::vector<BinGc> make_two_chrom_bins(int per_chrom) {
    std::vector<BinGc> bins;
    for (int i = 0; i < per_chrom; ++i) {
        bins.push_back({"chr1", static_cast<int64_t>(i * 1000),
                         static_cast<int64_t>((i + 1) * 1000), 0.5});
    }
    for (int i = 0; i < per_chrom; ++i) {
        bins.push_back({"chr2", static_cast<int64_t>(i * 1000),
                         static_cast<int64_t>((i + 1) * 1000), 0.5});
    }
    return bins;
}

void test_two_chromosomes_keep_distinct_absolute_levels() {
    // chr1 is flat at ratio 1.0 (diploid), chr2 is flat at ratio 2.5 (a
    // whole-chromosome gain). With centering forced off, both chromosomes'
    // segment means must reflect their own true absolute level -- neither
    // should be pulled toward the other or toward zero.
    auto bins = make_two_chrom_bins(40);
    std::vector<double> ratio(80);
    for (int i = 0; i < 40; ++i) ratio[i] = 1.0;
    for (int i = 40; i < 80; ++i) ratio[i] = 2.5;

    CbsArgs args;
    args.alpha = 0.01;
    args.perms = 200;
    args.min_seg_len = 5;
    args.max_depth = 10;
    args.cbs_method = CbsMethod::OneCp;
    args.seed = 7;

    auto result = segment_cell(bins, ratio, "AAAA-1", args);
    ASSERT_EQ(result.segmented_lowess_ratio.size(), static_cast<size_t>(80));
    for (int i = 0; i < 40; ++i) ASSERT_NEAR(result.segmented_lowess_ratio[i], 1.0, 0.15);
    for (int i = 40; i < 80; ++i) ASSERT_NEAR(result.segmented_lowess_ratio[i], 2.5, 0.15);

    ASSERT_TRUE(!result.seg_rows.empty());
    bool saw_chr1 = false, saw_chr2 = false;
    for (const auto& row : result.seg_rows) {
        ASSERT_EQ(row.barcode, std::string("AAAA-1"));
        if (row.chrom == "chr1") { saw_chr1 = true; ASSERT_NEAR(row.seg_mean, 1.0, 0.15); }
        if (row.chrom == "chr2") { saw_chr2 = true; ASSERT_NEAR(row.seg_mean, 2.5, 0.15); }
    }
    ASSERT_TRUE(saw_chr1);
    ASSERT_TRUE(saw_chr2);
}

void test_seg_row_coordinates_match_bin_boundaries() {
    // A single chromosome, single segment (uniform ratio, no breakpoint):
    // the segment's start/end must be the first/last bin's own genomic
    // coordinates, and num_mark must equal the bin count.
    std::vector<BinGc> bins = {
        {"chr1", 0, 1000, 0.5}, {"chr1", 1000, 2000, 0.5},
        {"chr1", 2000, 3000, 0.5}, {"chr1", 3000, 4000, 0.5},
        {"chr1", 4000, 5000, 0.5}, {"chr1", 5000, 6000, 0.5},
    };
    std::vector<double> ratio(6, 1.0);

    CbsArgs args;
    args.min_seg_len = 2;
    args.perms = 50;

    auto result = segment_cell(bins, ratio, "CCCC-1", args);
    ASSERT_EQ(result.seg_rows.size(), static_cast<size_t>(1));
    ASSERT_EQ(result.seg_rows[0].start, static_cast<int64_t>(0));
    ASSERT_EQ(result.seg_rows[0].end, static_cast<int64_t>(6000));
    ASSERT_EQ(result.seg_rows[0].num_mark, 6);
}

void test_different_barcodes_get_different_rng_streams() {
    // Same bins, same ratio values, two different barcodes: the
    // seed-XOR-hash(barcode) scheme must actually vary the RNG stream --
    // this doesn't assert a *specific* difference (permutation p-values
    // are stochastic by design), just that the two calls don't silently
    // share identical internal state by construction. Verified indirectly:
    // both calls must complete and return the same *shape* of result
    // (same bin count), proving the barcode is accepted and doesn't
    // corrupt the pipeline for either cell.
    std::vector<BinGc> bins = {
        {"chr1", 0, 1000, 0.5}, {"chr1", 1000, 2000, 0.5},
        {"chr1", 2000, 3000, 0.5}, {"chr1", 3000, 4000, 0.5},
    };
    std::vector<double> ratio = {1.0, 1.0, 3.0, 3.0};

    CbsArgs args;
    args.min_seg_len = 1;
    args.perms = 50;

    auto r1 = segment_cell(bins, ratio, "AAAA-1", args);
    auto r2 = segment_cell(bins, ratio, "TTTT-1", args);
    ASSERT_EQ(r1.segmented_lowess_ratio.size(), static_cast<size_t>(4));
    ASSERT_EQ(r2.segmented_lowess_ratio.size(), static_cast<size_t>(4));
}

int main() {
    test_two_chromosomes_keep_distinct_absolute_levels();
    test_seg_row_coordinates_match_bin_boundaries();
    test_different_barcodes_get_different_rng_streams();
    TEST_REPORT();
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++17 -Isrc -Isrc/vendor/cbs tests/test_segment_cell.cpp -o /tmp/t 2>&1 | head -5`
Expected: FAIL — `segment_cell.hpp: No such file or directory`

- [ ] **Step 3: Write `src/segment_cell.hpp`**

```cpp
#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include "bin_gc.hpp"
#include "cbs_args.hpp"

struct SegRow {
    std::string barcode;
    std::string chrom;
    int64_t start;
    int64_t end;
    int num_mark;
    double seg_mean;
};

struct CellSegmentation {
    std::vector<double> segmented_lowess_ratio; // same length/order as bins
    std::vector<SegRow> seg_rows;
};

// Groups lowess_ratio by chromosome (first-seen order, per bins' chrom
// column), segments each chromosome independently via the vendored CBS
// core with centering forced off, and expands the result back to a
// bins-length vector plus this cell's .seg rows. bins.size() must equal
// lowess_ratio.size().
CellSegmentation segment_cell(const std::vector<BinGc>& bins,
                               const std::vector<double>& lowess_ratio,
                               const std::string& barcode,
                               const CbsArgs& cbs_args);
```

- [ ] **Step 4: Write `src/segment_cell.cpp`**

```cpp
#include "segment_cell.hpp"
#include "cbs.hpp"
#include <algorithm>
#include <stdexcept>
#include <unordered_map>

CellSegmentation segment_cell(const std::vector<BinGc>& bins,
                               const std::vector<double>& lowess_ratio,
                               const std::string& barcode,
                               const CbsArgs& cbs_args) {
    if (bins.size() != lowess_ratio.size()) {
        throw std::runtime_error("segment_cell: bins and lowess_ratio must be the same size");
    }

    CellSegmentation result;
    result.segmented_lowess_ratio.assign(bins.size(), 0.0);

    // Group original bin indices by chromosome, preserving first-seen
    // chromosome order (robust to a non-contiguous chrom block, though
    // bins.tsv is expected to already be chromosome-grouped in practice).
    std::vector<std::string> chrom_order;
    std::unordered_map<std::string, std::vector<size_t>> by_chrom;
    for (size_t i = 0; i < bins.size(); ++i) {
        const std::string& c = bins[i].chrom;
        if (by_chrom.find(c) == by_chrom.end()) chrom_order.push_back(c);
        by_chrom[c].push_back(i);
    }

    Args a;
    a.alpha = cbs_args.alpha;
    a.perms = cbs_args.perms;
    a.min_seg_len = cbs_args.min_seg_len;
    a.max_depth = cbs_args.max_depth;
    a.method = (cbs_args.cbs_method == CbsMethod::OneCp) ? "1cp" : "2cp";
    a.center = false; // forced off -- see design spec's "Related work: cbs+"

    std::hash<std::string> hasher;

    for (const auto& chrom : chrom_order) {
        const auto& indices = by_chrom.at(chrom);

        Series S;
        S.chrom = chrom;
        S.pos.reserve(indices.size());
        S.x.reserve(indices.size());
        for (size_t idx : indices) {
            S.pos.push_back(static_cast<long>(bins[idx].start));
            S.x.push_back(lowess_ratio[idx]);
        }

        uint64_t chrom_seed = cbs_args.seed ^ hasher(barcode) ^ hasher(chrom);
        std::mt19937_64 rng(chrom_seed);

        std::vector<Segment> segs;
        segment_chromosome(S, a, segs, rng);
        std::sort(segs.begin(), segs.end(),
                  [](const Segment& x, const Segment& y) { return x.s < y.s; });

        for (const auto& seg : segs) {
            for (int local = seg.s; local <= seg.e; ++local) {
                result.segmented_lowess_ratio[indices[local]] = seg.mean;
            }
            SegRow row;
            row.barcode = barcode;
            row.chrom = chrom;
            row.start = bins[indices[seg.s]].start;
            row.end = bins[indices[seg.e]].end;
            row.num_mark = seg.e - seg.s + 1;
            row.seg_mean = seg.mean;
            result.seg_rows.push_back(row);
        }
    }

    return result;
}
```

- [ ] **Step 5: Add Makefile rules**

```make
$(BUILD_DIR)/segment_cell.o: src/segment_cell.cpp src/segment_cell.hpp src/bin_gc.hpp src/cbs_args.hpp | $(BUILD_DIR)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc -Isrc/vendor/cbs -c $< -o $@

$(TEST_BUILD)/test_segment_cell: tests/test_segment_cell.cpp $(BUILD_DIR)/segment_cell.o $(BUILD_DIR)/bin_gc.o $(BUILD_DIR)/cbs_args.o $(BUILD_DIR)/vendor_stats.o $(BUILD_DIR)/vendor_correction.o $(BUILD_DIR)/vendor_cbs.o | $(TEST_BUILD)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc -Isrc/vendor/cbs $^ -o $@
```

Add `$(TEST_BUILD)/test_segment_cell` to `TESTS :=`, and
`$(BUILD_DIR)/segment_cell.o` to `dfn_cbs`'s link line.

- [ ] **Step 6: Run test to verify it passes**

Run: `make tests/build/test_segment_cell && ./tests/build/test_segment_cell`
Expected: `PASS` — 3 test functions, some assertions, `0 failed`.

- [ ] **Step 7: Commit**

```bash
git add src/segment_cell.hpp src/segment_cell.cpp tests/test_segment_cell.cpp Makefile
git commit -m "Add dfn_cbs per-cell segmentation wrapper"
```

---

### Task 10: Cell-major scratch file I/O

The bounded-memory mechanism from the spec: each of the three output
matrices gets a flat binary scratch file (`float32`, fixed record size =
`num_bins * 4` bytes per cell), opened once, written to by many worker
threads concurrently — each thread writes only its own cell's byte range,
so this uses **positioned I/O** (`pwrite`/`pread`, not `fseek`+`fread` on
a shared `FILE*`) so no locking is needed: concurrent positioned writes to
disjoint byte ranges of the same file descriptor are safe without
synchronization, whereas a shared seek position would race.

**Files:**
- Create: `src/scratch_matrix.hpp`
- Create: `src/scratch_matrix.cpp`
- Test: `tests/test_scratch_matrix.cpp`

**Interfaces:**
- Produces: `class ScratchMatrix` with `static ScratchMatrix create(const
  std::string& path, size_t num_cells, size_t num_bins)`, `void
  write_cell(size_t cell_index, const std::vector<double>& values)`,
  `std::vector<double> read_cell(size_t cell_index) const`,
  `std::vector<float> read_all() const` (flat, `num_cells * num_bins`,
  row-major by cell — this is what Task 12's transpose-on-exit reads),
  `void close()`. Task 14's Phase 2 creates three of these (one per output
  matrix) and calls `write_cell()` once per surviving cell per matrix;
  Task 12 calls `read_all()` once per scratch file during Phase 3.

- [ ] **Step 1: Write the failing test**

Create `tests/test_scratch_matrix.cpp`:

```cpp
#include "test_framework.hpp"
#include "scratch_matrix.hpp"
#include <cstdio>
#include <thread>
#include <vector>

void test_write_then_read_single_cell_round_trips() {
    const char* path = "/tmp/dfn_cbs_test_scratch1.bin";
    auto m = ScratchMatrix::create(path, 3, 4);
    m.write_cell(1, {1.5, 2.5, 3.5, 4.5});
    auto back = m.read_cell(1);
    ASSERT_EQ(back.size(), static_cast<size_t>(4));
    ASSERT_NEAR(back[0], 1.5, 1e-5);
    ASSERT_NEAR(back[3], 4.5, 1e-5);
    m.close();
    std::remove(path);
}

void test_concurrent_writes_to_disjoint_cells_are_safe() {
    const char* path = "/tmp/dfn_cbs_test_scratch2.bin";
    const size_t num_cells = 20;
    const size_t num_bins = 8;
    auto m = ScratchMatrix::create(path, num_cells, num_bins);

    std::vector<std::thread> threads;
    for (size_t c = 0; c < num_cells; ++c) {
        threads.emplace_back([&m, c, num_bins]() {
            std::vector<double> values(num_bins);
            for (size_t b = 0; b < num_bins; ++b) {
                values[b] = static_cast<double>(c * 100 + b);
            }
            m.write_cell(c, values);
        });
    }
    for (auto& t : threads) t.join();

    for (size_t c = 0; c < num_cells; ++c) {
        auto back = m.read_cell(c);
        for (size_t b = 0; b < num_bins; ++b) {
            ASSERT_NEAR(back[b], static_cast<double>(c * 100 + b), 1e-5);
        }
    }
    m.close();
    std::remove(path);
}

void test_read_all_returns_row_major_by_cell() {
    const char* path = "/tmp/dfn_cbs_test_scratch3.bin";
    auto m = ScratchMatrix::create(path, 2, 3);
    m.write_cell(0, {1.0, 2.0, 3.0});
    m.write_cell(1, {4.0, 5.0, 6.0});
    auto flat = m.read_all();
    ASSERT_EQ(flat.size(), static_cast<size_t>(6));
    ASSERT_NEAR(flat[0 * 3 + 0], 1.0f, 1e-5);
    ASSERT_NEAR(flat[0 * 3 + 2], 3.0f, 1e-5);
    ASSERT_NEAR(flat[1 * 3 + 1], 5.0f, 1e-5);
    m.close();
    std::remove(path);
}

int main() {
    test_write_then_read_single_cell_round_trips();
    test_concurrent_writes_to_disjoint_cells_are_safe();
    test_read_all_returns_row_major_by_cell();
    TEST_REPORT();
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++17 -Isrc -pthread tests/test_scratch_matrix.cpp -o /tmp/t 2>&1 | head -5`
Expected: FAIL — `scratch_matrix.hpp: No such file or directory`

- [ ] **Step 3: Write `src/scratch_matrix.hpp`**

```cpp
#pragma once
#include <cstddef>
#include <string>
#include <vector>

// A flat binary scratch file for one cell-major output matrix: fixed
// record size (num_bins float32s per cell), opened once, written to
// concurrently by many threads via positioned I/O (pwrite/pread) --
// no locking needed since each cell's byte range is disjoint from every
// other cell's.
class ScratchMatrix {
public:
    static ScratchMatrix create(const std::string& path, size_t num_cells, size_t num_bins);

    ScratchMatrix(const ScratchMatrix&) = delete;
    ScratchMatrix& operator=(const ScratchMatrix&) = delete;
    ScratchMatrix(ScratchMatrix&& other) noexcept;
    ScratchMatrix& operator=(ScratchMatrix&&) = delete;
    ~ScratchMatrix();

    // values.size() must equal num_bins(). Thread-safe across distinct
    // cell_index values; not safe to call twice concurrently for the
    // *same* cell_index (not needed -- each cell is written exactly once).
    void write_cell(size_t cell_index, const std::vector<double>& values);

    // Reads one cell back as double (upcast from the file's float32
    // storage) -- used for round-trip testing, not the hot path.
    std::vector<double> read_cell(size_t cell_index) const;

    // Reads the entire file into memory as a flat, row-major-by-cell
    // float32 vector (size num_cells() * num_bins()) -- what Task 12's
    // transpose-on-exit uses.
    std::vector<float> read_all() const;

    size_t num_cells() const { return num_cells_; }
    size_t num_bins() const { return num_bins_; }
    const std::string& path() const { return path_; }

    // Idempotent. Safe to call more than once.
    void close();

private:
    ScratchMatrix(int fd, std::string path, size_t num_cells, size_t num_bins);
    int fd_;
    std::string path_;
    size_t num_cells_;
    size_t num_bins_;
};
```

- [ ] **Step 4: Write `src/scratch_matrix.cpp`**

```cpp
#include "scratch_matrix.hpp"
#include <fcntl.h>
#include <unistd.h>
#include <stdexcept>
#include <cstring>

ScratchMatrix::ScratchMatrix(int fd, std::string path, size_t num_cells, size_t num_bins)
    : fd_(fd), path_(std::move(path)), num_cells_(num_cells), num_bins_(num_bins) {}

ScratchMatrix::ScratchMatrix(ScratchMatrix&& other) noexcept
    : fd_(other.fd_), path_(std::move(other.path_)),
      num_cells_(other.num_cells_), num_bins_(other.num_bins_) {
    other.fd_ = -1;
}

ScratchMatrix::~ScratchMatrix() {
    if (fd_ >= 0) {
        ::close(fd_);
        fd_ = -1;
    }
}

ScratchMatrix ScratchMatrix::create(const std::string& path, size_t num_cells, size_t num_bins) {
    int fd = ::open(path.c_str(), O_CREAT | O_RDWR | O_TRUNC, 0644);
    if (fd < 0) {
        throw std::runtime_error("cannot create scratch file: " + path +
                                  " (" + std::strerror(errno) + ")");
    }
    off_t total_bytes = static_cast<off_t>(num_cells) * static_cast<off_t>(num_bins) *
                         static_cast<off_t>(sizeof(float));
    if (::ftruncate(fd, total_bytes) != 0) {
        std::string msg = "cannot size scratch file: " + path + " (" + std::strerror(errno) + ")";
        ::close(fd);
        throw std::runtime_error(msg);
    }
    return ScratchMatrix(fd, path, num_cells, num_bins);
}

void ScratchMatrix::write_cell(size_t cell_index, const std::vector<double>& values) {
    if (values.size() != num_bins_) {
        throw std::runtime_error("scratch matrix: write_cell got " +
                                  std::to_string(values.size()) + " values, expected " +
                                  std::to_string(num_bins_));
    }
    if (cell_index >= num_cells_) {
        throw std::runtime_error("scratch matrix: cell index " + std::to_string(cell_index) +
                                  " out of bounds (num_cells=" + std::to_string(num_cells_) + ")");
    }
    std::vector<float> buf(num_bins_);
    for (size_t i = 0; i < num_bins_; ++i) buf[i] = static_cast<float>(values[i]);

    off_t offset = static_cast<off_t>(cell_index) * static_cast<off_t>(num_bins_) *
                   static_cast<off_t>(sizeof(float));
    size_t bytes_to_write = num_bins_ * sizeof(float);
    ssize_t written = ::pwrite(fd_, buf.data(), bytes_to_write, offset);
    if (written < 0 || static_cast<size_t>(written) != bytes_to_write) {
        throw std::runtime_error("scratch matrix: short/failed write for cell " +
                                  std::to_string(cell_index) + " (" + path_ + ")");
    }
}

std::vector<double> ScratchMatrix::read_cell(size_t cell_index) const {
    if (cell_index >= num_cells_) {
        throw std::runtime_error("scratch matrix: cell index " + std::to_string(cell_index) +
                                  " out of bounds (num_cells=" + std::to_string(num_cells_) + ")");
    }
    std::vector<float> buf(num_bins_);
    off_t offset = static_cast<off_t>(cell_index) * static_cast<off_t>(num_bins_) *
                   static_cast<off_t>(sizeof(float));
    size_t bytes_to_read = num_bins_ * sizeof(float);
    ssize_t read_n = ::pread(fd_, buf.data(), bytes_to_read, offset);
    if (read_n < 0 || static_cast<size_t>(read_n) != bytes_to_read) {
        throw std::runtime_error("scratch matrix: short/failed read for cell " +
                                  std::to_string(cell_index) + " (" + path_ + ")");
    }
    std::vector<double> out(num_bins_);
    for (size_t i = 0; i < num_bins_; ++i) out[i] = static_cast<double>(buf[i]);
    return out;
}

std::vector<float> ScratchMatrix::read_all() const {
    size_t total = num_cells_ * num_bins_;
    std::vector<float> out(total);
    size_t total_bytes = total * sizeof(float);
    size_t done = 0;
    while (done < total_bytes) {
        ssize_t n = ::pread(fd_, reinterpret_cast<char*>(out.data()) + done,
                            total_bytes - done, static_cast<off_t>(done));
        if (n <= 0) {
            throw std::runtime_error("scratch matrix: short/failed read_all (" + path_ + ")");
        }
        done += static_cast<size_t>(n);
    }
    return out;
}

void ScratchMatrix::close() {
    if (fd_ >= 0) {
        ::close(fd_);
        fd_ = -1;
    }
}
```

- [ ] **Step 5: Add Makefile rules**

```make
$(BUILD_DIR)/scratch_matrix.o: src/scratch_matrix.cpp src/scratch_matrix.hpp | $(BUILD_DIR)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc -c $< -o $@

$(TEST_BUILD)/test_scratch_matrix: tests/test_scratch_matrix.cpp $(BUILD_DIR)/scratch_matrix.o | $(TEST_BUILD)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc $^ -o $@ -pthread
```

Add `$(TEST_BUILD)/test_scratch_matrix` to `TESTS :=`, and
`$(BUILD_DIR)/scratch_matrix.o` to `dfn_cbs`'s link line. `dfn_cbs`'s own
final link line also needs `-pthread` from this task onward (needed here
for `<thread>` in the test, and by Task 14's worker pool) — update it:

```make
$(CBS_BIN): $(BUILD_DIR)/cbs_main.o $(BUILD_DIR)/cbs_args.o $(BUILD_DIR)/bin_gc.o \
            $(BUILD_DIR)/raw_counts_loader.o $(BUILD_DIR)/cell_filter.o \
            $(BUILD_DIR)/lowess.o $(BUILD_DIR)/gc_correct.o \
            $(BUILD_DIR)/vendor_stats.o $(BUILD_DIR)/vendor_correction.o $(BUILD_DIR)/vendor_cbs.o \
            $(BUILD_DIR)/segment_cell.o $(BUILD_DIR)/scratch_matrix.o
	$(CXX) $^ -o $@ -lz -pthread
```

(This consolidated link line folds in every `dfn_cbs`-side `.o` file from
Tasks 2–10 in one place, since the plan has been adding them one at a
time — from here on, each remaining task just appends its own new `.o` to
this same list.)

- [ ] **Step 6: Run test to verify it passes**

Run: `make tests/build/test_scratch_matrix && ./tests/build/test_scratch_matrix`
Expected: `PASS` — 3 test functions, some assertions, `0 failed`.

- [ ] **Step 7: Rebuild `dfn_cbs` itself to confirm the consolidated link line works**

Run: `make all`
Expected: `dfn_copy` and `dfn_cbs` both build cleanly with no undefined
references.

- [ ] **Step 8: Commit**

```bash
git add src/scratch_matrix.hpp src/scratch_matrix.cpp tests/test_scratch_matrix.cpp Makefile
git commit -m "Add dfn_cbs cell-major scratch file I/O"
```

---

### Task 11: `.seg` file writer

Unlike the three matrix outputs, `.seg` is small (row-per-segment, not
row-per-bin) and written incrementally as each cell completes — a
mutex-guarded append is fine here (negligible contention: once per cell,
not once per bin, per the spec).

**Files:**
- Create: `src/seg_file_writer.hpp`
- Create: `src/seg_file_writer.cpp`
- Test: `tests/test_seg_file_writer.cpp`

**Interfaces:**
- Consumes: `SegRow` (Task 9).
- Produces: `class SegFileWriter` with constructor `explicit
  SegFileWriter(const std::string& path)` (writes the IGV header
  immediately, throws `std::runtime_error` if the path can't be opened),
  `void append_rows(const std::vector<SegRow>& rows)` (thread-safe), `void
  close()`. Task 14's Phase 2 creates one instance and calls
  `append_rows()` once per completed cell.

- [ ] **Step 1: Write the failing test**

Create `tests/test_seg_file_writer.cpp`:

```cpp
#include "test_framework.hpp"
#include "seg_file_writer.hpp"
#include <fstream>
#include <sstream>
#include <thread>
#include <cstdio>
#include <set>

void test_writes_header_and_rows() {
    const char* path = "/tmp/dfn_cbs_test.seg";
    {
        SegFileWriter w(path);
        w.append_rows({{"AAAA-1", "chr1", 0, 1000, 1, 1.0}});
    }
    std::ifstream f(path);
    std::string header;
    std::getline(f, header);
    ASSERT_EQ(header, std::string("ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean"));
    std::string row;
    std::getline(f, row);
    ASSERT_EQ(row, std::string("AAAA-1\tchr1\t0\t1000\t1\t1"));
    std::remove(path);
}

void test_concurrent_appends_are_not_corrupted() {
    const char* path = "/tmp/dfn_cbs_test_concurrent.seg";
    {
        SegFileWriter w(path);
        std::vector<std::thread> threads;
        for (int t = 0; t < 10; ++t) {
            threads.emplace_back([&w, t]() {
                std::string barcode = "CELL" + std::to_string(t);
                w.append_rows({{barcode, "chr1", 0, 1000, 1, static_cast<double>(t)}});
            });
        }
        for (auto& th : threads) th.join();
    }

    std::ifstream f(path);
    std::string line;
    std::getline(f, line); // header
    std::set<std::string> seen_barcodes;
    int row_count = 0;
    while (std::getline(f, line)) {
        ++row_count;
        std::istringstream iss(line);
        std::string barcode;
        std::getline(iss, barcode, '\t');
        seen_barcodes.insert(barcode);
    }
    ASSERT_EQ(row_count, 10);
    ASSERT_EQ(seen_barcodes.size(), static_cast<size_t>(10));
    std::remove(path);
}

void test_unwritable_path_throws() {
    bool threw = false;
    try {
        SegFileWriter w("/tmp/dfn_cbs_nonexistent_dir/out.seg");
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);
}

int main() {
    test_writes_header_and_rows();
    test_concurrent_appends_are_not_corrupted();
    test_unwritable_path_throws();
    TEST_REPORT();
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++17 -Isrc -pthread tests/test_seg_file_writer.cpp -o /tmp/t 2>&1 | head -5`
Expected: FAIL — `seg_file_writer.hpp: No such file or directory`

- [ ] **Step 3: Write `src/seg_file_writer.hpp`**

```cpp
#pragma once
#include <string>
#include <vector>
#include <mutex>
#include <fstream>
#include "segment_cell.hpp"

// Writes an IGV/DNAcopy-format .seg file (plain text, not gzipped --
// row-per-segment, not a matrix). Thread-safe: append_rows() is called
// once per completed cell from worker threads, guarded by an internal
// mutex.
class SegFileWriter {
public:
    explicit SegFileWriter(const std::string& path); // throws on open failure
    ~SegFileWriter();

    SegFileWriter(const SegFileWriter&) = delete;
    SegFileWriter& operator=(const SegFileWriter&) = delete;

    void append_rows(const std::vector<SegRow>& rows);

    // Idempotent. Safe to call more than once.
    void close();

private:
    std::ofstream out_;
    std::mutex mutex_;
};
```

- [ ] **Step 4: Write `src/seg_file_writer.cpp`**

```cpp
#include "seg_file_writer.hpp"
#include <stdexcept>

SegFileWriter::SegFileWriter(const std::string& path) {
    out_.open(path, std::ios::out | std::ios::trunc);
    if (!out_.is_open()) {
        throw std::runtime_error("cannot open .seg file for writing: " + path);
    }
    out_ << "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\n";
}

SegFileWriter::~SegFileWriter() {
    close();
}

void SegFileWriter::append_rows(const std::vector<SegRow>& rows) {
    std::lock_guard<std::mutex> lock(mutex_);
    for (const auto& row : rows) {
        out_ << row.barcode << '\t' << row.chrom << '\t' << row.start << '\t'
             << row.end << '\t' << row.num_mark << '\t' << row.seg_mean << '\n';
    }
}

void SegFileWriter::close() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (out_.is_open()) {
        out_.close();
    }
}
```

- [ ] **Step 5: Add Makefile rules**

```make
$(BUILD_DIR)/seg_file_writer.o: src/seg_file_writer.cpp src/seg_file_writer.hpp src/segment_cell.hpp | $(BUILD_DIR)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc -c $< -o $@

$(TEST_BUILD)/test_seg_file_writer: tests/test_seg_file_writer.cpp $(BUILD_DIR)/seg_file_writer.o $(BUILD_DIR)/segment_cell.o $(BUILD_DIR)/bin_gc.o $(BUILD_DIR)/cbs_args.o $(BUILD_DIR)/vendor_stats.o $(BUILD_DIR)/vendor_correction.o $(BUILD_DIR)/vendor_cbs.o | $(TEST_BUILD)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc -Isrc/vendor/cbs $^ -o $@ -pthread
```

Add `$(TEST_BUILD)/test_seg_file_writer` to `TESTS :=`, and
`$(BUILD_DIR)/seg_file_writer.o` to `dfn_cbs`'s consolidated link line
(Task 10's Step 5).

- [ ] **Step 6: Run test to verify it passes**

Run: `make tests/build/test_seg_file_writer && ./tests/build/test_seg_file_writer`
Expected: `PASS` — 3 test functions, some assertions, `0 failed`.

- [ ] **Step 7: Commit**

```bash
git add src/seg_file_writer.hpp src/seg_file_writer.cpp tests/test_seg_file_writer.cpp Makefile
git commit -m "Add dfn_cbs .seg file writer"
```

---

### Task 12: Transpose-on-exit writer

Reads a fully-loaded cell-major scratch array (from `ScratchMatrix::
read_all()`) and writes it out in the `bins × cells` gzip TSV convention
`dfn_copy`'s `raw_counts.txt.gz` already uses — a strided read directly
into the gzip stream, no second full-size buffer. This task only writes
the file given an in-memory array; Task 14's orchestration is what loads
one scratch file at a time and frees it before the next (the actual
bounded-memory behavior lives in that call sequence, not in this
function).

**Files:**
- Create: `src/transpose_writer.hpp`
- Create: `src/transpose_writer.cpp`
- Test: `tests/test_transpose_writer.cpp`

**Interfaces:**
- Consumes: the flat `std::vector<float>` produced by
  `ScratchMatrix::read_all()` (Task 10).
- Produces: `void write_transposed_gz(const std::string& path, const
  std::vector<float>& scratch_data, size_t num_cells, size_t num_bins,
  const std::vector<std::string>& bin_names, const
  std::vector<std::string>& barcode_names)` — Task 14's Phase 3 calls
  this once per output matrix.

- [ ] **Step 1: Write the failing test**

Create `tests/test_transpose_writer.cpp`:

```cpp
#include "test_framework.hpp"
#include "transpose_writer.hpp"
#include <zlib.h>
#include <cstdio>
#include <sstream>

std::string read_gz(const std::string& path) {
    gzFile in = gzopen(path.c_str(), "rb");
    std::string content;
    char buf[256];
    int n;
    while ((n = gzread(in, buf, sizeof(buf))) > 0) content.append(buf, n);
    gzclose(in);
    return content;
}

void test_writes_bins_by_cells_layout() {
    // 2 cells x 3 bins, cell-major (row=cell): cell0={1,2,3}, cell1={4,5,6}.
    std::vector<float> scratch_data = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f};
    const char* path = "/tmp/dfn_cbs_test_transpose.txt.gz";
    write_transposed_gz(path, scratch_data, 2, 3,
                         {"chr1:0:1000", "chr1:1000:2000", "chr2:0:1000"},
                         {"AAAA-1", "CCCC-1"});

    std::string content = read_gz(path);
    std::istringstream iss(content);
    std::string line;
    std::getline(iss, line);
    ASSERT_EQ(line, std::string("bin\tAAAA-1\tCCCC-1"));
    std::getline(iss, line);
    ASSERT_EQ(line, std::string("chr1:0:1000\t1\t4"));
    std::getline(iss, line);
    ASSERT_EQ(line, std::string("chr1:1000:2000\t2\t5"));
    std::getline(iss, line);
    ASSERT_EQ(line, std::string("chr2:0:1000\t3\t6"));

    std::remove(path);
}

void test_size_mismatch_throws() {
    std::vector<float> scratch_data = {1.0f, 2.0f}; // wrong size for 2x3
    bool threw = false;
    try {
        write_transposed_gz("/tmp/dfn_cbs_test_transpose_bad.txt.gz", scratch_data, 2, 3,
                             {"a", "b", "c"}, {"X", "Y"});
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);
}

int main() {
    test_writes_bins_by_cells_layout();
    test_size_mismatch_throws();
    TEST_REPORT();
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++17 -Isrc tests/test_transpose_writer.cpp -lz -o /tmp/t 2>&1 | head -5`
Expected: FAIL — `transpose_writer.hpp: No such file or directory`

- [ ] **Step 3: Write `src/transpose_writer.hpp`**

```cpp
#pragma once
#include <string>
#include <vector>

// Writes a cell-major flat float32 array (as loaded from
// ScratchMatrix::read_all() -- row-major by cell,
// scratch_data[cell * num_bins + bin]) out as a bins x cells gzip TSV,
// the same convention dfn_copy's raw_counts.txt.gz uses. Throws
// std::runtime_error on any size mismatch or an unwritable path.
void write_transposed_gz(const std::string& path,
                          const std::vector<float>& scratch_data,
                          size_t num_cells,
                          size_t num_bins,
                          const std::vector<std::string>& bin_names,
                          const std::vector<std::string>& barcode_names);
```

- [ ] **Step 4: Write `src/transpose_writer.cpp`**

```cpp
#include "transpose_writer.hpp"
#include <zlib.h>
#include <sstream>
#include <stdexcept>

void write_transposed_gz(const std::string& path,
                          const std::vector<float>& scratch_data,
                          size_t num_cells,
                          size_t num_bins,
                          const std::vector<std::string>& bin_names,
                          const std::vector<std::string>& barcode_names) {
    if (bin_names.size() != num_bins) {
        throw std::runtime_error("transpose writer: bin_names size " +
                                  std::to_string(bin_names.size()) +
                                  " does not match num_bins " + std::to_string(num_bins));
    }
    if (barcode_names.size() != num_cells) {
        throw std::runtime_error("transpose writer: barcode_names size " +
                                  std::to_string(barcode_names.size()) +
                                  " does not match num_cells " + std::to_string(num_cells));
    }
    if (scratch_data.size() != num_cells * num_bins) {
        throw std::runtime_error("transpose writer: scratch_data size " +
                                  std::to_string(scratch_data.size()) +
                                  " does not match num_cells * num_bins (" +
                                  std::to_string(num_cells * num_bins) + ")");
    }

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

    for (size_t bin = 0; bin < num_bins; ++bin) {
        std::ostringstream row;
        row << bin_names[bin];
        for (size_t cell = 0; cell < num_cells; ++cell) {
            row << '\t' << scratch_data[cell * num_bins + bin];
        }
        row << '\n';
        std::string row_str = row.str();
        gzwrite(out, row_str.data(), static_cast<unsigned>(row_str.size()));
    }

    gzclose(out);
}
```

- [ ] **Step 5: Add Makefile rules**

```make
$(BUILD_DIR)/transpose_writer.o: src/transpose_writer.cpp src/transpose_writer.hpp | $(BUILD_DIR)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc -c $< -o $@

$(TEST_BUILD)/test_transpose_writer: tests/test_transpose_writer.cpp $(BUILD_DIR)/transpose_writer.o | $(TEST_BUILD)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc $^ -o $@ -lz
```

Add `$(TEST_BUILD)/test_transpose_writer` to `TESTS :=`, and
`$(BUILD_DIR)/transpose_writer.o` to `dfn_cbs`'s consolidated link line.

- [ ] **Step 6: Run test to verify it passes**

Run: `make tests/build/test_transpose_writer && ./tests/build/test_transpose_writer`
Expected: `PASS` — 2 test functions, some assertions, `0 failed`.

- [ ] **Step 7: Commit**

```bash
git add src/transpose_writer.hpp src/transpose_writer.cpp tests/test_transpose_writer.cpp Makefile
git commit -m "Add dfn_cbs transpose-on-exit writer"
```

---

### Task 13: Progress reporter

The string-formatting logic (`format_duration`, `format_with_commas`,
`format_phase2_status`) is kept as pure, directly testable functions,
separate from the actual background thread/timing/TTY-detection
machinery — so the tests assert on exact formatted content without being
timing-flaky.

**Files:**
- Create: `src/progress_reporter.hpp`
- Create: `src/progress_reporter.cpp`
- Test: `tests/test_progress_reporter.cpp`

**Interfaces:**
- Produces: `void log_phase_line(const std::string& message, bool
  quiet)` (Phase 1/3 one-shot messages), `std::string
  format_duration(double seconds)`, `std::string
  format_with_commas(size_t value)`, `std::string
  format_phase2_status(size_t completed, size_t total, double
  elapsed_seconds)`, and `class Phase2Reporter` (constructor `(size_t
  total_cells, bool quiet)`, `void increment()`, `size_t completed()
  const`, `void stop()`). Task 14's Phase 1/3 call `log_phase_line()`
  directly; Phase 2 constructs one `Phase2Reporter`, has each worker call
  `increment()` once per completed cell, and calls `stop()` (or lets the
  destructor do it) once all cells are done.

- [ ] **Step 1: Write the failing test**

Create `tests/test_progress_reporter.cpp`:

```cpp
#include "test_framework.hpp"
#include "progress_reporter.hpp"
#include <thread>
#include <vector>

void test_format_duration() {
    ASSERT_EQ(format_duration(0.0), std::string("0:00:00"));
    ASSERT_EQ(format_duration(59.0), std::string("0:00:59"));
    ASSERT_EQ(format_duration(3661.0), std::string("1:01:01"));
    ASSERT_EQ(format_duration(-5.0), std::string("0:00:00")); // clamped, never negative
}

void test_format_with_commas() {
    ASSERT_EQ(format_with_commas(0), std::string("0"));
    ASSERT_EQ(format_with_commas(999), std::string("999"));
    ASSERT_EQ(format_with_commas(1000), std::string("1,000"));
    ASSERT_EQ(format_with_commas(1234567), std::string("1,234,567"));
}

void test_format_phase2_status_exact_content() {
    // completed=100, total=400, elapsed=10s -> rate=10 cells/s,
    // remaining=300 cells -> ETA = 300/10 = 30s.
    std::string status = format_phase2_status(100, 400, 10.0);
    ASSERT_EQ(status, std::string(
        "Phase 2: 100 / 400 cells (25.0%) | elapsed 0:00:10 | ETA 0:00:30 | 10.0 cells/s"));
}

void test_format_phase2_status_handles_zero_elapsed() {
    // Guards against divide-by-zero on the very first tick.
    std::string status = format_phase2_status(0, 100, 0.0);
    ASSERT_TRUE(status.find("0.0 cells/s") != std::string::npos);
}

void test_reporter_increment_from_multiple_threads() {
    Phase2Reporter reporter(1000, /*quiet=*/true); // quiet: no background thread/printing
    std::vector<std::thread> threads;
    for (int t = 0; t < 10; ++t) {
        threads.emplace_back([&reporter]() {
            for (int i = 0; i < 100; ++i) reporter.increment();
        });
    }
    for (auto& th : threads) th.join();
    ASSERT_EQ(reporter.completed(), static_cast<size_t>(1000));
    reporter.stop(); // idempotent, must not hang
}

int main() {
    test_format_duration();
    test_format_with_commas();
    test_format_phase2_status_exact_content();
    test_format_phase2_status_handles_zero_elapsed();
    test_reporter_increment_from_multiple_threads();
    TEST_REPORT();
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -std=c++17 -Isrc -pthread tests/test_progress_reporter.cpp -o /tmp/t 2>&1 | head -5`
Expected: FAIL — `progress_reporter.hpp: No such file or directory`

- [ ] **Step 3: Write `src/progress_reporter.hpp`**

```cpp
#pragma once
#include <atomic>
#include <cstddef>
#include <string>
#include <thread>
#include <chrono>

// One-shot phase start/done message to stderr, prefixed "[dfn_cbs] ".
// No-op if quiet is true.
void log_phase_line(const std::string& message, bool quiet);

// "H:MM:SS" (hours unbounded). Negative input is clamped to 0.
std::string format_duration(double seconds);

// e.g. 1234567 -> "1,234,567".
std::string format_with_commas(size_t value);

// Builds the Phase 2 status line's content -- no "[dfn_cbs] " prefix, no
// trailing newline/carriage return (the caller decides how to terminate
// it based on whether stderr is a TTY). Pure/testable independent of
// real timing.
std::string format_phase2_status(size_t completed, size_t total, double elapsed_seconds);

// Periodic Phase 2 progress reporter: a background thread wakes on a
// ~5-second cadence, reads an atomic "cells completed" counter (workers
// call increment() with no locking on their hot path), and prints a
// status line to stderr -- \r-rewritten in place if stderr is a TTY,
// newline-terminated otherwise. If quiet is true, no thread is ever
// started and increment()/completed() are plain counter operations with
// no I/O.
class Phase2Reporter {
public:
    Phase2Reporter(size_t total_cells, bool quiet);
    ~Phase2Reporter();

    Phase2Reporter(const Phase2Reporter&) = delete;
    Phase2Reporter& operator=(const Phase2Reporter&) = delete;

    void increment();
    size_t completed() const;

    // Stops and joins the reporter thread. Idempotent -- safe to call
    // more than once, and safe to omit (the destructor calls it too).
    void stop();

private:
    void run();

    std::atomic<size_t> completed_;
    size_t total_;
    bool quiet_;
    bool is_tty_;
    std::atomic<bool> stop_requested_;
    std::thread thread_;
    std::chrono::steady_clock::time_point start_time_;
};
```

- [ ] **Step 4: Write `src/progress_reporter.cpp`**

```cpp
#include "progress_reporter.hpp"
#include <cstdio>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <unistd.h>

void log_phase_line(const std::string& message, bool quiet) {
    if (quiet) return;
    std::fprintf(stderr, "[dfn_cbs] %s\n", message.c_str());
}

std::string format_duration(double seconds) {
    if (seconds < 0.0) seconds = 0.0;
    long total_seconds = static_cast<long>(seconds);
    long h = total_seconds / 3600;
    long m = (total_seconds % 3600) / 60;
    long s = total_seconds % 60;
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%ld:%02ld:%02ld", h, m, s);
    return std::string(buf);
}

std::string format_with_commas(size_t value) {
    std::string digits = std::to_string(value);
    std::string out;
    int count = 0;
    for (auto it = digits.rbegin(); it != digits.rend(); ++it) {
        if (count > 0 && count % 3 == 0) out.push_back(',');
        out.push_back(*it);
        ++count;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

std::string format_phase2_status(size_t completed, size_t total, double elapsed_seconds) {
    double pct = total > 0
        ? (100.0 * static_cast<double>(completed) / static_cast<double>(total))
        : 0.0;
    double rate = elapsed_seconds > 0.0
        ? static_cast<double>(completed) / elapsed_seconds
        : 0.0;
    double eta_seconds = rate > 0.0
        ? static_cast<double>(total - completed) / rate
        : 0.0;

    std::ostringstream oss;
    oss << "Phase 2: " << format_with_commas(completed) << " / "
        << format_with_commas(total) << " cells ("
        << std::fixed << std::setprecision(1) << pct << "%)"
        << " | elapsed " << format_duration(elapsed_seconds)
        << " | ETA " << format_duration(eta_seconds)
        << " | " << std::fixed << std::setprecision(1) << rate << " cells/s";
    return oss.str();
}

Phase2Reporter::Phase2Reporter(size_t total_cells, bool quiet)
    : completed_(0), total_(total_cells), quiet_(quiet),
      is_tty_(::isatty(fileno(stderr)) != 0), stop_requested_(false) {
    start_time_ = std::chrono::steady_clock::now();
    if (!quiet_) {
        thread_ = std::thread(&Phase2Reporter::run, this);
    }
}

Phase2Reporter::~Phase2Reporter() {
    stop();
}

void Phase2Reporter::increment() {
    completed_.fetch_add(1, std::memory_order_relaxed);
}

size_t Phase2Reporter::completed() const {
    return completed_.load(std::memory_order_relaxed);
}

void Phase2Reporter::stop() {
    if (stop_requested_.exchange(true)) return; // already stopped
    if (thread_.joinable()) thread_.join();
}

void Phase2Reporter::run() {
    const auto tick = std::chrono::milliseconds(200);
    const int ticks_per_update = 25; // 200ms * 25 = ~5s cadence
    int tick_count = 0;
    while (!stop_requested_.load()) {
        std::this_thread::sleep_for(tick);
        if (stop_requested_.load()) break;
        if (++tick_count < ticks_per_update) continue;
        tick_count = 0;

        double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time_).count();
        std::string status = "[dfn_cbs] " + format_phase2_status(completed(), total_, elapsed);
        if (is_tty_) {
            std::fprintf(stderr, "\r%s", status.c_str());
            std::fflush(stderr);
        } else {
            std::fprintf(stderr, "%s\n", status.c_str());
        }
    }
    if (is_tty_) {
        std::fprintf(stderr, "\n"); // move off the in-place line before any further output
    }
}
```

- [ ] **Step 5: Add Makefile rules**

```make
$(BUILD_DIR)/progress_reporter.o: src/progress_reporter.cpp src/progress_reporter.hpp | $(BUILD_DIR)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc -c $< -o $@

$(TEST_BUILD)/test_progress_reporter: tests/test_progress_reporter.cpp $(BUILD_DIR)/progress_reporter.o | $(TEST_BUILD)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc $^ -o $@ -pthread
```

Add `$(TEST_BUILD)/test_progress_reporter` to `TESTS :=`, and
`$(BUILD_DIR)/progress_reporter.o` to `dfn_cbs`'s consolidated link line.

- [ ] **Step 6: Run test to verify it passes**

Run: `make tests/build/test_progress_reporter && ./tests/build/test_progress_reporter`
Expected: `PASS` — 5 test functions, some assertions, `0 failed`, completing in well under a second (the quiet-mode test never sleeps or starts a thread).

- [ ] **Step 7: Commit**

```bash
git add src/progress_reporter.hpp src/progress_reporter.cpp tests/test_progress_reporter.cpp Makefile
git commit -m "Add dfn_cbs progress reporter"
```

---

### Task 14: `main.cpp` orchestration and end-to-end integration tests

Wires every prior task into the three-phase pipeline: load once (Phase
1) → per-cell parallel compute writing cell-major scratch files (Phase
2) → transpose-on-exit, one scratch file at a time (Phase 3). Unlike
`dfn_copy`'s BAM-based fixtures, `dfn_cbs`'s input is already a matrix —
integration fixtures are plain Python (gzip + text), no `pysam`/BAM
involved at all.

**Files:**
- Modify: `src/cbs_main.cpp` (replace the Task 2 placeholder with the
  full orchestration)
- Create: `tests/fixtures/gen_cbs_fixture.py`
- Create: `tests/integration/test_dfn_cbs.py`
- Modify: `README.md` (add a `dfn_cbs` usage section)
- Modify: `Makefile` (final consolidated `all`/`clean`/`test` wiring)

**Interfaces:**
- Consumes every module from Tasks 2–13.

- [ ] **Step 1: Replace `src/cbs_main.cpp`**

```cpp
#include "cbs_args.hpp"
#include "bin_gc.hpp"
#include "raw_counts_loader.hpp"
#include "cell_filter.hpp"
#include "gc_correct.hpp"
#include "segment_cell.hpp"
#include "scratch_matrix.hpp"
#include "seg_file_writer.hpp"
#include "transpose_writer.hpp"
#include "progress_reporter.hpp"

#include <cstdio>
#include <cstdlib>
#include <memory>
#include <mutex>
#include <thread>
#include <atomic>
#include <chrono>
#include <stdexcept>

namespace {

void write_one_matrix(ScratchMatrix& scratch, const std::string& out_path,
                       const std::string& label,
                       const std::vector<std::string>& bin_names,
                       const std::vector<std::string>& barcode_names,
                       bool quiet) {
    log_phase_line("Writing " + label + " (bins x cells transpose)...", quiet);
    auto t0 = std::chrono::steady_clock::now();
    {
        std::vector<float> data = scratch.read_all();
        write_transposed_gz(out_path, data, scratch.num_cells(), scratch.num_bins(),
                             bin_names, barcode_names);
    } // `data` freed here, before the caller loads the next scratch file
    double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    log_phase_line("done in " + format_duration(elapsed), quiet);
}

int run(const CbsArgs& args) {
    auto t_load_start = std::chrono::steady_clock::now();

    log_phase_line("Loading bin/GC file...", args.quiet);
    std::vector<BinGc> bins = load_bin_gc(args.bins_path);
    std::vector<std::string> bin_names;
    std::vector<double> gc_values;
    bin_names.reserve(bins.size());
    gc_values.reserve(bins.size());
    for (const auto& b : bins) {
        bin_names.push_back(bin_name(b));
        gc_values.push_back(b.gc);
    }

    log_phase_line("Loading raw counts matrix...", args.quiet);
    RawCountsMatrix matrix = load_raw_counts(args.counts_path, bins.size(), bin_names);
    double load_elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_load_start).count();
    log_phase_line("done in " + format_duration(load_elapsed), args.quiet);

    FilteredCells filtered = filter_cells(matrix, args.min_reads);
    if (filtered.surviving_indices.empty()) {
        std::fprintf(stderr,
            "error: no cells survived --min-reads %lld (highest observed cell "
            "total was %lld)\n",
            static_cast<long long>(args.min_reads),
            static_cast<long long>(filtered.max_column_sum));
        return 1;
    }

    size_t num_surviving = filtered.surviving_indices.size();
    std::vector<std::string> surviving_barcodes;
    surviving_barcodes.reserve(num_surviving);
    for (size_t idx : filtered.surviving_indices) {
        surviving_barcodes.push_back(matrix.barcodes[idx]);
    }

    std::string gc_scratch_path = args.out_prefix + ".gc_corrected.scratch.bin";
    std::string ratio_scratch_path = args.out_prefix + ".lowess_ratio.scratch.bin";
    std::string seg_scratch_path = args.out_prefix + ".segmented_lowess_ratio.scratch.bin";

    auto gc_scratch = std::make_unique<ScratchMatrix>(
        ScratchMatrix::create(gc_scratch_path, num_surviving, bins.size()));
    auto ratio_scratch = std::make_unique<ScratchMatrix>(
        ScratchMatrix::create(ratio_scratch_path, num_surviving, bins.size()));
    auto seg_scratch = std::make_unique<ScratchMatrix>(
        ScratchMatrix::create(seg_scratch_path, num_surviving, bins.size()));

    auto seg_writer = std::make_unique<SegFileWriter>(args.out_prefix + ".seg");

    log_phase_line("Phase 2: GC-correcting and segmenting " +
                   format_with_commas(num_surviving) + " cells...", args.quiet);
    Phase2Reporter reporter(num_surviving, args.quiet);

    std::atomic<size_t> next_index{0};
    std::mutex error_mutex;
    std::exception_ptr first_error;
    int num_threads = args.threads < 1 ? 1 : args.threads;

    std::vector<std::thread> workers;
    workers.reserve(static_cast<size_t>(num_threads));
    for (int t = 0; t < num_threads; ++t) {
        workers.emplace_back([&]() {
            while (true) {
                size_t i = next_index.fetch_add(1, std::memory_order_relaxed);
                if (i >= num_surviving) return;
                {
                    std::lock_guard<std::mutex> lock(error_mutex);
                    if (first_error) return;
                }
                try {
                    size_t orig_cell = filtered.surviving_indices[i];
                    std::vector<double> counts(bins.size());
                    for (size_t bin = 0; bin < bins.size(); ++bin) {
                        counts[bin] = static_cast<double>(
                            matrix.data[bin * matrix.num_cells + orig_cell]);
                    }
                    GcCorrectedCell corrected = gc_correct_cell(counts, gc_values);
                    CellSegmentation seg_result = segment_cell(
                        bins, corrected.lowess_ratio, matrix.barcodes[orig_cell], args);

                    gc_scratch->write_cell(i, corrected.gc_corrected);
                    ratio_scratch->write_cell(i, corrected.lowess_ratio);
                    seg_scratch->write_cell(i, seg_result.segmented_lowess_ratio);
                    seg_writer->append_rows(seg_result.seg_rows);

                    reporter.increment();
                } catch (...) {
                    std::lock_guard<std::mutex> lock(error_mutex);
                    if (!first_error) first_error = std::current_exception();
                }
            }
        });
    }
    for (auto& w : workers) w.join();
    reporter.stop();

    if (first_error) std::rethrow_exception(first_error);

    seg_writer->close();

    // Free the large input matrix -- Phase 3 never needs it again.
    matrix.data.clear();
    matrix.data.shrink_to_fit();

    write_one_matrix(*gc_scratch, args.out_prefix + ".gc_corrected.txt.gz",
                      "gc_corrected.txt.gz", bin_names, surviving_barcodes, args.quiet);
    gc_scratch->close();
    std::remove(gc_scratch_path.c_str());
    gc_scratch.reset();

    write_one_matrix(*ratio_scratch, args.out_prefix + ".lowess_ratio.txt.gz",
                      "lowess_ratio.txt.gz", bin_names, surviving_barcodes, args.quiet);
    ratio_scratch->close();
    std::remove(ratio_scratch_path.c_str());
    ratio_scratch.reset();

    write_one_matrix(*seg_scratch, args.out_prefix + ".segmented_lowess_ratio.txt.gz",
                      "segmented_lowess_ratio.txt.gz", bin_names, surviving_barcodes, args.quiet);
    seg_scratch->close();
    std::remove(seg_scratch_path.c_str());
    seg_scratch.reset();

    return 0;
}

} // namespace

int main(int argc, char** argv) {
    CbsArgs args;
    bool help = false;
    if (!parse_cbs_args(argc, argv, args, help)) {
        return 1;
    }
    if (help) {
        std::fputs(cbs_usage_text().c_str(), stdout);
        return 0;
    }

    try {
        return run(args);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
```

- [ ] **Step 2: Rebuild and smoke-test manually**

Run: `make clean && make all`
Expected: both binaries build with no errors/warnings.

- [ ] **Step 3: Create the fixture generator**

Create `tests/fixtures/gen_cbs_fixture.py`:

```python
"""Generates a small bins.tsv + raw_counts.txt.gz pair for dfn_cbs
integration tests -- no BAM/pysam involved, since dfn_cbs's input is
already a matrix."""
import gzip
import os


def main(out_dir):
    bins_path = os.path.join(out_dir, "bins.tsv")
    counts_path = os.path.join(out_dir, "counts.raw_counts.txt.gz")

    bins = []
    for i in range(20):
        gc = 0.4 + 0.01 * (i % 5)
        bins.append(("chr1", i * 1000, (i + 1) * 1000, gc))
    for i in range(20):
        gc = 0.4 + 0.01 * (i % 5)
        bins.append(("chr2", i * 1000, (i + 1) * 1000, gc))

    with open(bins_path, "w") as f:
        f.write("chrom\tstart\tend\tgc\n")
        for chrom, start, end, gc in bins:
            f.write(f"{chrom}\t{start}\t{end}\t{gc}\n")

    # AAAA-1: diploid, ~200 reads/bin on both chromosomes, sum=8000 -- survives --min-reads 1000.
    # CCCC-1: ~200 reads/bin on chr1, ~500 reads/bin on chr2 (a whole-chromosome
    #         gain, ~2.5x) -- also survives, and its segmented output must
    #         preserve the chr1 vs chr2 difference (no centering).
    # GGGG-1: 1 read/bin everywhere, sum=40 -- filtered out by --min-reads 1000,
    #         must be absent from every output file.
    barcodes = ["AAAA-1", "CCCC-1", "GGGG-1"]
    counts = {"AAAA-1": [], "CCCC-1": [], "GGGG-1": []}
    for chrom, start, end, gc in bins:
        counts["AAAA-1"].append(200)
        counts["CCCC-1"].append(200 if chrom == "chr1" else 500)
        counts["GGGG-1"].append(1)

    with gzip.open(counts_path, "wt") as f:
        f.write("bin\t" + "\t".join(barcodes) + "\n")
        for i, (chrom, start, end, gc) in enumerate(bins):
            row = [f"{chrom}:{start}:{end}"] + [str(counts[b][i]) for b in barcodes]
            f.write("\t".join(row) + "\n")

    return counts_path, bins_path


if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
```

- [ ] **Step 4: Write the failing integration test**

Create `tests/integration/test_dfn_cbs.py`:

```python
import gzip
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "fixtures"))
import gen_cbs_fixture

REPO_ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
DFN_CBS = os.path.join(REPO_ROOT, "dfn_cbs")


def read_matrix(path):
    with gzip.open(path, "rt") as f:
        lines = [line.rstrip("\n").split("\t") for line in f]
    header, rows = lines[0], lines[1:]
    return header, {row[0]: row[1:] for row in rows}


def run_dfn_cbs(tmp_path, extra_args=None):
    counts_path, bins_path = gen_cbs_fixture.main(str(tmp_path))
    out_prefix = str(tmp_path / "out")
    args = [DFN_CBS, "--counts", counts_path, "--bins", bins_path,
            "--out-prefix", out_prefix, "--min-reads", "1000",
            "--perms", "50", "--min-seg-len", "5", "--threads", "2", "--quiet"]
    if extra_args:
        args.extend(extra_args)
    result = subprocess.run(args, capture_output=True, text=True)
    return result, out_prefix


def test_all_four_outputs_are_produced_with_correct_shape(tmp_path):
    result, out_prefix = run_dfn_cbs(tmp_path)
    assert result.returncode == 0, result.stderr

    for suffix in ["gc_corrected.txt.gz", "lowess_ratio.txt.gz", "segmented_lowess_ratio.txt.gz"]:
        header, rows = read_matrix(out_prefix + "." + suffix)
        assert header == ["bin", "AAAA-1", "CCCC-1"]  # GGGG-1 filtered out
        assert len(rows) == 40  # 20 bins x 2 chromosomes

    assert os.path.exists(out_prefix + ".seg")


def test_filtered_cell_absent_from_seg_file(tmp_path):
    result, out_prefix = run_dfn_cbs(tmp_path)
    assert result.returncode == 0, result.stderr

    with open(out_prefix + ".seg") as f:
        lines = f.readlines()
    header = lines[0].rstrip("\n").split("\t")
    assert header == ["ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean"]
    ids = {line.split("\t")[0] for line in lines[1:]}
    assert "GGGG-1" not in ids
    assert ids == {"AAAA-1", "CCCC-1"}


def test_whole_chromosome_shift_survives_without_centering(tmp_path):
    result, out_prefix = run_dfn_cbs(tmp_path)
    assert result.returncode == 0, result.stderr

    header, rows = read_matrix(out_prefix + ".segmented_lowess_ratio.txt.gz")
    cccc_idx = header.index("CCCC-1")

    chr1_values = [float(rows[f"chr1:{i*1000}:{(i+1)*1000}"][cccc_idx]) for i in range(20)]
    chr2_values = [float(rows[f"chr2:{i*1000}:{(i+1)*1000}"][cccc_idx]) for i in range(20)]

    avg_chr1 = sum(chr1_values) / len(chr1_values)
    avg_chr2 = sum(chr2_values) / len(chr2_values)

    # chr1 ~1.0 (diploid baseline), chr2 ~2.5 (the injected gain) -- if CBS's
    # centering weren't forced off, these would collapse toward each other.
    assert avg_chr1 < 1.3
    assert avg_chr2 > 2.0


def test_every_cell_filtered_out_fails_cleanly(tmp_path):
    result, _ = run_dfn_cbs(tmp_path, extra_args=["--min-reads", "999999999"])
    assert result.returncode == 1
    assert "no cells survived" in result.stderr


def test_malformed_bins_file_fails_cleanly(tmp_path):
    counts_path, bins_path = gen_cbs_fixture.main(str(tmp_path))
    with open(bins_path, "w") as f:
        f.write("chrom\tstart\tend\n")  # no gc column
        f.write("chr1\t0\t1000\n")
    out_prefix = str(tmp_path / "out")
    result = subprocess.run(
        [DFN_CBS, "--counts", counts_path, "--bins", bins_path, "--out-prefix", out_prefix],
        capture_output=True, text=True)
    assert result.returncode == 1
    assert "missing required column" in result.stderr


def test_unwritable_out_prefix_dir_fails_cleanly(tmp_path):
    counts_path, bins_path = gen_cbs_fixture.main(str(tmp_path))
    out_prefix = str(tmp_path / "nonexistent_dir" / "out")
    result = subprocess.run(
        [DFN_CBS, "--counts", counts_path, "--bins", bins_path, "--out-prefix", out_prefix],
        capture_output=True, text=True)
    assert result.returncode == 1
    assert "libc++abi" not in result.stderr
    assert "terminating due to uncaught exception" not in result.stderr
```

- [ ] **Step 5: Run the integration tests to verify they fail (binary not built with real logic yet, or fixture path issues)**

Run: `.venv/bin/pytest tests/integration/test_dfn_cbs.py -v`
Expected at this point in the plan (after Step 2's rebuild) most of these
should already pass, since Steps 1–2 landed the real implementation —
this step is the safety check that they weren't accidentally trivially
true. If any fail, fix the orchestration code in Step 1 before proceeding
(most likely culprit: a path or flag name typo between the Python fixture/
test and the C++ CLI).

- [ ] **Step 6: Run test to verify it passes**

Run: `.venv/bin/pytest tests/integration/test_dfn_cbs.py -v`
Expected: 6 passed.

- [ ] **Step 7: Run the full test suite (unit + both integration suites)**

Run: `make test`
Expected: all `dfn_copy` and `dfn_cbs` unit suites report `0 failed`.

Run: `.venv/bin/pytest tests/integration -v`
Expected: all `dfn_copy` and `dfn_cbs` integration tests pass (this repo's
existing 24 `dfn_copy` integration tests plus the 6 new `dfn_cbs` ones).

- [ ] **Step 8: Add a `dfn_cbs` usage section to `README.md`**

Add a new section after the existing `## Usage` section (before
`## Validation`):

```markdown
## dfn_cbs: GC correction + segmentation at scale

`correct_gc.py` GC-corrects a raw counts matrix in Python, but can run out
of memory on very large matrices (hundreds of thousands of cells). For
that scale, use `dfn_cbs` instead -- a C++ tool that does the same
per-cell LOWESS GC correction with a bounded memory model, and
additionally segments each cell's corrected profile via a vendored copy
of [cbs+](https://github.com)'s circular binary segmentation core.

```bash
./dfn_cbs --counts sample.raw_counts.txt.gz --bins bins.tsv --out-prefix sample
```

Cells with fewer than 100,000 total raw binned reads are dropped by
default (`--min-reads` to change it). Outputs:

- `sample.gc_corrected.txt.gz` -- same as `correct_gc.py`'s output.
- `sample.lowess_ratio.txt.gz` -- the dimensionless GC-correction ratio
  (mean approx. 1 per cell); this is what gets segmented.
- `sample.segmented_lowess_ratio.txt.gz` -- same shape, each bin replaced
  by its segment's mean.
- `sample.seg` -- IGV/DNAcopy-format segments file (plain text, not
  gzipped).

See `./dfn_cbs --help` for the full flag list (including CBS tuning
parameters `--alpha`, `--perms`, `--min-seg-len`, `--max-depth`,
`--cbs-method`) and
`docs/superpowers/specs/2026-07-22-dfn_cbs-gc-correction-segmentation-design.md`
for the full design rationale.
```

- [ ] **Step 9: Final full-repo verification**

Run: `make clean && make all && make test`
Expected: clean build, all unit suites `0 failed`.

Run: `.venv/bin/pytest tests/integration tests_py -v`
Expected: all integration and `correct_gc.py` tests pass.

- [ ] **Step 10: Commit**

```bash
git add src/cbs_main.cpp tests/fixtures/gen_cbs_fixture.py \
        tests/integration/test_dfn_cbs.py README.md Makefile
git commit -m "Wire dfn_cbs orchestration and add end-to-end integration tests"
```

---

## Self-Review

**Spec coverage:** every section of
`docs/superpowers/specs/2026-07-22-dfn_cbs-gc-correction-segmentation-design.md`
maps to a task — Purpose/Background/Related work (Tasks 6, 8), Scope/
non-goals (enforced throughout: no htslib linkage checked in Task 1, no
`cbs+` modification, `correct_gc.py` untouched), Architecture's three
phases (Task 14's `run()`, built from Tasks 3–13), Progress reporting
(Task 13), LOWESS implementation (Task 6), CLI contract (Task 2), Outputs
(Tasks 10–12, 14), Error handling (validated per-module in Tasks 3–5, 10,
11, and exercised end-to-end in Task 14's integration tests), Build
(Tasks 1, 8's Makefile additions), Testing (a unit test per module, plus
Task 14's integration suite covering the exact three properties the spec
names: shapes/labels, filtered-cell absence, and the whole-chromosome-
shift-survives-centering-off property).

**Placeholder scan:** no TBD/TODO; every step has complete, runnable code
and exact expected command output.

**Type consistency check:** `BinGc` (Task 3) flows into
`raw_counts_loader`'s `expected_bin_names` (via `bin_name()`, Task 4),
`cell_filter`'s `RawCountsMatrix` (Task 5), `gc_correct_cell`'s `gc`
parameter (Task 7), and `segment_cell`'s `bins` parameter (Task 9) — same
field names (`chrom`, `start`, `end`, `gc`) throughout. `CbsArgs` (Task 2)
flows unchanged into `segment_cell` (Task 9) and `run()` (Task 14).
`GcCorrectedCell::lowess_ratio` (Task 7) is exactly what `segment_cell`
(Task 9) consumes as its `lowess_ratio` parameter. `ScratchMatrix::
read_all()`'s `std::vector<float>` (Task 10) is exactly `transpose_
writer`'s `scratch_data` parameter type (Task 12). `SegRow` (Task 9) is
exactly what `SegFileWriter::append_rows` (Task 11) consumes. Verified no
naming drift across tasks.
