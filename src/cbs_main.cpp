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
