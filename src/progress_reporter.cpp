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
