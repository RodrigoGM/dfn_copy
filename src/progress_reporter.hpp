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
