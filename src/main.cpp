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
