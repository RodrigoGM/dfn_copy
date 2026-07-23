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
