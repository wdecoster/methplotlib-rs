# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
make ci                 # fmt-check + clippy + test, what CI runs
make all                # fmt, clippy, test, build
cargo test              # unit tests, tests against the Python version, end-to-end runs
python scripts/make_test_data.py   # regenerate test-data/ (needs the Python methplotlib)
```

`test-data/` is synthetic (the annotation is from the Python version's examples).
`test-data/expected/*.json` is what the Python methplotlib parses from it, written by
`scripts/make_test_data.py`; `tests/python_reference.rs` compares the Rust parsers against it.
The repository is public: test data must stay synthetic or public.

Clippy is run as `cargo clippy --all-targets --all-features -- -D warnings`; warnings fail the build.

## What this is

A Rust port of the Python [methplotlib](https://github.com/wdecoster/methplotlib), which stays maintained.
The binary is named `methplotlib`, like the Python command.

Scope, as decided when starting the port:

- Input: bam/cram with MM/ML tags, bedMethyl (modkit and `modbam2bed --extended`), bedgraph; gtf/gff and bed annotation.
- Not ported: nanopolish and nanocompore formats, the QC report, the differential/allele-specific scripts.
- noodles for bam/cram and tabix, to stay pure Rust; kuva (`TrackStack`) for plotting.
- kuva has no zoom/pan. Output is static svg/png/pdf for now.

Plotting (`src/plot.rs`): each track is a kuva `Track`. `Panel` lets kuva draw the axes and title
through `render_multiple` (with an invisible scatter, as kuva only draws axes for a plot), and draws
the data itself in the coordinates of the same `ComputedLayout`, clipped to the plot area.

When porting behaviour, the Python version is the reference: match its results (parsed positions,
frequencies, smoothing, read-row assignment) unless a difference is a deliberate, documented decision.
