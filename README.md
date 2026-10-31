# methplotlib-rs

A Rust port of [methplotlib](https://github.com/wdecoster/methplotlib): plot nucleotide modifications
(e.g. 5mC and 5hmC) from long-read sequencing in a genomic window, together with gene and region annotation.

**Status: early.** Reading all supported formats and plotting work, and the parsed data is tested against
the [Python version](https://github.com/wdecoster/methplotlib), which remains maintained.

**This is not a drop-in replacement for the Python methplotlib.** It makes static images rather than
interactive html, supports fewer input formats, has no QC report, and some options differ.
Both versions are maintained, for different uses.

## Rust or Python?

Use the **Python version** to explore a region: its html output can be zoomed and panned, and shows details
when hovering. It also reads nanopolish and nanocompore output, makes a QC report, and comes with the
`differential_modification` and `allele_specific_modification` scripts.

Use the **Rust version** to make figures, in particular many of them, e.g. in a pipeline:

- it is a single binary without dependencies, rather than a Python environment with pysam, pandas,
  pyranges and plotly (and kaleido with a browser for static images)
- it starts in a fraction of a second, where Python spends a second or more on imports before reading
  any data, which adds up when plotting hundreds of regions
- it uses less memory: for 1 Mb of 50x per-read data about 0.3-0.5 GB rather than 1 GB, which matters when
  running many in parallel on a shared server

It is not much faster for a single plot: on that same 1 Mb of 50x per-read data it takes 3-4 s, against
5-8 s for the current Python version, and both take about a second for a typical window.

| | Python | Rust |
|---|---|---|
| output | interactive html, optionally a static image | png (default), svg or pdf |
| zoom, pan, hover | yes | no |
| bam/cram, bedMethyl, bedgraph | yes | yes |
| nanopolish, nanocompore | yes | no |
| QC report, differential and allele-specific scripts | yes | no |
| gtf/gff and bed annotation, built-in MANE | yes | yes |
| installation | Python environment | single binary |

## Scope

Input:

- bam/cram with MM/ML tags
- bedMethyl from `modkit pileup` (with or without `--only-tabs`) and from `modbam2bed --extended`
- bedgraph, including the 5-column output of `modkit pileup --bedgraph`
- annotation from gtf/gff and bed, or `-g mane` for the [MANE Select](https://www.ncbi.nlm.nih.gov/refseq/MANE/)
  transcripts of GRCh38 (v1.5), downloaded once to `~/.cache/methplotlib` (or `$XDG_CACHE_HOME/methplotlib`)

Output: a png, svg or pdf genome-browser view, with a shared genomic x-axis across a track per dataset
and the annotation.

Differences in the options and output compared to the Python version:

- the format is set by the extension of `-o`: png (default), svg or pdf, rather than html;
  an svg of many per-read calls gets large, as every call is an element
- no `--static`, `--qcfile`, `--binary`, `--store` or `--example`
- transcript and gene names are written in the annotation track, as a static image has no hovering

## Design

- [noodles](https://github.com/zaeleus/noodles) to read bam/cram and tabix-indexed files, which keeps the
  tool pure Rust and easy to ship as a static binary
- [kuva](https://github.com/Psy-Fer/kuva) to plot, using its `TrackStack` genome-browser layout

## Installation

From source, with a [Rust toolchain](https://rustup.rs):

```bash
cargo install --git https://github.com/wdecoster/methplotlib-rs
```

The binary is called `methplotlib`, like the Python version's command. Don't install both in the same
environment.

## Usage

```text
methplotlib -m sample1.bam sample2.bedmethyl.gz -n sample1 sample2 -w chr7:5,525,542-5,543,028 -g genes.gtf.gz
```

Run `methplotlib --help` for all options.

## License

MIT
