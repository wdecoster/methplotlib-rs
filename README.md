# methplotlib-rs

A Rust port of [methplotlib](https://github.com/wdecoster/methplotlib): plot nucleotide modifications
(e.g. 5mC and 5hmC) from long-read sequencing in a genomic window, together with gene and region annotation.

**Status: under construction.** The command-line interface is in place, plotting is not implemented yet.
Use the [Python version](https://github.com/wdecoster/methplotlib) in the meantime, which remains maintained.

## Scope

Input:

- bam/cram with MM/ML tags
- bedMethyl from `modkit pileup` (with or without `--only-tabs`) and from `modbam2bed --extended`
- bedgraph, including the 5-column output of `modkit pileup --bedgraph`
- annotation from gtf/gff and bed

Output: an svg, png or pdf genome-browser view, with a shared genomic x-axis across a track per dataset
and the annotation.

Compared to the Python version, the nanopolish and nanocompore formats, the QC report and the
`differential_modification` / `allele_specific_modification` scripts are not ported.

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
