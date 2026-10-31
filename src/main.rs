use anyhow::Result;
use clap::Parser;
use log::info;

mod region;

use region::Region;

/// Plot nucleotide modifications from long-read sequencing
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Modification data: bam/cram with MM/ML tags, bedMethyl (modkit or modbam2bed --extended) or bedgraph
    #[arg(short, long, num_args = 1.., required = true)]
    methylation: Vec<String>,

    /// Names of the datasets in --methylation, in the same order
    #[arg(short, long, num_args = 1.., required = true)]
    names: Vec<String>,

    /// Window (region) to plot, e.g. chr7:5,525,542-5,543,028, or an entire chromosome with --fasta
    #[arg(short, long)]
    window: String,

    /// Add annotation from a gtf or gff file
    #[arg(short, long)]
    gtf: Option<String>,

    /// Add annotation from a bed file
    #[arg(short, long)]
    bed: Option<String>,

    /// Reference fasta, required when --window is an entire chromosome and to decode cram
    #[arg(short, long)]
    fasta: Option<String>,

    /// Simplify the annotation track to show genes rather than transcripts
    #[arg(long)]
    simplify: bool,

    /// Split, rather than overlay, the modification tracks
    #[arg(long)]
    split: bool,

    /// Rolling window size for averaging frequency values, counted in records
    #[arg(long, default_value_t = 5)]
    smooth: usize,

    /// Comma separated list of modifications to restrict to, e.g. C+m,C+h
    #[arg(long)]
    mods: Option<String>,

    /// Size of the dots in the per read plots
    #[arg(long, default_value_t = 4)]
    dotsize: u32,

    /// Minimal modification probability to show, on the 0-255 scale of the ML tag
    #[arg(long, default_value_t = 20)]
    minqual: u8,

    /// Output file, its extension (svg, png or pdf) sets the format.
    /// Use {region} as a shorthand for {chr}_{start}_{end}.
    /// Default: methylation_browser_{region}.svg
    #[arg(short, long)]
    outfile: Option<String>,
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();
    if cli.names.len() != cli.methylation.len() {
        anyhow::bail!(
            "expecting the same number of --names ({}) as --methylation datasets ({})",
            cli.names.len(),
            cli.methylation.len()
        );
    }
    let region = Region::parse(&cli.window)?;
    info!("Plotting {region:?}");
    anyhow::bail!("methplotlib-rs is under construction, plotting is not implemented yet")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cli_is_consistent() {
        use clap::CommandFactory;
        Cli::command().debug_assert();
    }
}
