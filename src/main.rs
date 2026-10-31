use anyhow::{Result, bail};
use clap::Parser;
use log::info;

use methplotlib_rs::annotation::{parse_annotation, parse_bed, resolve_annotation};
use methplotlib_rs::frequency::{parse_bedgraph, parse_bedmethyl};
use methplotlib_rs::plot::{Dataset, PlotOptions, make_scene, write_scene};
use methplotlib_rs::reads::parse_alignments;
use methplotlib_rs::region::{Region, Window};
use methplotlib_rs::sniff::{FileType, file_type};

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

    /// Add annotation from a gtf or gff file,
    /// or 'mane' for the MANE Select transcripts of GRCh38 (downloaded once)
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
    /// (both strands are separate records in bedMethyl)
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

    /// Output file, its extension (png, svg or pdf) sets the format; svg gets large with many reads.
    /// Use {region} as a shorthand for {chr}_{start}_{end}.
    /// Default: methylation_browser_{region}.png
    #[arg(short, long)]
    outfile: Option<String>,
}

fn read_dataset(path: &str, name: &str, window: &Window, cli: &Cli) -> Result<Vec<Dataset>> {
    let kind = file_type(path)?;
    info!("File {path} is of type {kind:?}");
    let mods = cli.mods.as_deref();
    Ok(match kind {
        FileType::Bam | FileType::Cram => parse_alignments(
            path,
            name,
            window,
            mods,
            cli.fasta.as_deref(),
            kind == FileType::Cram,
        )?
        .into_iter()
        .map(Dataset::Reads)
        .collect(),
        FileType::BedMethylModkit | FileType::BedMethylModbam2bed => {
            parse_bedmethyl(path, name, window, cli.smooth, kind, mods)?
                .into_iter()
                .map(Dataset::Frequencies)
                .collect()
        }
        FileType::Bedgraph => parse_bedgraph(path, name, window)?
            .into_iter()
            .map(Dataset::Bedgraph)
            .collect(),
    })
}

fn output_path(cli: &Cli, window: &Window) -> Result<String> {
    let path = match &cli.outfile {
        Some(outfile) => outfile.replace("{region}", &window.label),
        None => format!("methylation_browser_{}.png", window.label),
    };
    let lower = path.to_lowercase();
    if ![".svg", ".png", ".pdf"]
        .iter()
        .any(|ext| lower.ends_with(ext))
    {
        bail!("unsupported output format for {path}: use a .svg, .png or .pdf extension");
    }
    Ok(path)
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("warn")).init();
    let cli = Cli::parse();
    if cli.names.len() != cli.methylation.len() {
        bail!(
            "expecting the same number of --names ({}) as --methylation datasets ({})",
            cli.names.len(),
            cli.methylation.len()
        );
    }
    let windows = Region::parse(&cli.window)?.windows(cli.fasta.as_deref())?;
    let gtf = cli.gtf.as_deref().map(resolve_annotation).transpose()?;
    let options = PlotOptions {
        split: cli.split,
        dotsize: f64::from(cli.dotsize),
        minqual: cli.minqual,
    };
    for window in &windows {
        info!(
            "Processing {}:{}-{}",
            window.chromosome, window.begin, window.end
        );
        let outfile = output_path(&cli, window)?;
        let mut datasets = vec![];
        for (path, name) in cli.methylation.iter().zip(&cli.names) {
            datasets.extend(read_dataset(path, name, window, &cli)?);
        }
        if datasets.is_empty() {
            bail!(
                "no data found in {}:{}-{} for any of the input files!\n\
                 Please check your --window and the chromosome names in your input.",
                window.chromosome,
                window.begin,
                window.end
            );
        }
        let transcripts = match &gtf {
            Some(gtf) => parse_annotation(gtf, window, cli.simplify)?,
            None => vec![],
        };
        let regions = match &cli.bed {
            Some(bed) => parse_bed(bed, window)?,
            None => vec![],
        };
        let scene = make_scene(datasets, transcripts, regions, window, &options);
        write_scene(&scene, &outfile)?;
        eprintln!("Wrote {outfile}");
    }
    Ok(())
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
