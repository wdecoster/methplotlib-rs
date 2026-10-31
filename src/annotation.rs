//! Gene annotation from gtf/gff, and regions from bed

use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use std::path::{Path, PathBuf};
use std::time::Duration;

use anyhow::{Context, Result, bail};
use log::{info, warn};

use crate::region::Window;
use crate::text::{is_data_line, open_text};

/// The MANE Select transcripts of GRCh38, available as the built-in annotation "mane"
pub const MANE_VERSION: &str = "1.5";

fn mane_file() -> String {
    format!("MANE.GRCh38.v{MANE_VERSION}.ensembl_genomic.gtf.gz")
}

fn mane_url() -> String {
    format!(
        "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_{MANE_VERSION}/{}",
        mane_file()
    )
}

/// plotly's DEFAULT_PLOTLY_COLORS, which the Python version cycles through per gene
pub const GENE_COLORS: [&str; 10] = [
    "rgb(31, 119, 180)",
    "rgb(255, 127, 14)",
    "rgb(44, 160, 44)",
    "rgb(214, 39, 40)",
    "rgb(148, 103, 189)",
    "rgb(140, 86, 75)",
    "rgb(227, 119, 194)",
    "rgb(127, 127, 127)",
    "rgb(188, 189, 34)",
    "rgb(23, 190, 207)",
];

/// A transcript, or a gene with --simplify, with its exons
#[derive(Debug, Clone, PartialEq)]
pub struct Transcript {
    pub transcript: String,
    pub gene: String,
    /// sorted by begin, in the coordinates of the annotation file
    pub exons: Vec<(u64, u64)>,
    pub strand: String,
    pub begin: u64,
    pub end: u64,
    pub color: String,
}

/// A region of a bed file
#[derive(Debug, Clone, PartialEq)]
pub struct BedRegion {
    pub begin: u64,
    pub end: u64,
    pub name: String,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum AnnotationType {
    Gtf,
    Gff,
}

#[derive(Debug, Clone)]
struct Feature {
    is_gene: bool,
    begin: u64,
    end: u64,
    strand: String,
    gene: Option<String>,
    transcript: Option<String>,
}

/// The directory to keep downloaded annotation in, following the XDG convention
fn cache_dir() -> Result<PathBuf> {
    let base = match std::env::var_os("XDG_CACHE_HOME").filter(|v| !v.is_empty()) {
        Some(dir) => PathBuf::from(dir),
        None => PathBuf::from(
            std::env::var_os("HOME").context("neither XDG_CACHE_HOME nor HOME is set")?,
        )
        .join(".cache"),
    };
    Ok(base.join("methplotlib"))
}

/// The path of the annotation to use for --gtf.
///
/// "mane" is the built-in MANE Select annotation of GRCh38, downloaded once and kept in the
/// cache directory, unless a file named "mane" exists.
pub fn resolve_annotation(gtf: &str) -> Result<String> {
    if !gtf.eq_ignore_ascii_case("mane") || Path::new(gtf).exists() {
        return Ok(gtf.to_string());
    }
    let path = cache_dir()?.join(mane_file());
    if !path.is_file() {
        download(&mane_url(), &path).with_context(|| {
            format!("could not download the MANE annotation from {}", mane_url())
        })?;
    }
    Ok(path.to_string_lossy().into_owned())
}

/// Download to a temporary name first, so that an interrupted download is not used later
fn download(url: &str, path: &Path) -> Result<()> {
    eprintln!(
        "Downloading the MANE v{MANE_VERSION} annotation to {}",
        path.display()
    );
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir)
            .with_context(|| format!("failed to create {}", dir.display()))?;
    }
    let partial = path.with_extension("gz.part");
    // the connect timeout is shared by the addresses of the server, of which an unreachable
    // IPv6 address may be tried before IPv4
    let agent: ureq::Agent = ureq::Agent::config_builder()
        .timeout_connect(Some(Duration::from_secs(20)))
        .timeout_recv_body(Some(Duration::from_secs(600)))
        .build()
        .into();
    let result = (|| -> Result<()> {
        let response = agent.get(url).call()?;
        let mut out = std::fs::File::create(&partial)
            .with_context(|| format!("failed to create {}", partial.display()))?;
        std::io::copy(&mut response.into_body().into_reader(), &mut out)?;
        Ok(())
    })();
    if result.is_err() {
        let _ = std::fs::remove_file(&partial);
    }
    result?;
    std::fs::rename(&partial, path)
        .with_context(|| format!("failed to move the download to {}", path.display()))
}

fn annotation_type(path: &str) -> Result<AnnotationType> {
    if path.ends_with(".gtf") || path.ends_with(".gtf.gz") {
        Ok(AnnotationType::Gtf)
    } else if path.ends_with(".gff")
        || path.ends_with(".gff.gz")
        || path.ends_with(".gff3")
        || path.ends_with(".gff3.gz")
    {
        Ok(AnnotationType::Gff)
    } else {
        bail!(
            "unrecognized extension of the annotation file {path}.\n\
             Supported are gtf, gtf.gz, gff, gff.gz, gff3 and gff3.gz"
        )
    }
}

/// The gene name and transcript id of the attributes, or the locus tag for both if there is no gene name.
///
/// As the Python version: only the part after the key up to the next delimiter is the value,
/// with quotes and a trailing ';' removed.
fn parse_attributes(attributes: &str, kind: AnnotationType) -> (Option<String>, Option<String>) {
    let (attribute_delimiter, kv_delimiter) = match kind {
        AnnotationType::Gtf => ("; ", ' '),
        AnnotationType::Gff => (";", '='),
    };
    let mut info: HashMap<&str, String> = HashMap::new();
    for item in attributes.split(attribute_delimiter) {
        if !(item.starts_with("gene_name")
            || item.starts_with("transcript_id")
            || item.starts_with("locus_tag"))
        {
            continue;
        }
        let mut parts = item.split(kv_delimiter);
        let key = parts.next().unwrap_or_default();
        let value = parts.next().unwrap_or_default();
        info.insert(
            key,
            value.replace('"', "").trim_end_matches(';').to_string(),
        );
    }
    if info.contains_key("gene_name") {
        (
            info.get("gene_name").cloned(),
            info.get("transcript_id").cloned(),
        )
    } else {
        (
            info.get("locus_tag").cloned(),
            info.get("locus_tag").cloned(),
        )
    }
}

/// The transcripts (or genes with simplify) with an exon or gene record overlapping the window.
///
/// As the Python version: exon records are used, and the gene record of a gene without exon
/// records, which would otherwise be drawn as one exon over the introns. With simplify, a gene is
/// a meta-transcript of the merged exons of its transcripts.
/// The begin and end are as in the file, and compared inclusively with the window.
pub fn parse_annotation(path: &str, window: &Window, simplify: bool) -> Result<Vec<Transcript>> {
    let kind = annotation_type(path)?;
    info!("Parsing {kind:?} file {path}");
    let mut chromosomes_seen = HashSet::new();
    let mut features = vec![];
    for line in open_text(path)?.lines() {
        let line = line.with_context(|| format!("failed to read {path}"))?;
        if !is_data_line(&line) {
            continue;
        }
        let fields: Vec<&str> = line.splitn(4, '\t').collect();
        if fields.len() > 2
            && fields[0] == window.chromosome
            && (fields[2] == "exon" || fields[2] == "gene")
        {
            let all: Vec<&str> = line.split('\t').collect();
            if all.len() != 9 {
                bail!("expected 9 columns in {path} for line\n{line}");
            }
            let (gene, transcript) = parse_attributes(all[8].trim_end(), kind);
            features.push(Feature {
                is_gene: fields[2] == "gene",
                begin: all[3]
                    .parse()
                    .with_context(|| format!("invalid begin in {path}:\n{line}"))?,
                end: all[4]
                    .parse()
                    .with_context(|| format!("invalid end in {path}:\n{line}"))?,
                strand: all[6].to_string(),
                gene,
                transcript,
            });
        } else {
            chromosomes_seen.insert(fields[0].to_string());
        }
    }
    if features.is_empty() && !chromosomes_seen.is_empty() {
        let mut examples: Vec<&String> = chromosomes_seen.iter().collect();
        examples.sort();
        warn!(
            "no exons or genes of chromosome '{}' in {path}, which contains {} chromosomes, such as {}. \
             Do the chromosome names in --window match those in your annotation (e.g. '15' versus 'chr15')?",
            window.chromosome,
            chromosomes_seen.len(),
            examples
                .iter()
                .take(5)
                .map(|s| s.as_str())
                .collect::<Vec<_>>()
                .join(", ")
        );
    }
    let genes_with_exons: HashSet<String> = features
        .iter()
        .filter(|f| !f.is_gene)
        .filter_map(|f| f.gene.clone())
        .collect();
    features.retain(|f| {
        !(f.is_gene
            && f.gene
                .as_ref()
                .is_some_and(|g| genes_with_exons.contains(g)))
    });
    let overlapping = |f: &Feature| f.begin <= window.end && f.end >= window.begin;
    let mut result = vec![];
    if simplify {
        for gene in unique_in_order(features.iter().filter(|f| overlapping(f)).map(|f| &f.gene)) {
            let records: Vec<&Feature> = features
                .iter()
                .filter(|f| f.gene.as_ref() == Some(gene))
                .collect();
            let mut meta = make_transcript(gene, gene, &records);
            meta.exons = merge_intervals(&meta.exons);
            result.push(meta);
        }
        eprintln!("Found {} gene(s) in the region.", result.len());
    } else {
        for transcript in unique_in_order(
            features
                .iter()
                .filter(|f| overlapping(f))
                .map(|f| &f.transcript),
        ) {
            let records: Vec<&Feature> = features
                .iter()
                .filter(|f| f.transcript.as_ref() == Some(transcript))
                .collect();
            let gene = records[0].gene.clone().unwrap_or_default();
            result.push(make_transcript(transcript, &gene, &records));
        }
        eprintln!("Found {} transcript(s) in the region.", result.len());
    }
    assign_colors_to_genes(&mut result);
    Ok(result)
}

/// The distinct values in order of appearance, without the records missing it
fn unique_in_order<'a>(values: impl Iterator<Item = &'a Option<String>>) -> Vec<&'a String> {
    let mut seen = HashSet::new();
    values
        .flatten()
        .filter(|v| seen.insert(v.as_str()))
        .collect()
}

fn make_transcript(transcript: &str, gene: &str, records: &[&Feature]) -> Transcript {
    let mut exons: Vec<(u64, u64)> = records.iter().map(|f| (f.begin, f.end)).collect();
    exons.sort_by_key(|e| e.0);
    Transcript {
        transcript: transcript.to_string(),
        gene: gene.to_string(),
        strand: records[0].strand.clone(),
        begin: exons
            .iter()
            .flat_map(|e| [e.0, e.1])
            .min()
            .unwrap_or_default(),
        end: exons
            .iter()
            .flat_map(|e| [e.0, e.1])
            .max()
            .unwrap_or_default(),
        exons,
        color: String::new(),
    }
}

/// Merge overlapping and adjacent intervals, of which the begin and end are inclusive
fn merge_intervals(intervals: &[(u64, u64)]) -> Vec<(u64, u64)> {
    let mut sorted = intervals.to_vec();
    sorted.sort();
    let mut merged: Vec<(u64, u64)> = vec![];
    for (begin, end) in sorted {
        match merged.last_mut() {
            Some(last) if begin <= last.1 + 1 => last.1 = last.1.max(end),
            _ => merged.push((begin, end)),
        }
    }
    merged
}

fn assign_colors_to_genes(transcripts: &mut [Transcript]) {
    let mut colors: HashMap<String, &str> = HashMap::new();
    for t in transcripts.iter_mut() {
        let n = colors.len();
        let color = *colors
            .entry(t.gene.clone())
            .or_insert(GENE_COLORS[n % GENE_COLORS.len()]);
        t.color = color.to_string();
    }
}

/// The regions of a bed file overlapping the window, named "noname" if the file has no names
pub fn parse_bed(path: &str, window: &Window) -> Result<Vec<BedRegion>> {
    info!("Parsing bed file {path}");
    let mut regions = vec![];
    for line in open_text(path)?.lines() {
        let line = line.with_context(|| format!("failed to read {path}"))?;
        if !is_data_line(&line) {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 || fields[0] != window.chromosome {
            continue;
        }
        let begin: u64 = fields[1]
            .trim()
            .parse()
            .with_context(|| format!("invalid begin in {path}:\n{line}"))?;
        let end: u64 = fields[2]
            .trim()
            .parse()
            .with_context(|| format!("invalid end in {path}:\n{line}"))?;
        if begin < window.end && end > window.begin {
            regions.push(BedRegion {
                begin,
                end,
                name: fields.get(3).map_or("noname", |n| n.trim()).to_string(),
            });
        }
    }
    regions.sort_by_key(|r| (r.begin, r.end));
    Ok(regions)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::region::Region;

    #[test]
    fn intervals_are_merged() {
        assert_eq!(
            merge_intervals(&[(5, 8), (1, 3), (4, 4), (10, 12), (11, 20)]),
            vec![(1, 8), (10, 20)]
        );
    }

    #[test]
    fn simplify_merges_exons_and_ignores_gene_records() {
        let dir = std::env::temp_dir().join(format!("methplotlib-gtf-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("test.gtf");
        let gene = r#"gene_id "G1"; gene_name "GENE";"#;
        let records = [
            ("gene", 100, 1000, gene.to_string()),
            ("exon", 100, 200, format!(r#"{gene} transcript_id "T1";"#)),
            ("exon", 500, 600, format!(r#"{gene} transcript_id "T1";"#)),
            ("exon", 150, 250, format!(r#"{gene} transcript_id "T2";"#)),
            ("exon", 900, 1000, format!(r#"{gene} transcript_id "T2";"#)),
            (
                "gene",
                2000,
                3000,
                r#"gene_id "G2"; gene_name "NOEXONS";"#.to_string(),
            ),
        ];
        let gtf: String = records
            .iter()
            .map(|(feature, begin, end, attributes)| {
                format!("chr1\ttest\t{feature}\t{begin}\t{end}\t.\t+\t.\t{attributes}\n")
            })
            .collect();
        std::fs::write(&path, gtf).unwrap();
        let path = path.to_str().unwrap();
        let window = Region::parse("chr1:1-5000")
            .unwrap()
            .windows(None)
            .unwrap()
            .remove(0);
        let genes = parse_annotation(path, &window, true).unwrap();
        let genes: Vec<(&str, &[(u64, u64)])> = genes
            .iter()
            .map(|g| (g.gene.as_str(), g.exons.as_slice()))
            .collect();
        assert_eq!(
            genes,
            vec![
                ("GENE", &[(100, 250), (500, 600), (900, 1000)][..]),
                ("NOEXONS", &[(2000, 3000)][..]),
            ]
        );
        let transcripts = parse_annotation(path, &window, false).unwrap();
        let transcripts: Vec<(&str, &[(u64, u64)])> = transcripts
            .iter()
            .map(|t| (t.transcript.as_str(), t.exons.as_slice()))
            .collect();
        assert_eq!(
            transcripts,
            vec![
                ("T1", &[(100, 200), (500, 600)][..]),
                ("T2", &[(150, 250), (900, 1000)][..]),
            ]
        );
        std::fs::remove_dir_all(dir).unwrap();
    }
}
