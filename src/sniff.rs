//! Recognize the type of an input file, as the Python version's file_sniffer

use std::fs::File;
use std::io::{BufRead, Read};

use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;

use crate::text::{is_data_line, open_text};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileType {
    Cram,
    Bam,
    /// from `modkit pileup`: 18 columns, separated by tabs and spaces unless --only-tabs was used
    BedMethylModkit,
    /// from `modbam2bed --extended`: 14 columns
    BedMethylModbam2bed,
    /// 4 columns, or 5 from `modkit pileup --bedgraph`
    Bedgraph,
}

pub fn file_type(path: &str) -> Result<FileType> {
    if !std::path::Path::new(path).is_file() {
        bail!("file {path} does not exist, please check the path");
    }
    let mut magic = [0u8; 4];
    let n = File::open(path)?.read(&mut magic)?;
    if n == 4 && &magic == b"CRAM" {
        return Ok(FileType::Cram);
    }
    if n >= 2 && magic[..2] == [0x1f, 0x8b] {
        let mut decompressed = [0u8; 4];
        let mut decoder = MultiGzDecoder::new(File::open(path)?);
        if decoder.read_exact(&mut decompressed).is_ok() && &decompressed == b"BAM\x01" {
            return Ok(FileType::Bam);
        }
    }
    let line = first_data_line(path)?;
    if line.contains("log_lik_methylated")
        || line.contains("num_motifs_in_group")
        || line.contains("GMM_anova_pvalue")
        || line.contains("GMM_logit_pvalue")
    {
        bail!(
            "{path} is in a nanopolish or nanocompore format, which methplotlib-rs does not support. \
             Please use the Python version of methplotlib for these."
        );
    }
    let fields: Vec<&str> = line.split_ascii_whitespace().collect();
    if is_bedmethyl(&fields) {
        match fields.len() {
            18 => return Ok(FileType::BedMethylModkit),
            14 => return Ok(FileType::BedMethylModbam2bed),
            _ => {}
        }
    }
    if is_bedgraph(&line.split('\t').collect::<Vec<_>>()) {
        return Ok(FileType::Bedgraph);
    }
    bail!(unrecognized_file_error(path, &line))
}

fn first_data_line(path: &str) -> Result<String> {
    for line in open_text(path)?.lines() {
        let line = line.with_context(|| format!("failed to read {path}"))?;
        if is_data_line(&line) {
            return Ok(line.trim_end_matches('\r').to_string());
        }
    }
    bail!("file {path} does not contain any data")
}

fn is_digits(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit())
}

/// BED9 with a modification code in column 4, followed by counts
fn is_bedmethyl(fields: &[&str]) -> bool {
    fields.len() > 9
        && is_digits(fields[1])
        && is_digits(fields[2])
        && ["+", "-", "."].contains(&fields[5])
}

/// chromosome, start, end and a numeric value, optionally followed by the coverage
fn is_bedgraph(fields: &[&str]) -> bool {
    [4, 5].contains(&fields.len())
        && is_digits(fields[1])
        && is_digits(fields[2])
        && fields[3].trim().parse::<f64>().is_ok()
}

fn unrecognized_file_error(path: &str, line: &str) -> String {
    let short: String = line.chars().take(200).collect();
    format!(
        "input file {path} was not recognized!\n\n\
         Its first data line has {} tab separated fields ({} when splitting on any whitespace):\n\
         {short}\n\n\
         methplotlib accepts:\n \
         - bedgraph: 4 columns (chromosome, start, end, value), or the 5 columns of `modkit pileup --bedgraph`\n \
         - bedMethyl of `modkit pileup`: 18 columns, with or without --only-tabs\n \
         - bedMethyl of `modbam2bed --extended`: 14 columns\n \
         - bam or cram with MM/ML tags\n\n\
         Frequent causes are a header line, a file that is comma or space separated throughout, \
         or extra columns added by other tools.",
        line.split('\t').count(),
        line.split_ascii_whitespace().count()
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recognizes_the_test_data() {
        for (path, expected) in [
            ("test-data/reads.bam", FileType::Bam),
            ("test-data/reads.cram", FileType::Cram),
            ("test-data/modkit.bed.gz", FileType::BedMethylModkit),
            ("test-data/modkit_only_tabs.bed", FileType::BedMethylModkit),
            ("test-data/modbam2bed.bed.gz", FileType::BedMethylModbam2bed),
            ("test-data/values.bedgraph.gz", FileType::Bedgraph),
            ("test-data/modkit.bedgraph", FileType::Bedgraph),
        ] {
            assert_eq!(file_type(path).unwrap(), expected, "{path}");
        }
    }

    #[test]
    fn rejects_other_files() {
        assert!(file_type("test-data/ref.fa").is_err());
        assert!(file_type("test-data/does-not-exist.bed").is_err());
    }
}
