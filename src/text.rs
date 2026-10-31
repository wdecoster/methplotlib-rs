//! Reading the lines of a tab separated file (bedMethyl, bedgraph) in a window,
//! with tabix when the file has an index, or by reading the entire file otherwise.

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use log::{info, warn};
use noodles::core::{Position, Region as TabixRegion};
use noodles::csi::BinningIndex;
use noodles::tabix;

use crate::region::Window;

/// Open a file for reading lines, transparently decompressing gzip (and so bgzip)
pub fn open_text(path: &str) -> Result<Box<dyn BufRead>> {
    let mut file = File::open(path).with_context(|| format!("failed to open {path}"))?;
    let mut magic = [0u8; 2];
    let n = file.read(&mut magic)?;
    let file = File::open(path)?;
    if n == 2 && magic == [0x1f, 0x8b] {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Lines that carry no records: empty lines, comments and the track/browser lines of a UCSC file
pub fn is_data_line(line: &str) -> bool {
    !(line.trim().is_empty()
        || line.starts_with('#')
        || line.starts_with("track")
        || line.starts_with("browser"))
}

/// The lines of the file on the chromosome of the window.
///
/// With a tabix index only the lines around the window are read, otherwise the entire file.
/// The caller still selects the records overlapping the window: tabix returns those touching it,
/// and without an index all lines of the chromosome are returned.
pub fn lines_in_window(path: &str, window: &Window, index_suggestion: &str) -> Result<Vec<String>> {
    if Path::new(&format!("{path}.tbi")).is_file() {
        info!("Reading {path} with its tabix index");
        let mut reader = tabix::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .with_context(|| format!("failed to open {path} with its tabix index"))?;
        let known = reader.index().header().is_some_and(|h| {
            h.reference_sequence_names()
                .contains(window.chromosome.as_bytes())
        });
        if !known {
            // tabix itself silently returns nothing for a chromosome that is not in the index
            return Ok(vec![]);
        }
        // tabix regions are 1-based and inclusive: this covers the 0-based window [begin, end)
        let region = TabixRegion::new(
            window.chromosome.as_str(),
            Position::try_from(window.begin as usize + 1)?
                ..=Position::try_from(window.end as usize)?,
        );
        let mut lines = vec![];
        for record in reader.query(&region)? {
            lines.push(record?.as_ref().to_string());
        }
        Ok(lines)
    } else {
        info!("Reading {path} entirely, as it has no tabix index");
        eprintln!("Reading {path} would be faster with bgzip and '{index_suggestion}'.");
        let mut lines = vec![];
        let mut chromosomes_seen = std::collections::BTreeSet::new();
        for line in open_text(path)?.lines() {
            let line = line.with_context(|| format!("failed to read {path}"))?;
            if !is_data_line(&line) {
                continue;
            }
            let chromosome = line.split('\t').next().unwrap_or_default();
            if chromosome == window.chromosome {
                lines.push(line);
            } else if !chromosomes_seen.contains(chromosome) {
                chromosomes_seen.insert(chromosome.to_string());
            }
        }
        if lines.is_empty() && !chromosomes_seen.is_empty() {
            let examples: Vec<&str> = chromosomes_seen
                .iter()
                .take(5)
                .map(String::as_str)
                .collect();
            warn!(
                "chromosome '{}' not found in {path}, which contains {} chromosomes, such as {}. \
                 Do the chromosome names in --window match those in your file (e.g. '1' versus 'chr1')?",
                window.chromosome,
                chromosomes_seen.len(),
                examples.join(", ")
            );
        }
        Ok(lines)
    }
}
