//! Per-read modification calls from the MM/ML tags of bam or cram

use std::collections::BTreeMap;

use anyhow::{Context, Result, bail};
use noodles::core::{Position, Region as HtsRegion};
use noodles::sam::alignment::Record;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::value::Array;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use noodles::{bam, cram, fasta};

use crate::region::Window;

/// A primary alignment overlapping the window
#[derive(Debug, Clone)]
pub struct Read {
    pub name: String,
    /// 0-based start and exclusive end on the reference
    pub start: u64,
    pub end: u64,
    pub reverse: bool,
}

/// One modification call of a read
#[derive(Debug, Clone, PartialEq)]
pub struct Call {
    /// index of the read in [`ReadCalls::reads`]
    pub read: usize,
    /// the reference position, None for a base without one (softclip, insertion)
    pub pos: Option<u64>,
    /// the probability of the modification, on the 0-255 scale of the ML tag
    pub quality: u8,
}

/// The calls of one modification, with all reads in the window
#[derive(Debug, Clone)]
pub struct ReadCalls {
    pub name: String,
    pub reads: Vec<Read>,
    /// sorted by read name, then position
    pub calls: Vec<Call>,
}

/// The reference position (if any) and probability of a call
type PositionCall = (Option<u64>, u8);

const LEGACY_MODIFICATIONS: Tag = Tag::new(b'M', b'm');
const LEGACY_PROBABILITIES: Tag = Tag::new(b'M', b'l');

/// The calls of every modification in a bam or cram, for the primary alignments in the window.
/// A cram is decoded with the fasta, if one is given.
pub fn parse_alignments(
    path: &str,
    name: &str,
    window: &Window,
    mods: Option<&str>,
    fasta: Option<&str>,
    is_cram: bool,
) -> Result<Vec<ReadCalls>> {
    // the query is 1-based and inclusive: this covers the 0-based window [begin, end)
    let region = HtsRegion::new(
        window.chromosome.as_str(),
        Position::try_from(window.begin as usize + 1)?..=Position::try_from(window.end as usize)?,
    );
    let mut reads = vec![];
    // per modification: (read index, position, probability)
    let mut per_mod: BTreeMap<String, Vec<Call>> = BTreeMap::new();
    let mut add = |record: &dyn Record| -> Result<()> {
        let flags = record.flags()?;
        if flags.is_secondary() || flags.is_supplementary() || flags.is_unmapped() {
            return Ok(());
        }
        let read = Read {
            name: record
                .name()
                .map(|n| n.to_string())
                .unwrap_or_else(|| "*".to_string()),
            start: record
                .alignment_start()
                .context("missing alignment start")??
                .get() as u64
                - 1,
            end: record
                .alignment_end()
                .context("missing alignment end")??
                .get() as u64,
            reverse: flags.is_reverse_complemented(),
        };
        let index = reads.len();
        for (modification, calls) in modified_positions(record, &read)? {
            per_mod
                .entry(modification)
                .or_default()
                .extend(calls.into_iter().map(|(pos, quality)| Call {
                    read: index,
                    pos,
                    quality,
                }));
        }
        reads.push(read);
        Ok(())
    };
    let query_error = || {
        format!(
            "failed to query {path} for {}:{}-{}, does it have an index and this chromosome?",
            window.chromosome, window.begin, window.end
        )
    };
    if is_cram {
        let repository = match fasta {
            Some(fasta) => fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(
                fasta::io::indexed_reader::Builder::default()
                    .build_from_path(fasta)
                    .with_context(|| format!("failed to open {fasta}"))?,
            )),
            None => fasta::Repository::default(),
        };
        let mut reader = cram::io::indexed_reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_path(path)
            .with_context(|| format!("failed to open {path} with its index"))?;
        let header = reader.read_header()?;
        for record in reader
            .query(&header, &region)
            .with_context(query_error)?
            .records()
        {
            let record = record.with_context(|| {
                format!("failed to decode {path}, a cram needs its reference, given with --fasta")
            })?;
            add(&record)?;
        }
    } else {
        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .with_context(|| format!("failed to open {path} with its index"))?;
        let header = reader.read_header()?;
        for record in reader
            .query(&header, &region)
            .with_context(query_error)?
            .records()
        {
            add(&record?)?;
        }
    }
    if let Some(mods) = mods {
        let wanted: Vec<&str> = mods.split(',').collect();
        let found: Vec<String> = per_mod.keys().cloned().collect();
        per_mod.retain(|m, _| wanted.contains(&m.as_str()));
        if per_mod.is_empty() {
            bail!(
                "no records left for {path} after selecting --mods!\nModifications in this file: {}",
                found.join(", ")
            );
        }
    }
    Ok(per_mod
        .into_iter()
        .map(|(modification, mut calls)| {
            calls.sort_by(|a, b| {
                reads[a.read]
                    .name
                    .cmp(&reads[b.read].name)
                    .then_with(|| match (a.pos, b.pos) {
                        // positions without a reference coordinate last, as pandas sorts NaN
                        (Some(x), Some(y)) => x.cmp(&y),
                        (Some(_), None) => std::cmp::Ordering::Less,
                        (None, Some(_)) => std::cmp::Ordering::Greater,
                        (None, None) => std::cmp::Ordering::Equal,
                    })
            });
            ReadCalls {
                name: format!("{name}_{modification}"),
                reads: reads.clone(),
                calls,
            }
        })
        .collect())
}

fn get_string(record: &dyn Record, tags: [Tag; 2]) -> Result<Option<String>> {
    let data = record.data();
    for tag in tags {
        if let Some(value) = data.get(&tag) {
            return match value? {
                Value::String(s) => Ok(Some(s.to_string())),
                _ => bail!("the {tag:?} tag is not a string"),
            };
        }
    }
    Ok(None)
}

fn get_probabilities(record: &dyn Record, tags: [Tag; 2]) -> Result<Option<Vec<u8>>> {
    let data = record.data();
    for tag in tags {
        if let Some(value) = data.get(&tag) {
            return match value? {
                Value::Array(Array::UInt8(values)) => {
                    Ok(Some(values.iter().collect::<std::io::Result<Vec<u8>>>()?))
                }
                _ => bail!("the {tag:?} tag is not an array of 8-bit unsigned integers"),
            };
        }
    }
    Ok(None)
}

/// Per modification of the MM tag, the reference position and probability of every call.
///
/// As in the Python version, the modification is named as in the tag without the '?' (e.g. C+m),
/// a missing ML tag means a probability of 255, and calls without a probability in a short ML tag
/// are dropped. Unlike the Python version, a combined code such as C+mh is split in C+m and C+h,
/// of which the probabilities are interleaved in the ML tag.
fn modified_positions(
    record: &dyn Record,
    read: &Read,
) -> Result<Vec<(String, Vec<PositionCall>)>> {
    let Some(mm) = get_string(record, [Tag::BASE_MODIFICATIONS, LEGACY_MODIFICATIONS])? else {
        return Ok(vec![]);
    };
    let probabilities = get_probabilities(
        record,
        [Tag::BASE_MODIFICATION_PROBABILITIES, LEGACY_PROBABILITIES],
    )?;
    // the read as it came from the sequencer, the direction in which the MM tag is specified
    let mut sequence: Vec<u8> = record.sequence().iter().collect();
    let mut refpos = reference_positions(record)?;
    if read.reverse {
        sequence.reverse();
        for base in sequence.iter_mut() {
            *base = complement(*base);
        }
        refpos.reverse();
    }
    let mut result = vec![];
    let mut offset = 0;
    for context in mm.split(';').filter(|c| !c.is_empty()) {
        let mut parts = context.split(',');
        let basemod = parts.next().unwrap_or_default().replace('?', "");
        if basemod.contains('-') {
            bail!(
                "modifications on the negative strand are currently unsupported (read {}).\n\
                 Please open an issue if this would be of interest for you.",
                read.name
            );
        }
        let Some((base, code)) = basemod.split_once('+') else {
            bail!(
                "unexpected modification '{basemod}' in the MM tag of read {}",
                read.name
            );
        };
        // the code below does not work for the ambiguous N, and would search for a literal N
        if base == "N" {
            bail!(
                "modifications of N nucleotides are currently unsupported (read {}).\n\
                 Please open an issue if this would be of interest for you.",
                read.name
            );
        }
        let base = base.as_bytes().first().copied().unwrap_or_default();
        // the positions are encoded as the number of occurrences of that base to skip
        let deltas = parts
            .map(|d| d.trim().parse::<usize>())
            .collect::<std::result::Result<Vec<_>, _>>()
            .with_context(|| format!("malformed MM tag of read {}: {context}", read.name))?;
        if deltas.is_empty() {
            continue;
        }
        // a code of multiple letters is multiple modifications, a ChEBI number is one,
        // and a trailing '.' (unmodified bases are skipped) is kept in the name, as in Python
        let (code, suffix) = match code.strip_suffix('.') {
            Some(c) => (c, "."),
            None => (code, ""),
        };
        let codes: Vec<String> = if code.bytes().all(|b| b.is_ascii_digit()) {
            vec![code.to_string()]
        } else {
            code.chars().map(String::from).collect()
        };
        let base_index: Vec<usize> = sequence
            .iter()
            .enumerate()
            .filter(|(_, b)| b.eq_ignore_ascii_case(&base))
            .map(|(i, _)| i)
            .collect();
        // the read positions with a modification call
        let mut location = 0;
        let mut positions = Vec::with_capacity(deltas.len());
        for (i, delta) in deltas.iter().enumerate() {
            location += delta + usize::from(i > 0);
            let position = base_index.get(location).with_context(|| {
                format!(
                    "the MM tag of read {} refers to more {} bases than there are in the read",
                    read.name, base as char
                )
            })?;
            positions.push(*position);
        }
        // the probabilities are in one array for all contexts, interleaved per modification code
        let n_codes = codes.len();
        let available = match &probabilities {
            Some(p) => (p.len().saturating_sub(offset) / n_codes).min(deltas.len()),
            None => deltas.len(),
        };
        for (k, code) in codes.iter().enumerate() {
            let calls = positions[..available]
                .iter()
                .enumerate()
                .map(|(i, &p)| {
                    let quality = match &probabilities {
                        Some(probabilities) => probabilities[offset + i * n_codes + k],
                        None => 255,
                    };
                    (refpos[p], quality)
                })
                .collect();
            result.push((format!("{}+{code}{suffix}", base as char), calls));
        }
        offset += deltas.len() * n_codes;
    }
    Ok(result)
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        other => other,
    }
}

/// The reference position of every base of the read in alignment direction,
/// None where there is none (insertions, softclips),
/// as pysam's get_reference_positions(full_length=True)
fn reference_positions(record: &dyn Record) -> Result<Vec<Option<u64>>> {
    let mut pos = record
        .alignment_start()
        .context("missing alignment start")??
        .get() as u64
        - 1;
    let mut positions = Vec::with_capacity(record.sequence().len());
    for op in record.cigar().iter() {
        let op = op?;
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                positions.extend((pos..pos + len as u64).map(Some));
                pos += len as u64;
            }
            Kind::Insertion | Kind::SoftClip => positions.extend(std::iter::repeat_n(None, len)),
            Kind::Deletion | Kind::Skip => pos += len as u64,
            Kind::HardClip | Kind::Pad => {}
        }
    }
    Ok(positions)
}

/// The row of every read in the per-read plot, as assign_y_height_per_read in the Python version:
/// reads sorted by start and then longest first, each in the first row where it fits after the
/// read last placed in that row. Reads not fitting in max_rows get None and are hidden.
pub fn assign_rows(reads: &[Read], max_rows: usize) -> Vec<Option<usize>> {
    let mut order: Vec<usize> = (0..reads.len()).collect();
    order.sort_by(|&a, &b| {
        reads[a]
            .start
            .cmp(&reads[b].start)
            .then(reads[b].end.cmp(&reads[a].end))
    });
    let mut rows: Vec<u64> = vec![];
    let mut assigned = vec![None; reads.len()];
    for i in order {
        if let Some(row) = rows.iter().position(|&last| reads[i].start > last) {
            rows[row] = reads[i].end;
            assigned[i] = Some(row + 1);
        } else if rows.len() < max_rows {
            rows.push(reads[i].end);
            assigned[i] = Some(rows.len());
        }
    }
    assigned
}
