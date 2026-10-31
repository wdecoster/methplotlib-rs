//! bedMethyl and bedgraph: aggregated modification frequencies or values per position or interval

use std::collections::BTreeMap;

use anyhow::{Context, Result, bail};
use log::warn;

use crate::region::Window;
use crate::sniff::FileType;
use crate::text::lines_in_window;

/// The modification frequency at the start of every record, for one modification
#[derive(Debug, Clone)]
pub struct Frequencies {
    pub name: String,
    pub start: Vec<u64>,
    /// smoothed, and NaN where unknown
    pub frequency: Vec<f64>,
}

/// The value of every interval of a bedgraph
#[derive(Debug, Clone)]
pub struct Bedgraph {
    pub name: String,
    pub start: Vec<u64>,
    pub end: Vec<u64>,
    pub value: Vec<f64>,
}

/// Records overlapping the half-open window, as pyranges selects them
fn overlaps(start: u64, end: u64, window: &Window) -> bool {
    start < window.end && end > window.begin
}

fn field<'a>(fields: &[&'a str], i: usize, path: &str, line: &str) -> Result<&'a str> {
    fields.get(i).copied().with_context(|| {
        format!(
            "{path}: expected at least {} columns in line\n{line}",
            i + 1
        )
    })
}

fn number<T: std::str::FromStr>(s: &str, path: &str, line: &str) -> Result<T> {
    s.trim()
        .parse()
        .ok()
        .with_context(|| format!("{path}: '{s}' is not a number in line\n{line}"))
}

/// The modification code as in the bedMethyl, from the --mods as in the MM tag: C+m is m.
/// Without a '+', the name is used as is, e.g. 5mC for modbam2bed.
fn bedmethyl_mod_name(m: &str) -> &str {
    match m.split_once('+') {
        Some((_, code)) => code,
        None => m,
    }
}

pub fn parse_bedmethyl(
    path: &str,
    name: &str,
    window: &Window,
    smooth: usize,
    file_type: FileType,
    mods: Option<&str>,
) -> Result<Vec<Frequencies>> {
    // per modification: start and frequency of every record in the window
    let mut per_mod: BTreeMap<String, Vec<(u64, f64)>> = BTreeMap::new();
    for line in lines_in_window(path, window, "tabix -p bed")? {
        // modkit separates the columns after BED9 by spaces, unless run with --only-tabs
        let fields: Vec<&str> = line.split_ascii_whitespace().collect();
        let start: u64 = number(field(&fields, 1, path, &line)?, path, &line)?;
        let end: u64 = number(field(&fields, 2, path, &line)?, path, &line)?;
        if field(&fields, 0, path, &line)? != window.chromosome || !overlaps(start, end, window) {
            continue;
        }
        let modification = field(&fields, 3, path, &line)?;
        let frequency = match file_type {
            FileType::BedMethylModkit => {
                number::<f64>(field(&fields, 10, path, &line)?, path, &line)? / 100.0
            }
            FileType::BedMethylModbam2bed => {
                let canonical: f64 = number(field(&fields, 11, path, &line)?, path, &line)?;
                let modified: f64 = number(field(&fields, 12, path, &line)?, path, &line)?;
                modified / (canonical + modified)
            }
            _ => unreachable!("not a bedMethyl file type: {file_type:?}"),
        };
        per_mod
            .entry(modification.to_string())
            .or_default()
            .push((start, frequency));
    }
    if per_mod.is_empty() {
        no_records_warning(path, window);
        return Ok(vec![]);
    }
    if let Some(mods) = mods {
        let wanted: Vec<&str> = mods.split(',').map(bedmethyl_mod_name).collect();
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
        .map(|(modification, mut records)| {
            records.sort_by_key(|r| r.0);
            let (start, frequency): (Vec<u64>, Vec<f64>) = records.into_iter().unzip();
            Frequencies {
                name: format!("{name}_{modification}"),
                start,
                frequency: centered_rolling_mean(&frequency, smooth),
            }
        })
        .collect())
}

pub fn parse_bedgraph(path: &str, name: &str, window: &Window) -> Result<Option<Bedgraph>> {
    let mut records: Vec<(u64, u64, f64)> = vec![];
    for line in lines_in_window(path, window, "tabix -p bed")? {
        // a bedgraph from `modkit pileup --bedgraph` has a fifth column with the coverage,
        // which is ignored here
        let fields: Vec<&str> = line.split('\t').collect();
        let start: u64 = number(field(&fields, 1, path, &line)?, path, &line)?;
        let end: u64 = number(field(&fields, 2, path, &line)?, path, &line)?;
        if field(&fields, 0, path, &line)? != window.chromosome || !overlaps(start, end, window) {
            continue;
        }
        records.push((
            start,
            end,
            number(field(&fields, 3, path, &line)?, path, &line)?,
        ));
    }
    if records.is_empty() {
        no_records_warning(path, window);
        return Ok(None);
    }
    records.sort_by_key(|r| r.0);
    Ok(Some(Bedgraph {
        name: name.to_string(),
        start: records.iter().map(|r| r.0).collect(),
        end: records.iter().map(|r| r.1).collect(),
        value: records.iter().map(|r| r.2).collect(),
    }))
}

/// Not fatal, as other files may still have data in the window
fn no_records_warning(path: &str, window: &Window) {
    warn!(
        "No records for {path} in {}:{}-{}! This dataset is skipped.",
        window.chromosome, window.begin, window.end
    );
}

/// The rolling mean of pandas' `rolling(window=size, center=True, min_periods=1).mean()`:
/// for an even size the window has one more record before than after the center.
/// Averaging over fewer records at the edges avoids losing those records, and NaN is skipped.
pub fn centered_rolling_mean(values: &[f64], size: usize) -> Vec<f64> {
    let size = size.max(1);
    let after = (size - 1) / 2;
    let before = size - 1 - after;
    (0..values.len())
        .map(|i| {
            let lo = i.saturating_sub(before);
            let hi = (i + after).min(values.len() - 1);
            let (sum, n) = values[lo..=hi]
                .iter()
                .filter(|v| !v.is_nan())
                .fold((0.0, 0), |(s, n), v| (s + v, n + 1));
            if n == 0 { f64::NAN } else { sum / n as f64 }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rolling_mean_matches_pandas() {
        // pandas.Series([1, 2, 3, 4, 5, 6]).rolling(4, center=True, min_periods=1).mean()
        let v = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        assert_eq!(
            centered_rolling_mean(&v, 4),
            vec![1.5, 2.0, 2.5, 3.5, 4.5, 5.0]
        );
        // .rolling(5, center=True, min_periods=1)
        assert_eq!(
            centered_rolling_mean(&v, 5),
            vec![2.0, 2.5, 3.0, 4.0, 4.5, 5.0]
        );
        assert_eq!(centered_rolling_mean(&v, 1), v.to_vec());
        let with_nan = centered_rolling_mean(&[f64::NAN, 1.0, f64::NAN], 1);
        assert!(with_nan[0].is_nan() && with_nan[1] == 1.0 && with_nan[2].is_nan());
    }

    #[test]
    fn mod_names() {
        assert_eq!(bedmethyl_mod_name("C+m"), "m");
        assert_eq!(bedmethyl_mod_name("5mC"), "5mC");
    }
}
