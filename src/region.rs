use anyhow::{Result, bail};

/// The window to plot, as given with --window
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Region {
    /// An interval, as chromosome:begin-end
    Interval {
        chromosome: String,
        begin: u64,
        end: u64,
    },
    /// An entire chromosome or contig, of which the length comes from the fasta index
    Contig(String),
}

impl Region {
    /// Parse a window such as chr7:5,525,542-5,543,028 or chr7.
    /// Thousands separators are allowed, as when copied from a genome browser.
    pub fn parse(window: &str) -> Result<Self> {
        let window = window.replace(',', "");
        let Some((chromosome, interval)) = window.rsplit_once(':') else {
            if window.is_empty() {
                bail!("--window is empty");
            }
            return Ok(Region::Contig(window));
        };
        let malformed = || {
            anyhow::anyhow!(
                "--window '{window}' is malformed, \
                 examples of accepted formats are 'chr5:150200605-150423790' or 'chr5'"
            )
        };
        let (begin, end) = interval.split_once('-').ok_or_else(malformed)?;
        let begin: u64 = begin.trim().parse().map_err(|_| malformed())?;
        let end: u64 = end.trim().parse().map_err(|_| malformed())?;
        if chromosome.is_empty() {
            return Err(malformed());
        }
        if begin >= end {
            bail!("--window '{window}': the begin of the interval has to be smaller than its end");
        }
        Ok(Region::Interval {
            chromosome: chromosome.to_string(),
            begin,
            end,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interval_with_thousands_separators() {
        assert_eq!(
            Region::parse("chr7:5,525,542-5,543,028").unwrap(),
            Region::Interval {
                chromosome: "chr7".to_string(),
                begin: 5_525_542,
                end: 5_543_028
            }
        );
    }

    #[test]
    fn contig_name_with_colon() {
        // HLA contigs in GRCh38 contain a colon, so only the last one separates the interval
        assert_eq!(
            Region::parse("HLA-A*01:01:01:01:1-100").unwrap(),
            Region::Interval {
                chromosome: "HLA-A*01:01:01:01".to_string(),
                begin: 1,
                end: 100
            }
        );
    }

    #[test]
    fn entire_contig() {
        assert_eq!(
            Region::parse("chr21").unwrap(),
            Region::Contig("chr21".to_string())
        );
    }

    #[test]
    fn begin_after_end_is_an_error() {
        assert!(Region::parse("chr1:200-100").is_err());
        assert!(Region::parse("chr1:100-100").is_err());
    }

    #[test]
    fn malformed_windows_are_errors() {
        for window in ["chr1:100", "chr1:a-b", ":1-100", "", "chr1:-100"] {
            assert!(Region::parse(window).is_err(), "{window} should fail");
        }
    }
}
