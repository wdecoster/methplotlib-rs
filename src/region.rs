use anyhow::{Context, Result, bail};

/// Windows larger than this are split in chunks, each plotted separately
const MAX_WINDOW_SIZE: u64 = 1_000_000;

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

/// A concrete window to plot: 0-based, half-open coordinates on a chromosome
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Window {
    pub chromosome: String,
    pub begin: u64,
    pub end: u64,
    /// Used in the default output file name, as {region}
    pub label: String,
}

impl Window {
    fn new(chromosome: &str, begin: u64, end: u64) -> Self {
        Window {
            chromosome: chromosome.to_string(),
            begin,
            end,
            label: format!("{chromosome}_{begin}_{end}"),
        }
    }
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

    /// The windows to plot: the region itself, or chunks of at most 1 Mb when it is larger.
    /// An entire contig needs the fasta, for its length.
    pub fn windows(&self, fasta: Option<&str>) -> Result<Vec<Window>> {
        let (chromosome, begin, end) = match self {
            Region::Interval {
                chromosome,
                begin,
                end,
            } => (chromosome.as_str(), *begin, *end),
            Region::Contig(name) => {
                let Some(fasta) = fasta else {
                    bail!(
                        "a fasta reference (--fasta) is required when --window is an entire chromosome, \
                         contig or transcript"
                    );
                };
                (name.as_str(), 0, contig_length(fasta, name)?)
            }
        };
        let size = end - begin;
        if size <= MAX_WINDOW_SIZE {
            let mut window = Window::new(chromosome, begin, end);
            if let Region::Contig(name) = self {
                window.label = name.clone();
            }
            return Ok(vec![window]);
        }
        let chunks = size.div_ceil(MAX_WINDOW_SIZE);
        let chunk_size = size.div_ceil(chunks);
        // unlike the Python version, the last chunk ends at the end of the region rather than beyond it
        Ok((0..chunks)
            .map(|i| {
                let chunk_begin = begin + i * chunk_size;
                Window::new(chromosome, chunk_begin, (chunk_begin + chunk_size).min(end))
            })
            .collect())
    }
}

/// The length of a contig in the fasta, from its index, which is made if there is none
fn contig_length(fasta: &str, name: &str) -> Result<u64> {
    let fai = format!("{fasta}.fai");
    let index = if std::path::Path::new(&fai).is_file() {
        noodles::fasta::fai::fs::read(&fai).with_context(|| format!("failed to read {fai}"))?
    } else {
        noodles::fasta::fs::index(fasta).with_context(|| format!("failed to index {fasta}"))?
    };
    index
        .as_ref()
        .iter()
        .find(|record| record.name() == name.as_bytes())
        .map(|record| record.length())
        .with_context(|| format!("contig '{name}' is not in {fasta}"))
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

    #[test]
    fn small_window_is_one_window() {
        let windows = Region::parse("chr7:5,525,542-5,543,028")
            .unwrap()
            .windows(None)
            .unwrap();
        assert_eq!(windows, vec![Window::new("chr7", 5_525_542, 5_543_028)]);
        assert_eq!(windows[0].label, "chr7_5525542_5543028");
    }

    #[test]
    fn large_window_is_split_as_in_python() {
        // utils.make_windows("chr1:0-2500000"): 3 chunks of ceil(2500000 / 3)
        let windows = Region::parse("chr1:0-2500000")
            .unwrap()
            .windows(None)
            .unwrap();
        let bounds: Vec<(u64, u64)> = windows.iter().map(|w| (w.begin, w.end)).collect();
        assert_eq!(
            bounds,
            vec![(0, 833_334), (833_334, 1_666_668), (1_666_668, 2_500_000)]
        );
    }

    #[test]
    fn entire_contig_needs_fasta() {
        assert!(Region::parse("chrS").unwrap().windows(None).is_err());
        let windows = Region::parse("chrS")
            .unwrap()
            .windows(Some("test-data/ref.fa"))
            .unwrap();
        assert_eq!(windows.len(), 1);
        assert_eq!((windows[0].begin, windows[0].end), (0, 30_000));
        assert_eq!(windows[0].label, "chrS");
    }
}
