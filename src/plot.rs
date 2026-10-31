//! The genome-browser view: a track per dataset (or one overlaid track) and the annotation,
//! sharing one genomic x-axis, drawn with kuva's TrackStack

use std::sync::atomic::{AtomicUsize, Ordering};

use anyhow::{Context, Result, bail};
use kuva::plot::ScatterPlot;
use kuva::plot::legend::{LegendEntry, LegendShape};
use kuva::render::color::Color;
use kuva::render::layout::{ComputedLayout, Layout, TickFormat};
use kuva::render::render::{PathData, Primitive, Scene, TextAnchor, render_multiple};
use kuva::render::track_stack::{AxisSpec, Track, TrackCtx, TrackHeight, TrackStack};
use log::warn;

use crate::annotation::{BedRegion, Transcript};
use crate::frequency::{Bedgraph, Frequencies};
use crate::reads::{ReadCalls, assign_rows};
use crate::region::Window;

/// Reads beyond this number of rows are hidden, as in the Python version
const MAX_ROWS: usize = 100;
const WIDTH: f64 = 1400.0;

/// plotly's default colorway, which the Python version's frequency lines are colored with
const COLORWAY: [&str; 10] = [
    "#636efa", "#EF553B", "#00cc96", "#ab63fa", "#FFA15A", "#19d3f3", "#FF6692", "#B6E880",
    "#FF97FF", "#FECB52",
];

/// The ColorBrewer Reds scale, as plotly's "Reds" with which the Python version colors the calls
const REDS: [(f64, (u8, u8, u8)); 9] = [
    (0.0, (255, 245, 240)),
    (0.125, (254, 224, 210)),
    (0.25, (252, 187, 161)),
    (0.375, (252, 146, 114)),
    (0.5, (251, 106, 74)),
    (0.625, (239, 59, 44)),
    (0.75, (203, 24, 29)),
    (0.875, (165, 15, 21)),
    (1.0, (103, 0, 13)),
];

pub enum Dataset {
    Reads(ReadCalls),
    Frequencies(Frequencies),
    Bedgraph(Bedgraph),
}

impl Dataset {
    pub fn name(&self) -> &str {
        match self {
            Dataset::Reads(d) => &d.name,
            Dataset::Frequencies(d) => &d.name,
            Dataset::Bedgraph(d) => &d.name,
        }
    }
}

pub struct PlotOptions {
    pub split: bool,
    pub dotsize: f64,
    pub minqual: u8,
}

fn reds(t: f64) -> String {
    let t = t.clamp(0.0, 1.0);
    let i = REDS
        .iter()
        .position(|(stop, _)| *stop >= t)
        .unwrap_or(REDS.len() - 1)
        .max(1);
    let (t0, c0) = REDS[i - 1];
    let (t1, c1) = REDS[i];
    let f = if t1 > t0 { (t - t0) / (t1 - t0) } else { 0.0 };
    let mix = |a: u8, b: u8| (f64::from(a) + f * (f64::from(b) - f64::from(a))).round() as u8;
    format!(
        "rgb({},{},{})",
        mix(c0.0, c1.0),
        mix(c0.1, c1.1),
        mix(c0.2, c1.2)
    )
}

/// A triangle pointing right (forward strand) or left, centered on (x, y)
fn triangle(x: f64, y: f64, size: f64, forward: bool, fill: &str) -> Primitive {
    let h = size / 2.0;
    let d = if forward {
        format!("M{},{}L{},{}L{},{}Z", x - h, y - h, x + h, y, x - h, y + h)
    } else {
        format!("M{},{}L{},{}L{},{}Z", x + h, y - h, x - h, y, x + h, y + h)
    };
    Primitive::Path(Box::new(PathData {
        d,
        fill: Some(Color::from(fill)),
        stroke: Color::from("black"),
        stroke_width: 0.5,
        opacity: None,
        stroke_dasharray: None,
    }))
}

fn line(x1: f64, y1: f64, x2: f64, y2: f64, color: &str, width: f64) -> Primitive {
    Primitive::Line {
        x1,
        y1,
        x2,
        y2,
        stroke: Color::from(color),
        stroke_width: width,
        stroke_dasharray: None,
    }
}

fn text(
    x: f64,
    y: f64,
    content: &str,
    size: u32,
    anchor: TextAnchor,
    rotate: Option<f64>,
) -> Primitive {
    Primitive::Text {
        x,
        y,
        content: content.to_string(),
        size,
        anchor,
        rotate,
        bold: false,
        color: None,
    }
}

/// Every clip region needs its own id in the svg
fn clip_id() -> String {
    static COUNTER: AtomicUsize = AtomicUsize::new(0);
    format!(
        "methplotlib-clip-{}",
        COUNTER.fetch_add(1, Ordering::Relaxed)
    )
}

type Draw = Box<dyn FnOnce(&ComputedLayout, &mut Vec<Primitive>)>;

/// A track with a y-axis drawn by kuva, and data drawn by `draw` in the coordinates of its layout,
/// clipped to the plot area
struct Panel {
    title: Option<String>,
    y_label: String,
    y_range: (f64, f64),
    y_integer_step: Option<f64>,
    height: f64,
    legend: Vec<LegendEntry>,
    /// the range of the probabilities, for a colorbar right of the panel
    colorbar: Option<(u8, u8)>,
    draw: Draw,
}

const COLORBAR_MARGIN: f64 = 160.0;

impl Panel {
    fn layout(&self) -> Layout {
        let mut layout = Layout::new((0.0, 1.0), self.y_range)
            .with_y_axis_min(self.y_range.0)
            .with_y_axis_max(self.y_range.1)
            .with_y_label(self.y_label.clone());
        if let Some(step) = self.y_integer_step {
            layout = layout
                .with_y_tick_step(step)
                .with_y_tick_format(TickFormat::Integer);
        }
        if let Some(title) = &self.title {
            layout = layout.with_title(title.clone()).with_title_size(13);
        }
        layout
    }
}

impl Track for Panel {
    fn height(&self) -> TrackHeight {
        TrackHeight::Fixed(self.height)
    }

    fn left_margin(&self) -> f64 {
        ComputedLayout::from_layout(&self.layout()).margin_left
    }

    fn right_margin(&self) -> f64 {
        let margin = ComputedLayout::from_layout(&self.layout()).margin_right;
        if self.colorbar.is_some() {
            margin.max(COLORBAR_MARGIN)
        } else {
            margin
        }
    }

    fn legend_entries(&self) -> Vec<LegendEntry> {
        self.legend.clone()
    }

    fn legend_group_title(&self) -> Option<String> {
        None
    }

    fn render(self: Box<Self>, cx: &TrackCtx<'_>, scene: &mut Scene) {
        let (x_min, x_max) = cx.x.x_range();
        let mut layout = self
            .layout()
            .with_width(cx.width)
            .with_height(cx.height)
            .with_x_axis_min(x_min)
            .with_x_axis_max(x_max)
            .with_force_margins(cx.x.px_left(), cx.width - cx.x.px_right());
        // the stack draws the one shared x-axis and the legend
        layout.suppress_x_ticks = true;
        layout.show_legend = false;
        layout.x_range = (x_min, x_max);
        let computed = ComputedLayout::from_layout(&layout);
        // kuva only draws axes for a plot, so an invisible one spanning the ranges
        let anchor = ScatterPlot::new()
            .with_data(vec![(x_min, self.y_range.0), (x_max, self.y_range.1)])
            .with_size(0.0);
        let axes = render_multiple(vec![anchor.into()], layout);
        scene.defs.extend(axes.defs);
        scene.add(Primitive::GroupStart {
            transform: Some(format!("translate(0,{})", cx.y_top)),
            title: None,
            extra_attrs: None,
        });
        for element in axes.elements {
            scene.add(element);
        }
        scene.add(Primitive::ClipStart {
            x: computed.margin_left,
            y: computed.margin_top,
            width: computed.plot_width(),
            height: computed.plot_height(),
            id: clip_id(),
        });
        let mut primitives = vec![];
        (self.draw)(&computed, &mut primitives);
        for p in primitives {
            scene.add(p);
        }
        scene.add(Primitive::ClipEnd);
        if let Some((low, high)) = self.colorbar {
            draw_colorbar(scene, &computed, cx.x.px_right(), low, high);
        }
        scene.add(Primitive::GroupEnd);
    }
}

/// A vertical colorbar of the probabilities, as the Python version shows for the first per-read panel
fn draw_colorbar(scene: &mut Scene, computed: &ComputedLayout, px_right: f64, low: u8, high: u8) {
    let x = px_right + 20.0;
    let top = computed.margin_top;
    let height = computed.plot_height();
    let steps = 50;
    for i in 0..steps {
        let t = f64::from(i) / f64::from(steps - 1);
        scene.add(Primitive::Rect {
            x,
            y: top + height * (1.0 - (f64::from(i) + 1.0) / f64::from(steps)),
            width: 14.0,
            height: height / f64::from(steps) + 0.5,
            fill: Color::from(reds(t).as_str()),
            stroke: None,
            stroke_width: None,
            opacity: None,
        });
    }
    scene.add(text(
        x + 18.0,
        top + 10.0,
        "Likely modified",
        10,
        TextAnchor::Start,
        None,
    ));
    scene.add(text(
        x + 18.0,
        top + height,
        "Likely unmodified",
        10,
        TextAnchor::Start,
        None,
    ));
    scene.add(text(
        x + 18.0,
        top + 24.0,
        &format!("({high})"),
        9,
        TextAnchor::Start,
        None,
    ));
    scene.add(text(
        x + 18.0,
        top + height - 14.0,
        &format!("({low})"),
        9,
        TextAnchor::Start,
        None,
    ));
    scene.add(text(
        x + 112.0,
        top + height / 2.0,
        "Modification probability",
        11,
        TextAnchor::Middle,
        Some(90.0),
    ));
}

fn reads_panel(data: ReadCalls, options: &PlotOptions, with_colorbar: bool) -> Panel {
    let rows = assign_rows(&data.reads, MAX_ROWS);
    let mut has_calls = vec![false; data.reads.len()];
    for call in &data.calls {
        has_calls[call.read] = true;
    }
    let hidden = (0..data.reads.len())
        .filter(|&i| has_calls[i] && rows[i].is_none())
        .count();
    if hidden > 0 {
        warn!(
            "{}: hiding {hidden} reads because coverage above {MAX_ROWS}x.",
            data.name
        );
    }
    let n_rows = rows.iter().flatten().max().copied().unwrap_or(1);
    let shown: Vec<(u64, usize, u8)> = data
        .calls
        .iter()
        .filter(|c| c.quality > options.minqual)
        .filter_map(|c| Some((c.pos?, rows[c.read]?, c.quality)))
        .collect();
    let low = shown.iter().map(|c| c.2).min().unwrap_or(0);
    let high = shown.iter().map(|c| c.2).max().unwrap_or(255);
    let dotsize = options.dotsize;
    let reads = data.reads;
    let draw: Draw = Box::new(move |layout, out| {
        for (read, row) in reads.iter().zip(&rows) {
            let Some(row) = row else { continue };
            let y = layout.map_y(*row as f64);
            let (x1, x2) = (
                layout.map_x(read.start as f64),
                layout.map_x(read.end as f64),
            );
            out.push(line(x1, y, x2, y, "lightgrey", 1.0));
            out.push(triangle(x1, y, 8.0, !read.reverse, "black"));
            out.push(triangle(x2, y, 8.0, !read.reverse, "black"));
        }
        // one batch of circles per probability, as each has its own color
        let mut per_quality: Vec<(Vec<f64>, Vec<f64>)> = vec![(vec![], vec![]); 256];
        for (pos, row, quality) in shown {
            per_quality[quality as usize]
                .0
                .push(layout.map_x(pos as f64));
            per_quality[quality as usize]
                .1
                .push(layout.map_y(row as f64));
        }
        for (quality, (cx, cy)) in per_quality.into_iter().enumerate() {
            if cx.is_empty() {
                continue;
            }
            let t = if high > low {
                (quality as f64 - f64::from(low)) / f64::from(high - low)
            } else {
                0.5
            };
            out.push(Primitive::CircleBatch {
                cx,
                cy,
                r: dotsize / 2.0,
                fill: Color::from(reds(t).as_str()),
                fill_opacity: None,
                stroke: None,
                stroke_width: None,
            });
        }
    });
    Panel {
        title: Some(data.name),
        y_label: "Reads".to_string(),
        y_range: (0.0, n_rows as f64 + 1.0),
        y_integer_step: Some(((n_rows as f64 / 5.0).ceil()).max(1.0)),
        height: (n_rows as f64 * 5.0).clamp(80.0, 500.0) + 50.0,
        legend: vec![],
        colorbar: with_colorbar.then_some((low, high)),
        draw,
    }
}

/// The points of a line, split where the value is unknown, as plotly leaves a gap there
fn segments(points: &[(f64, f64)]) -> Vec<Vec<(f64, f64)>> {
    points
        .split(|p| p.1.is_nan())
        .filter(|s| !s.is_empty())
        .map(<[_]>::to_vec)
        .collect()
}

/// Lines (with markers for frequencies) of one or more frequency and bedgraph datasets
fn lines_panel(datasets: Vec<(Dataset, usize)>, title: Option<String>, overlay: bool) -> Panel {
    let y_label = match datasets.first().map(|d| &d.0) {
        Some(Dataset::Bedgraph(_)) => "Value",
        _ => "Modified frequency",
    };
    let only_frequencies = datasets
        .iter()
        .all(|d| matches!(d.0, Dataset::Frequencies(_)));
    let y_range = if only_frequencies {
        (0.0, 1.0)
    } else {
        let values = datasets.iter().flat_map(|(d, _)| match d {
            Dataset::Bedgraph(b) => b.value.clone(),
            Dataset::Frequencies(f) => f.frequency.clone(),
            Dataset::Reads(_) => vec![],
        });
        let (low, high) = values
            .filter(|v| !v.is_nan())
            .fold((0.0_f64, f64::MIN), |(lo, hi), v| (lo.min(v), hi.max(v)));
        let high = if high > low { high } else { low + 1.0 };
        (low, high + 0.05 * (high - low))
    };
    let legend = if overlay {
        datasets
            .iter()
            .map(|(d, i)| LegendEntry {
                label: d.name().to_string(),
                color: COLORWAY[i % COLORWAY.len()].to_string(),
                shape: LegendShape::Line,
                dasharray: None,
            })
            .collect()
    } else {
        vec![]
    };
    let draw: Draw = Box::new(move |layout, out| {
        for (dataset, i) in datasets {
            let color = COLORWAY[i % COLORWAY.len()];
            let (points, markers): (Vec<(f64, f64)>, bool) = match dataset {
                Dataset::Frequencies(f) => (
                    f.start
                        .iter()
                        .map(|&s| s as f64)
                        .zip(f.frequency.iter().copied())
                        .collect(),
                    true,
                ),
                // a step for every interval, connected to the next interval as in the Python version
                Dataset::Bedgraph(b) => (
                    b.start
                        .iter()
                        .zip(&b.end)
                        .zip(&b.value)
                        .flat_map(|((&s, &e), &v)| [(s as f64, v), (e as f64, v)])
                        .collect(),
                    false,
                ),
                Dataset::Reads(_) => unreachable!("reads are drawn in their own panel"),
            };
            for segment in segments(&points) {
                let pixels: Vec<(f64, f64)> = segment
                    .iter()
                    .map(|&(x, y)| (layout.map_x(x), layout.map_y(y)))
                    .collect();
                if markers {
                    out.push(Primitive::CircleBatch {
                        cx: pixels.iter().map(|p| p.0).collect(),
                        cy: pixels.iter().map(|p| p.1).collect(),
                        r: 3.0,
                        fill: Color::from(color),
                        fill_opacity: None,
                        stroke: None,
                        stroke_width: None,
                    });
                }
                out.push(Primitive::PolyLine {
                    points: pixels,
                    stroke: Color::from(color),
                    stroke_width: 2.0,
                    stroke_dasharray: None,
                });
            }
        }
    });
    Panel {
        title,
        y_label: y_label.to_string(),
        y_range,
        y_integer_step: None,
        height: if overlay { 300.0 } else { 200.0 },
        legend,
        colorbar: None,
        draw,
    }
}

/// Transcripts (or genes) as a line with exon blocks and strand arrows, a row each,
/// and bed regions as grey blocks in the bottom row. Hovering shows the name in an svg viewer.
struct AnnotationTrack {
    transcripts: Vec<Transcript>,
    regions: Vec<BedRegion>,
    window: Window,
}

const ROW_PX: f64 = 22.0;

impl AnnotationTrack {
    fn rows(&self) -> usize {
        self.transcripts.len() + usize::from(!self.regions.is_empty())
    }
}

impl Track for AnnotationTrack {
    fn height(&self) -> TrackHeight {
        TrackHeight::Fixed(self.rows() as f64 * ROW_PX + 10.0)
    }

    fn label(&self) -> Option<&str> {
        None
    }

    fn render(self: Box<Self>, cx: &TrackCtx<'_>, scene: &mut Scene) {
        let (left, right) = (cx.x.px_left(), cx.x.px_right());
        scene.add(Primitive::ClipStart {
            x: left,
            y: cx.y_top,
            width: right - left,
            height: cx.height,
            id: clip_id(),
        });
        // the first transcript in the row above the bed regions, the last at the top
        let row_y = |row: usize| cx.y_bottom() - 5.0 - (row as f64 + 0.5) * ROW_PX;
        let first_transcript_row = usize::from(!self.regions.is_empty());
        for region in &self.regions {
            let y = row_y(0);
            scene.add(Primitive::GroupStart {
                transform: None,
                title: Some(format!("{} ({}-{})", region.name, region.begin, region.end)),
                extra_attrs: None,
            });
            scene.add(line(
                cx.x.map(region.begin as f64),
                y,
                cx.x.map(region.end as f64),
                y,
                "grey",
                16.0,
            ));
            scene.add(Primitive::GroupEnd);
        }
        for (i, transcript) in self.transcripts.iter().enumerate() {
            let y = row_y(first_transcript_row + i);
            let color = transcript.color.as_str();
            let forward = transcript.strand != "-";
            let begin = transcript.begin.max(self.window.begin) as f64;
            let end = transcript.end.min(self.window.end) as f64;
            let name = if transcript.gene == transcript.transcript || transcript.gene.is_empty() {
                transcript.transcript.clone()
            } else {
                format!("{} ({})", transcript.gene, transcript.transcript)
            };
            scene.add(Primitive::GroupStart {
                transform: None,
                title: Some(name.clone()),
                extra_attrs: None,
            });
            scene.add(line(cx.x.map(begin), y, cx.x.map(end), y, color, 2.0));
            // exons overlapping the window, rather than contained in it,
            // so that an exon larger than the window is shown as well
            for &(exon_begin, exon_end) in &transcript.exons {
                if exon_begin > self.window.end || exon_end < self.window.begin {
                    continue;
                }
                let (x1, x2) = (cx.x.map(exon_begin as f64), cx.x.map(exon_end as f64));
                scene.add(line(x1, y, x2, y, color, 8.0));
                scene.add(triangle(x1, y, 8.0, forward, color));
                scene.add(triangle(x2, y, 8.0, forward, color));
            }
            scene.add(Primitive::GroupEnd);
            // the name, as a static image has no hovering
            scene.add(Primitive::Text {
                x: cx.x.map(begin).max(left) + 2.0,
                y: y - 6.0,
                content: name,
                size: 10,
                anchor: TextAnchor::Start,
                rotate: None,
                bold: false,
                color: Some(Color::from(color)),
            });
        }
        scene.add(Primitive::ClipEnd);
    }
}

/// The figure for one window: a track per dataset with --split or when there are per-read data,
/// otherwise the datasets overlaid in one track, and the annotation at the bottom
pub fn make_scene(
    datasets: Vec<Dataset>,
    transcripts: Vec<Transcript>,
    regions: Vec<BedRegion>,
    window: &Window,
    options: &PlotOptions,
) -> Scene {
    let split = options.split || datasets.iter().any(|d| matches!(d, Dataset::Reads(_)));
    let mut stack = TrackStack::new()
        .with_title("Nucleotide modifications")
        .x_range(window.begin as f64, window.end as f64)
        .spacing(6.0);
    if split {
        let mut colorbar_drawn = false;
        for (i, dataset) in datasets.into_iter().enumerate() {
            match dataset {
                Dataset::Reads(data) => {
                    stack = stack.track(reads_panel(data, options, !colorbar_drawn));
                    colorbar_drawn = true;
                }
                other => {
                    let title = Some(other.name().to_string());
                    stack = stack.track(lines_panel(vec![(other, i)], title, false));
                }
            }
        }
    } else {
        let indexed = datasets
            .into_iter()
            .enumerate()
            .map(|(i, d)| (d, i))
            .collect();
        stack = stack.track(lines_panel(indexed, None, true));
    }
    if !transcripts.is_empty() || !regions.is_empty() {
        stack = stack.track(AnnotationTrack {
            transcripts,
            regions,
            window: window.clone(),
        });
    }
    stack
        .x_axis_with(AxisSpec::genomic(window.chromosome.clone()))
        .render(WIDTH)
}

/// Write the figure as svg, png or pdf, as set by the extension of the path
pub fn write_scene(scene: &Scene, path: &str) -> Result<()> {
    let lower = path.to_lowercase();
    let bytes = if lower.ends_with(".svg") {
        kuva::backend::svg::SvgBackend::new()
            .render_scene(scene)
            .into_bytes()
    } else if lower.ends_with(".png") {
        kuva::backend::png::PngBackend::new()
            .with_scale(2.0)
            .render_scene(scene)
            .map_err(anyhow::Error::msg)?
    } else if lower.ends_with(".pdf") {
        kuva::backend::pdf::PdfBackend::new()
            .render_scene(scene)
            .map_err(anyhow::Error::msg)?
    } else {
        bail!("unsupported output format for {path}: use a .svg, .png or .pdf extension");
    };
    if let Some(parent) = std::path::Path::new(path).parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("failed to create the directory {}", parent.display()))?;
    }
    std::fs::write(path, bytes).with_context(|| format!("failed to write {path}"))
}
