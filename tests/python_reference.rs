//! The parsed data must match what the Python methplotlib parses from the same files,
//! as stored in test-data/expected by scripts/make_test_data.py

use std::collections::HashMap;

use methplotlib_rs::annotation::{parse_annotation, parse_bed};
use methplotlib_rs::frequency::{parse_bedgraph, parse_bedmethyl};
use methplotlib_rs::reads::{assign_rows, parse_alignments};
use methplotlib_rs::region::{Region, Window};
use methplotlib_rs::sniff::file_type;
use serde_json::Value;

const WINDOW: &str = "chrS:5000-25000";

fn expected(name: &str) -> Value {
    let path = format!("test-data/expected/{name}.json");
    serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap()
}

fn window(region: &str) -> Window {
    Region::parse(region)
        .unwrap()
        .windows(None)
        .unwrap()
        .remove(0)
}

/// NaN is null in the json
fn assert_floats(actual: &[f64], expected: &Value, context: &str) {
    let expected = expected.as_array().unwrap();
    assert_eq!(actual.len(), expected.len(), "{context}: length");
    for (i, (a, e)) in actual.iter().zip(expected).enumerate() {
        match e.as_f64() {
            None => assert!(a.is_nan(), "{context}[{i}]: {a} should be NaN"),
            Some(e) => assert!((a - e).abs() < 1e-9, "{context}[{i}]: {a} != {e}"),
        }
    }
}

fn u64s(v: &Value) -> Vec<u64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_u64().unwrap())
        .collect()
}

#[test]
fn bedmethyl_matches_python() {
    let expected = expected("bedmethyl");
    let cases = expected.as_object().unwrap();
    assert_eq!(cases.len(), 18);
    for (key, datasets) in cases {
        let parts: Vec<&str> = key.split('|').collect();
        let path = format!("test-data/{}", parts[0]);
        let smooth: usize = parts[1].trim_start_matches("smooth=").parse().unwrap();
        let mods = parts[2].trim_start_matches("mods=");
        let mods = (mods != "None").then_some(mods);
        let parsed = parse_bedmethyl(
            &path,
            "sample",
            &window(WINDOW),
            smooth,
            file_type(&path).unwrap(),
            mods,
        )
        .unwrap();
        let datasets = datasets.as_array().unwrap();
        assert_eq!(parsed.len(), datasets.len(), "{key}: number of datasets");
        for (p, e) in parsed.iter().zip(datasets) {
            assert_eq!(p.name, e["name"].as_str().unwrap(), "{key}");
            assert_eq!(p.start, u64s(&e["start"]), "{key}: start");
            assert_floats(&p.frequency, &e["frequency"], &format!("{key} {}", p.name));
        }
    }
}

#[test]
fn bedgraph_matches_python() {
    for (file, datasets) in expected("bedgraph").as_object().unwrap() {
        let path = format!("test-data/{file}");
        let parsed = parse_bedgraph(&path, "sample", &window(WINDOW))
            .unwrap()
            .unwrap();
        let e = &datasets[0];
        assert_eq!(parsed.name, e["name"].as_str().unwrap());
        assert_eq!(parsed.start, u64s(&e["start"]), "{file}: start");
        assert_eq!(parsed.end, u64s(&e["end"]), "{file}: end");
        assert_floats(&parsed.value, &e["value"], file);
    }
}

#[test]
fn bam_and_cram_match_python() {
    let expected = expected("reads");
    for path in ["test-data/reads.bam", "test-data/reads.cram"] {
        for (mods_key, case) in expected.as_object().unwrap() {
            let mods = (mods_key != "all").then_some(mods_key.as_str());
            let context = format!("{path} mods={mods_key}");
            let parsed = parse_alignments(
                path,
                "reads",
                &window(WINDOW),
                mods,
                Some("test-data/ref.fa"),
                path.ends_with(".cram"),
            )
            .unwrap();
            let datasets = case["datasets"].as_array().unwrap();
            assert_eq!(
                parsed.len(),
                datasets.len(),
                "{context}: number of datasets"
            );
            for (p, e) in parsed.iter().zip(datasets) {
                assert_eq!(p.name, e["name"].as_str().unwrap(), "{context}");
                let rows = e["rows"].as_array().unwrap();
                assert_eq!(
                    p.calls.len(),
                    rows.len(),
                    "{context} {}: number of calls",
                    p.name
                );
                for (call, row) in p.calls.iter().zip(rows) {
                    let read = &p.reads[call.read];
                    assert_eq!(read.name, row[0].as_str().unwrap(), "{context}: read name");
                    assert_eq!(
                        if read.reverse { "-" } else { "+" },
                        row[1].as_str().unwrap()
                    );
                    assert_eq!(
                        call.pos,
                        row[2].as_f64().map(|p| p as u64),
                        "{context}: {}",
                        read.name
                    );
                    assert_eq!(
                        f64::from(call.quality),
                        row[3].as_f64().unwrap(),
                        "{context}: {}",
                        read.name
                    );
                }
            }
            // all primary reads in the window, and the rows they are plotted in
            let reads = &parsed[0].reads;
            let start_end: Vec<(String, u64, u64, &str)> = reads
                .iter()
                .map(|r| {
                    let strand = if r.reverse { "-" } else { "+" };
                    (r.name.clone(), r.start, r.end, strand)
                })
                .collect();
            let expected_start_end: Vec<(String, u64, u64, &str)> = case["start_end"]
                .as_array()
                .unwrap()
                .iter()
                .map(|r| {
                    (
                        r[0].as_str().unwrap().to_string(),
                        r[1].as_u64().unwrap(),
                        r[2].as_u64().unwrap(),
                        r[3].as_str().unwrap(),
                    )
                })
                .collect();
            assert_eq!(start_end, expected_start_end, "{context}: reads");
            let rows: HashMap<String, usize> = reads
                .iter()
                .zip(assign_rows(reads, 100))
                .filter_map(|(r, row)| row.map(|row| (r.name.clone(), row)))
                .collect();
            let expected_rows: HashMap<String, usize> = case["heights"]
                .as_object()
                .unwrap()
                .iter()
                .map(|(k, v)| (k.clone(), v.as_u64().unwrap() as usize))
                .collect();
            assert_eq!(rows, expected_rows, "{context}: rows");
        }
    }
}

#[test]
fn annotation_matches_python() {
    for (key, expected) in expected("annotation").as_object().unwrap() {
        let parts: Vec<&str> = key.split('|').collect();
        let path = format!("test-data/{}", parts[0]);
        let window = window(parts[1]);
        if parts[0].ends_with(".bed.gz") {
            let parsed = parse_bed(&path, &window).unwrap();
            let mut expected: Vec<(u64, u64, String)> = expected
                .as_array()
                .unwrap()
                .iter()
                .map(|r| {
                    (
                        r[0].as_u64().unwrap(),
                        r[1].as_u64().unwrap(),
                        r[2].as_str().unwrap().to_string(),
                    )
                })
                .collect();
            expected.sort();
            let parsed: Vec<(u64, u64, String)> = parsed
                .into_iter()
                .map(|r| (r.begin, r.end, r.name))
                .collect();
            assert_eq!(parsed, expected, "{key}");
            continue;
        }
        let simplify = parts[2] == "simplify=True";
        let parsed = parse_annotation(&path, &window, simplify).unwrap();
        let expected = expected.as_array().unwrap();
        assert_eq!(parsed.len(), expected.len(), "{key}: number of transcripts");
        for (p, e) in parsed.iter().zip(expected) {
            assert_eq!(p.transcript, e["transcript"].as_str().unwrap(), "{key}");
            assert_eq!(p.gene, e["gene"].as_str().unwrap_or_default(), "{key}");
            assert_eq!(p.strand, e["strand"].as_str().unwrap(), "{key}");
            assert_eq!(p.color, e["color"].as_str().unwrap(), "{key}");
            assert_eq!(
                (p.begin, p.end),
                (e["begin"].as_u64().unwrap(), e["end"].as_u64().unwrap()),
                "{key}"
            );
            // pandas sorts the exons by begin with an unstable sort, so ties may be in any order
            let mut exons = p.exons.clone();
            exons.sort();
            let mut expected_exons: Vec<(u64, u64)> = e["exons"]
                .as_array()
                .unwrap()
                .iter()
                .map(|x| (x[0].as_u64().unwrap(), x[1].as_u64().unwrap()))
                .collect();
            expected_exons.sort();
            assert_eq!(exons, expected_exons, "{key} {}", p.transcript);
        }
    }
}
