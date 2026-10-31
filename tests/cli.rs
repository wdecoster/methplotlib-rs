//! End-to-end runs of the binary on the test data

use std::path::PathBuf;
use std::process::{Command, Output};

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_methplotlib"))
        .args(args)
        .output()
        .expect("failed to run methplotlib")
}

fn out_dir(name: &str) -> PathBuf {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(name);
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn assert_success(output: &Output) {
    assert!(
        output.status.success(),
        "methplotlib failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn all_input_types_in_every_output_format() {
    let dir = out_dir("all_formats");
    for ext in ["svg", "png", "pdf"] {
        let outfile = dir.join(format!("{{region}}/browser.{ext}"));
        let output = run(&[
            "-m",
            "test-data/reads.bam",
            "test-data/reads.cram",
            "test-data/modkit.bed.gz",
            "test-data/modbam2bed.bed.gz",
            "test-data/values.bedgraph.gz",
            "-n",
            "bam",
            "cram",
            "modkit",
            "modbam2bed",
            "bedgraph",
            "-w",
            "chrS:5,000-25,000",
            "-f",
            "test-data/ref.fa",
            "--mods",
            "C+m,5mC",
            "-o",
            outfile.to_str().unwrap(),
        ]);
        assert_success(&output);
        let written = dir.join(format!("chrS_5000_25000/browser.{ext}"));
        assert!(
            std::fs::metadata(&written).unwrap().len() > 1000,
            "{written:?} is too small"
        );
    }
    let svg = std::fs::read_to_string(dir.join("chrS_5000_25000/browser.svg")).unwrap();
    for title in [
        "bam_C+m",
        "cram_C+m",
        "modkit_m",
        "modbam2bed_5mC",
        "bedgraph",
    ] {
        assert!(svg.contains(title), "{title} missing from the svg");
    }
}

#[test]
fn overlay_with_annotation_and_entire_contig() {
    let dir = out_dir("overlay");
    let output = run(&[
        "-m",
        "test-data/modkit_only_tabs.bed",
        "test-data/modkit.bedgraph",
        "-n",
        "frequency",
        "values",
        "-w",
        "chrS",
        "-f",
        "test-data/ref.fa",
        "-o",
        dir.join("{region}.svg").to_str().unwrap(),
    ]);
    assert_success(&output);
    assert!(dir.join("chrS.svg").is_file());
}

#[test]
fn errors_are_reported() {
    let dir = out_dir("errors");
    let outfile = dir.join("out.svg");
    let outfile = outfile.to_str().unwrap();
    let cases: [(&[&str], &str); 5] = [
        (
            &[
                "-m",
                "test-data/reads.bam",
                "-n",
                "a",
                "b",
                "-w",
                "chrS:1-100",
            ],
            "same number of --names",
        ),
        (
            &[
                "-m",
                "test-data/reads.bam",
                "-n",
                "a",
                "-w",
                "chrS:1-100",
                "-o",
                "out.html",
            ],
            "unsupported output format",
        ),
        (
            &[
                "-m",
                "test-data/ref.fa",
                "-n",
                "a",
                "-w",
                "chrS:1-100",
                "-o",
                outfile,
            ],
            "was not recognized",
        ),
        (
            &[
                "-m",
                "test-data/reads.bam",
                "-n",
                "a",
                "-w",
                "chrS:5000-25000",
                "--mods",
                "C+x",
                "-o",
                outfile,
            ],
            "C+h, C+h., C+m, C+m.",
        ),
        (
            &[
                "-m",
                "test-data/values.bedgraph.gz",
                "-n",
                "a",
                "-w",
                "chrX:1-100",
                "-o",
                outfile,
            ],
            "no data found",
        ),
    ];
    for (args, message) in cases {
        let output = run(args);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(!output.status.success(), "{args:?} should fail");
        assert!(
            stderr.contains(message),
            "{args:?}: expected '{message}' in\n{stderr}"
        );
    }
}
