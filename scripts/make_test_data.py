"""
Generate the synthetic test data in test-data/, and the results of the Python methplotlib for it
in test-data/expected/, which the Rust tests compare against.

The Python methplotlib is the reference for the port: the expected results are what it parses
from these files, so a Rust test failing on them is a difference in behaviour.

Run from the repository root, with the Python methplotlib (and so pysam) installed:
    python scripts/make_test_data.py
All data is synthetic, apart from the annotation copied from the examples of the Python methplotlib.
"""

import array
import gzip
import json
import math
import random
import shutil
from pathlib import Path

import numpy as np
import pysam
from methplotlib import annotation, plots
from methplotlib.import_methylation import parse_bedgraph, parse_bedmethyl, parse_cram
from methplotlib.utils import Region

OUT = Path("test-data")
EXPECTED = OUT / "expected"
CONTIG = "chrS"
CONTIG_LENGTH = 30_000
WINDOW = f"{CONTIG}:5000-25000"
rng = random.Random(1)
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def random_sequence(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def clean(value):
    """numpy scalars and NaN to plain json"""
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return None if math.isnan(value) else float(value)
    return value


def write_json(name, data):
    with open(EXPECTED / f"{name}.json", "w") as handle:
        json.dump(data, handle, separators=(",", ":"), allow_nan=False)


# --- reference -----------------------------------------------------------------------------------
REFERENCE = random_sequence(CONTIG_LENGTH)


def write_reference():
    with open(OUT / "ref.fa", "w") as handle:
        handle.write(f">{CONTIG}\n")
        for i in range(0, CONTIG_LENGTH, 60):
            handle.write(REFERENCE[i : i + 60] + "\n")
    pysam.faidx(str(OUT / "ref.fa"))


# --- reads ---------------------------------------------------------------------------------------
def make_alignment(start, length, complex_cigar):
    """A query sequence and CIGAR, aligned from start over about length reference bases"""
    ops = []
    query = []
    if rng.random() < 0.3:
        ops.append((5, rng.randrange(1, 30)))  # hard clip
    clip = rng.randrange(0, 150)
    if clip:
        ops.append((4, clip))
        query.append(random_sequence(clip))
    pos = start
    end = start + length
    while pos < end:
        match = min(rng.randrange(20, 400) if complex_cigar else end - pos, end - pos)
        ops.append((0, match))
        query.append(REFERENCE[pos : pos + match])
        pos += match
        if complex_cigar and pos < end - 10:
            kind = rng.choice([1, 2, 2, 3])  # insertion, deletion, skip
            n = rng.randrange(1, 8) if kind != 3 else rng.randrange(20, 60)
            ops.append((kind, n))
            if kind == 1:
                query.append(random_sequence(n))
            else:
                pos += n
    clip = rng.randrange(0, 150)
    if clip:
        ops.append((4, clip))
        query.append(random_sequence(clip))
    if rng.random() < 0.3:
        ops.append((5, rng.randrange(1, 30)))
    return "".join(query), ops


def modification_tags(forward, style):
    """MM/ML tags for every CpG C of the read as sequenced"""
    c_index = [i for i, base in enumerate(forward) if base == "C"]
    cpg = [k for k, i in enumerate(c_index) if forward[i + 1 : i + 2] == "G"]
    if not cpg:
        return None
    deltas = ",".join(str(d) for d in np.diff([-1] + cpg) - 1)
    probabilities = [rng.randrange(0, 256) for _ in range(2 * len(cpg))]
    if style == "one_mod":
        return f"C+m?,{deltas};", probabilities[: len(cpg)]
    if style == "dot":
        return f"C+m.,{deltas};C+h.,{deltas};", probabilities
    if style == "combined":
        # both modifications in one context, with their probabilities interleaved per position
        return f"C+mh?,{deltas};", probabilities
    return f"C+m?,{deltas};C+h?,{deltas};", probabilities


def make_reads():
    reads = []
    specs = [(rng.randrange(0, CONTIG_LENGTH - 3000), rng.randrange(500, 3000)) for _ in range(40)]
    # a stack of short reads, deeper than the 100 rows that are shown,
    # with equal starts and different lengths, as the longest read is placed first
    specs += [(20_000 + rng.randrange(0, 50), rng.randrange(100, 250)) for _ in range(120)]
    # reads starting where another ends, which go in another row
    specs += [(12_000, 1000), (13_000, 1000), (14_000, 500)]
    for i, (start, length) in enumerate(specs):
        query, ops = make_alignment(start, length, complex_cigar=i < 40)
        reverse = rng.random() < 0.5
        forward = query[::-1].translate(COMPLEMENT) if reverse else query
        style = rng.choice(["two_mods"] * 6 + ["one_mod", "dot", "combined", "legacy_tags", "no_ml", "no_mm"])
        tags = None if style == "no_mm" else modification_tags(forward, style)
        reads.append(dict(name=f"read{i:03d}", start=start, query=query, ops=ops,
                          reverse=reverse, tags=tags, style=style, flag=0))
    # secondary and supplementary alignments are skipped, so a copy of a read as each
    for flag in (256, 2048):
        copy = dict(reads[0])
        copy["flag"] = flag
        reads.append(copy)
    return sorted(reads, key=lambda r: r["start"])


def remove_reference_path(cram_path):
    """
    htslib writes the absolute path of the reference in the header of a cram (UR),
    which is local to the machine that made the test data, and not needed with the M5 checksum
    """
    with pysam.AlignmentFile(str(cram_path), "rc", reference_filename=str(OUT / "ref.fa")) as cram:
        lines = str(cram.header).splitlines()
    header_path = OUT / "header.sam"
    header_path.write_text(
        "".join("\t".join(f for f in line.split("\t") if not f.startswith("UR:")) + "\n" for line in lines)
    )
    pysam.samtools.reheader("-i", str(header_path), str(cram_path), catch_stdout=False)
    header_path.unlink()


def write_alignments(reads):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": CONTIG, "LN": CONTIG_LENGTH}]}
    bam_path, cram_path = OUT / "reads.bam", OUT / "reads.cram"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam, pysam.AlignmentFile(
        str(cram_path), "wc", header=header, reference_filename=str(OUT / "ref.fa")
    ) as cram:
        for r in reads:
            a = pysam.AlignedSegment()
            a.query_name = r["name"]
            a.query_sequence = r["query"]
            a.flag = r["flag"] | (16 if r["reverse"] else 0)
            a.reference_id = 0
            a.reference_start = r["start"]
            a.mapping_quality = 60
            a.cigartuples = r["ops"]
            a.query_qualities = pysam.qualitystring_to_array("I" * len(r["query"]))
            if r["tags"]:
                mm, ml = r["tags"]
                modtag, qualtag = ("Mm", "Ml") if r["style"] == "legacy_tags" else ("MM", "ML")
                a.set_tag(modtag, mm)
                if r["style"] != "no_ml":
                    a.set_tag(qualtag, array.array("B", ml))
            bam.write(a)
            cram.write(a)
    pysam.index(str(bam_path))
    remove_reference_path(cram_path)
    pysam.index(str(cram_path))
    return bam_path


def expected_alignments(bam_path):
    window = Region(WINDOW)
    result = {}
    for mods in (None, "C+h"):
        datasets = parse_cram(str(bam_path), "bam", "reads", window, mods_of_interest=mods)
        heights = plots.assign_y_height_per_read(datasets[0].start_end_table, max_coverage=100)
        result[mods or "all"] = dict(
            datasets=[
                dict(
                    name=d.name,
                    called_sites=int(d.called_sites),
                    rows=[[r.read_name, r.strand, clean(r.pos), clean(r.quality)]
                          for r in d.table.itertuples()],
                )
                for d in datasets
            ],
            start_end=[[name, int(r.posmin), int(r.posmax), r.strand]
                       for name, r in datasets[0].start_end_table.iterrows()],
            heights={read: int(h) for read, h in heights["height"].items()},
        )
    write_json("reads", result)


# --- bedMethyl and bedgraph ----------------------------------------------------------------------
def cpg_positions():
    """every fourth CpG, to keep the files small"""
    return [i for i in range(CONTIG_LENGTH - 1) if REFERENCE[i : i + 2] == "CG"][::4]


def bedmethyl_records(flavor):
    """Records for both strands of every CpG, for 5mC and 5hmC, sorted as tabix needs them"""
    records = []
    for contig, positions in ((CONTIG, cpg_positions()), (CONTIG + "_alt", [10, 200, 9000])):
        for pos in positions:
            for strand, start in (("+", pos), ("-", pos + 1)):
                for mod in ("h", "m"):
                    valid = rng.randrange(3, 30)
                    modified = rng.randrange(0, valid + 1)
                    canonical = valid - modified
                    if flavor == "modkit":
                        fraction = 100 * modified / valid
                        bed9 = [contig, start, start + 1, mod, valid, strand, start, start + 1, "255,0,0"]
                        rest = [valid, f"{fraction:.2f}", modified, canonical, 0, 0, 0, 0, 0]
                    else:
                        bed9 = [contig, start, start + 1, "5mC" if mod == "m" else "5hmC", valid, strand,
                                start, start + 1, "0,0,0"]
                        rest = [valid, f"{100 * modified / valid:.2f}", canonical, modified, 0]
                    records.append((contig, start, bed9, rest))
    return [(b, r) for _, _, b, r in sorted(records, key=lambda x: (x[0], x[1]))]


def write_text(path, lines, compression):
    """compression is None, 'gzip' (not indexable) or 'bgzip' (indexed with tabix later)"""
    text = "".join(line + "\n" for line in lines)
    if compression == "gzip":
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        with open(path, "w") as handle:
            handle.write(text)


def bgzip_and_index(path):
    pysam.tabix_compress(str(path), str(path) + ".gz", force=True)
    path.unlink()
    pysam.tabix_index(str(path) + ".gz", preset="bed", force=True)


def make_bedmethyl():
    window = Region(WINDOW)
    modkit = bedmethyl_records("modkit")
    modbam2bed = bedmethyl_records("modbam2bed")
    files = {
        # modkit separates the columns after BED9 by spaces, unless run with --only-tabs
        "modkit.bed": ([("\t".join(map(str, b)) + "\t" + " ".join(map(str, r))) for b, r in modkit], "bgzip"),
        # with a comment, which is not data
        "modkit_only_tabs.bed": (["#chrom\tstart\tend"] + ["\t".join(map(str, b + r)) for b, r in modkit], None),
        "modbam2bed.bed.gz": (["\t".join(map(str, b + r)) for b, r in modbam2bed], "gzip"),
    }
    expected = {}
    for filename, (lines, compression) in files.items():
        path = OUT / filename
        write_text(path, lines, compression)
        flavor = "modbam2bed" if filename.startswith("modbam2bed") else "modkit"
        final = filename + ".gz" if compression == "bgzip" else filename
        for smooth in (1, 4, 5):
            # modbam2bed names the modification in full, modkit with the code of the MM tag
            for mods in (None, "5mC" if flavor == "modbam2bed" else "C+m"):
                datasets = parse_bedmethyl(str(path), "sample", window, smoothen=smooth, flavor=flavor,
                                           mods_of_interest=mods)
                expected[f"{final}|smooth={smooth}|mods={mods}"] = [
                    dict(name=d.name, called_sites=int(d.called_sites),
                         start=[int(s) for s in d.table["Start"]],
                         frequency=[clean(f) for f in d.table["modified_frequency"]])
                    for d in datasets
                ]
        if compression == "bgzip":
            bgzip_and_index(path)
    write_json("bedmethyl", expected)


def make_bedgraph():
    window = Region(WINDOW)
    records = []
    for contig in (CONTIG, CONTIG + "_alt"):
        start = 0
        while start < (CONTIG_LENGTH if contig == CONTIG else 3000):
            end = start + rng.randrange(1, 200)
            records.append((contig, start, end, round(rng.random(), 4), rng.randrange(1, 40)))
            start = end + rng.randrange(0, 100)
    files = {
        "values.bedgraph": (["\t".join(map(str, r[:4])) for r in records], "bgzip"),
        # `modkit pileup --bedgraph` writes the coverage as a fifth column
        # and with a track line, which is not data
        "modkit.bedgraph": (["track type=bedGraph"] + ["\t".join(map(str, r)) for r in records], None),
    }
    expected = {}
    for filename, (lines, compression) in files.items():
        path = OUT / filename
        write_text(path, lines, compression)
        datasets = parse_bedgraph(str(path), "sample", window)
        expected[filename + ".gz" if compression == "bgzip" else filename] = [
            dict(name=d.name, called_sites=int(d.called_sites),
                 start=[int(s) for s in d.table["Start"]], end=[int(e) for e in d.table["End"]],
                 value=[clean(v) for v in d.table["Value"]])
            for d in datasets
        ]
        if compression == "bgzip":
            bgzip_and_index(path)
    write_json("bedgraph", expected)


# --- annotation ----------------------------------------------------------------------------------
def make_annotation(python_examples):
    for name in ("GRCh38-ACTB-locus.gtf.gz", "DNase_cluster_ACTB.bed.gz"):
        shutil.copy(python_examples / name, OUT / name)
    # a gff3 version of the gtf, as the attributes are parsed differently
    with gzip.open(OUT / "GRCh38-ACTB-locus.gtf.gz", "rt") as gtf, gzip.open(
        OUT / "GRCh38-ACTB-locus.gff.gz", "wt"
    ) as gff:
        for line in gtf:
            fields = line.rstrip("\n").split("\t")
            attributes = [a.strip().split(" ", 1) for a in fields[8].rstrip(";").split(";") if a.strip()]
            fields[8] = ";".join(f"{k}={v.strip(chr(34))}" for k, v in attributes)
            gff.write("\t".join(fields) + "\n")
    expected = {}
    for window in ("chr7:5,525,542-5,543,028", "chr7:5,300,000-5,600,000", "chr1:1-100000"):
        for filename in ("GRCh38-ACTB-locus.gtf.gz", "GRCh38-ACTB-locus.gff.gz"):
            for simplify in (False, True):
                transcripts = annotation.parse_annotation(str(OUT / filename), Region(window), simplify)
                expected[f"{filename}|{window}|simplify={simplify}"] = [
                    dict(transcript=t.transcript, gene=t.gene, strand=t.strand, color=t.color,
                         begin=int(t.begin), end=int(t.end),
                         exons=[[int(b), int(e)] for b, e in t.exon_tuples])
                    for t in transcripts
                ]
        expected[f"DNase_cluster_ACTB.bed.gz|{window}"] = [
            [int(b), int(e), str(n)]
            for b, e, n in annotation.parse_bed(str(OUT / "DNase_cluster_ACTB.bed.gz"), Region(window))
        ]
    write_json("annotation", expected)


def main():
    import methplotlib

    python_examples = Path(methplotlib.__file__).parents[1] / "examples"
    OUT.mkdir(exist_ok=True)
    EXPECTED.mkdir(exist_ok=True)
    write_reference()
    expected_alignments(write_alignments(make_reads()))
    make_bedmethyl()
    make_bedgraph()
    make_annotation(python_examples)


if __name__ == "__main__":
    main()
