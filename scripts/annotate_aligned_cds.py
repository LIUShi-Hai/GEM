# annotate_aligned_cds.py
"""
Annotate aligned CDS regions with Prokka for GEM output.

REQUIRES:
- Python >=3.8
- prokka (via conda)
- biopython
- bcbio-gff

To install:
    conda create -n gem_cds_env python=3.10 prokka biopython -c bioconda -c conda-forge
    conda activate gem_cds_env
    conda install -c bioconda -c conda-forge bcbio-gf
    
Usage:
    python annotate_aligned_cds.py --input-csv-dir ./gem_output --output-dir ./gem_output --novel-fasta ./test/novel.fasta --threads 4

nohup python annotate_aligned_cds.py --input-csv-dir ./gem_output --output-dir ./gem_output --novel-fasta ./test/novel.fasta --threads 4 > nohup.out 2>&1 &


This script:
- Finds unique qseqids in GEM output across Species_link_Genetic_Exchange_Prediction_d{d}.csv files
- Extracts matching sequences from novel.fasta
- Runs Prokka once per unique sequence
- Parses overlapping CDS features
- Outputs annotated results with pair number, gene name, and product
"""

import os
import csv
import sys
import argparse
import subprocess
from pathlib import Path
from Bio import SeqIO
from BCBio import GFF
from concurrent.futures import ThreadPoolExecutor

# Setup logging
log_file = open("annotate_aligned_cds.log", "w", buffering=1)
class TeeLogger:
    def write(self, msg):
        sys.__stdout__.write(msg)
        log_file.write(msg)
    def flush(self):
        sys.__stdout__.flush()
        log_file.flush()
sys.stdout = TeeLogger()
sys.stderr = TeeLogger()


def extract_product_annotations(gff_file, qstart, qend):
    results = []
    with open(gff_file) as handle:
        for rec in GFF.parse(handle):
            for feat in rec.features:
                if feat.type == "CDS":
                    start, end = int(feat.location.start), int(feat.location.end)
                    if not (end < qstart or start > qend):  # overlapping range
                        strand = feat.strand
                        strand_symbol = "+" if strand == 1 else "-" if strand == -1 else "?"
                        gene = feat.qualifiers.get("gene", [""])[0]
                        product = feat.qualifiers.get("product", ["unknown"])[0]
                        results.append({
                            "CDS_start": start,
                            "CDS_end": end,
                            "strand": strand_symbol,
                            "gene": gene,
                            "product": product
                        })
    return results


def run_prokka(fasta_path, out_dir, cpus=1):
    out_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run([
        "prokka",
        "--outdir", str(out_dir),
        "--prefix", out_dir.name,
        "--cpus", str(cpus),
        "--force",
        "--kingdom", "Bacteria",
        "--gcode", "11",
        str(fasta_path)
    ], check=True)


def main():
    parser = argparse.ArgumentParser(description="Annotate aligned CDS using Prokka")
    parser.add_argument("--input-csv-dir", required=True, help="Directory of Species_link_Genetic_Exchange_Prediction_d*.csv")
    parser.add_argument("--output-dir", required=True, help="Directory to store annotations and results")
    parser.add_argument("--novel-fasta", required=True, help="FASTA file containing all novel host sequences")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for Prokka")
    args = parser.parse_args()

    input_csv_dir = Path(args.input_csv_dir)
    output_dir = Path(args.output_dir)
    novel_fasta = Path(args.novel_fasta)
    prokka_dir = output_dir / "prokka_cds"
    split_dir = prokka_dir / "split_novel_fasta"
    prokka_dir.mkdir(parents=True, exist_ok=True)
    split_dir.mkdir(exist_ok=True)

    # Step 1: Extract all qseqids needed
    qseqids = set()
    for csv_file in input_csv_dir.glob("Species_link_Genetic_Exchange_Prediction_d*.csv"):
        with open(csv_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                qseqids.add(row["qseqid"])

    # Step 2: Split novel.fasta by required qseqids
    all_seqs = SeqIO.to_dict(SeqIO.parse(novel_fasta, "fasta"))
    for qid in qseqids:
        if qid in all_seqs:
            out_fasta = split_dir / f"{qid}.fasta"
            SeqIO.write([all_seqs[qid]], out_fasta, "fasta")
        else:
            print(f"[WARNING] qseqid {qid} not found in novel.fasta")

    # Step 3: Annotate each qseqid once with Prokka
    def annotate(fasta_path):
        qid = fasta_path.stem
        out_path = prokka_dir / qid
        gff = out_path / f"{qid}.gff"
        if not gff.exists():
            print(f"[INFO] Annotating {qid} with Prokka...")
            run_prokka(fasta_path, out_path, cpus=1)
        else:
            print(f"[SKIP] {qid} already annotated")

    fasta_files = list(split_dir.glob("*.fasta"))
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        pool.map(annotate, fasta_files)

    # Step 4: Parse CDS from each d-value CSV
    for csv_file in input_csv_dir.glob("Species_link_Genetic_Exchange_Prediction_d*.csv"):
        d_val = csv_file.stem.split("_d")[-1]
        out_csv = prokka_dir / f"Aligned_CDS_products_d{d_val}.csv"
        annotated_rows = []

        with open(csv_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                qid = row["qseqid"]
                sid = row["sseqid"]
                try:
                    qstart = int(row["qstart"])
                    qend = int(row["qend"])
                except ValueError:
                    continue

                pair_num = row.get("pair_num", "")
                novel_host = row.get("novel host", "")
                known_host = row.get("known host", "")

                gff_file = prokka_dir / qid / f"{qid}.gff"
                if not gff_file.exists():
                    print(f"[WARNING] Missing GFF for {qid}, skipping.")
                    continue

                cds_matches = extract_product_annotations(gff_file, qstart, qend)
                for cds in cds_matches:
                    annotated_rows.append({
                        "pair_num": pair_num,
                        "qseqid": qid,
                        "sseqid": sid,
                        "qstart": qstart,
                        "qend": qend,
                        "CDS_start": cds["CDS_start"],
                        "CDS_end": cds["CDS_end"],
                        "strand": cds["strand"],
                        "gene": cds["gene"],
                        "product": cds["product"],
                        "novel host": novel_host,
                        "known host": known_host
                    })

        if annotated_rows:
            with open(out_csv, "w", newline="") as outf:
                writer = csv.DictWriter(outf, fieldnames=annotated_rows[0].keys())
                writer.writeheader()
                writer.writerows(annotated_rows)
            print(f"[✓] Wrote {len(annotated_rows)} CDS hits to {out_csv}")
        else:
            print(f"[!] No overlapping CDS found for d={d_val}")


if __name__ == "__main__":
    main()
