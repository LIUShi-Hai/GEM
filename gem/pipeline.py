import os
import csv
import subprocess
import shutil
from collections import defaultdict
from Bio import SeqIO, Entrez
from Bio.SeqRecord import SeqRecord

def filter_sequences(input_fasta, output_fasta, min_len):
    records = [rec for rec in SeqIO.parse(input_fasta, "fasta") if len(rec.seq) >= min_len]
    if not records:
        raise RuntimeError(f"No sequences ≥{min_len} bp in {input_fasta}")
    SeqIO.write(records, output_fasta, "fasta")
    print(f"Filtered {len(records)} sequences ≥{min_len} bp into {output_fasta}")

def build_db(ref_fasta, db_name):
    subprocess.run(["makeblastdb", "-in", ref_fasta, "-dbtype", "nucl", "-out", db_name], check=True)

def run_blast(query_fasta, db_name, output_file, evalue, threads=1):
    header = [
        "qseqid", "sseqid", "pident", "length",
        "mismatch", "gapopen", "qstart", "qend",
        "sstart", "send", "evalue", "bitscore"
    ]
    with open(output_file, 'w', newline='') as outf:
        outf.write(','.join(header) + '\n')
    with open(output_file, 'a', newline='') as outf:
        subprocess.run([
            "blastn", "-query", query_fasta, "-db", db_name,
            "-outfmt", "10", "-evalue", str(evalue),
            "-num_threads", str(threads)
        ], stdout=outf, check=True)

def parse_blast_results(blast_file, fasta_file):
    seq_lengths = {rec.id: len(rec.seq) for rec in SeqIO.parse(fasta_file, "fasta")}
    hits = []
    with open(blast_file, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            acc = row["qseqid"]
            if acc in seq_lengths:
                hits.append({
                    "accession": acc,
                    "target_gene_start": int(row["qstart"]),
                    "target_gene_end": int(row["qend"]),
                    "sequence_length": seq_lengths[acc]
                })
    return hits

def write_expanded_coords(hits, d, seg_size, output_file):
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            "accession", "target_gene_start", "target_gene_end", "sequence_length",
            "hypo_mob_seg_start", "hypo_mob_seg_end",
            "upstream_start", "upstream_end",
            "downstream_start", "downstream_end"
        ])
        writer.writeheader()
        for h in hits:
            tg_start, tg_end, seq_len = h["target_gene_start"], h["target_gene_end"], h["sequence_length"]
            seg_start = max(0, min(tg_start, tg_end) - d)
            seg_end = min(seq_len, max(tg_start, tg_end) + d)
            writer.writerow({
                **h,
                "hypo_mob_seg_start": seg_start,
                "hypo_mob_seg_end": seg_end,
                "upstream_start": max(0, seg_start - seg_size),
                "upstream_end": seg_start,
                "downstream_start": seg_end,
                "downstream_end": min(seq_len, seg_end + seg_size)
            })

def extract_segments(fasta_file, d_values, output_dir):
    seqs = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    for d in d_values:
        coord_file = os.path.join(output_dir, f"target_gene_coordinates_d{d}.csv")
        merged_file = os.path.join(output_dir, f"Genetic_contexts_merged_d{d}.fasta")
        if not os.path.exists(coord_file):
            continue
        merged = []
        with open(coord_file) as f:
            for row in csv.DictReader(f):
                acc = row["accession"]
                if acc not in seqs:
                    continue
                full = seqs[acc].seq
                up = full[int(row["upstream_start"]):int(row["upstream_end"])]
                dn = full[int(row["downstream_start"]):int(row["downstream_end"])]
                merged.append(SeqRecord(up + dn, id=f"{acc}_merged_d{d}", description=""))
        if merged:
            SeqIO.write(merged, merged_file, "fasta")
        else:
            print(f"[WARNING] No merged sequences written for d={d}; skipping this context.")

def extract_accession(s):
    return s.split(".")[0].split(":")[0]

def get_species(acc, cache):
    if acc in cache:
        return cache[acc]
    try:
        record = Entrez.read(Entrez.esummary(db="nuccore", id=acc))
        title = record[0]['Title']
        species = " ".join(title.split()[:2])
    except:
        species = "Unknown"
    cache[acc] = species
    return species

def link_and_filter(novel_fasta, d_values, thresholds, output_dir, threads=1):
    from collections import defaultdict
    import csv
    import os
    from Bio import SeqIO, Entrez

    summary = []

    for d in d_values:
        context = os.path.join(output_dir, f"Genetic_contexts_merged_d{d}.fasta")
        if not os.path.exists(context):
            print(f"[INFO] Skipping d={d}: {context} does not exist.")
            continue

        try:
            records = list(SeqIO.parse(context, "fasta"))
        except Exception as e:
            print(f"[ERROR] Could not parse {context}: {e}")
            continue

        if not records:
            print(f"[WARNING] No sequences found in {context}; skipping d={d}.")
            continue

        db_dir = os.path.join(output_dir, f"blast_results_d{d}")
        os.makedirs(db_dir, exist_ok=True)
        db = os.path.join(db_dir, "merged_context_db")
        subprocess.run(["makeblastdb", "-in", context, "-dbtype", "nucl", "-out", db], check=True)
        blast_out = os.path.join(db_dir, f"blast_result_d{d}.csv")
        subprocess.run([
            "blastn", "-query", novel_fasta, "-db", db,
            "-evalue", str(thresholds['evalue']),
            "-outfmt", "10 qseqid qstart qend sseqid sstart send pident length evalue",
            "-num_threads", str(threads),
            "-out", blast_out
        ], check=True)

        filtered_out = os.path.join(output_dir, f"Species_link_Genetic_Exchange_Prediction_d{d}.csv")
        top_hits = defaultdict(list)
        with open(blast_out) as f:
            for row in csv.reader(f):
                if len(row) != 9:
                    continue
                hit = dict(zip(["qseqid", "qstart", "qend", "sseqid", "sstart", "send", "pident", "length", "evalue"], row))
                top_hits[(hit['qseqid'], hit['sseqid'])].append(hit)

        cache, unique = {}, 0
        unique_host_links = set()

        with open(filtered_out, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=["pair_num"] + list(hit.keys()) + ["novel host", "known host"])
            writer.writeheader()
            for k, v in top_hits.items():
                v = sorted(v, key=lambda x: float(x['evalue']))

                if len(v) >= 2:
                    top2 = v[:2]
                    lengths = [int(h['length']) for h in top2]
                    if (sum(lengths) < thresholds['coverage']) or any(l < thresholds['coverage'] / 2 for l in lengths):
                        continue
                    v = top2
                elif len(v) == 1:
                    if int(v[0]['length']) < thresholds['coverage']:
                        continue
                else:
                    continue

                if float(v[0]['pident']) < thresholds['identity']:
                    continue
                if float(v[0]['evalue']) > thresholds['evalue']:
                    continue

                q, s = extract_accession(k[0]), extract_accession(k[1])
                q_sp, s_sp = get_species(q, cache), get_species(s, cache)
                unique += 1
                unique_host_links.add((q_sp, s_sp))

                for h in v:
                    h['pair_num'] = unique
                    h['novel host'] = q_sp
                    h['known host'] = s_sp
                    writer.writerow(h)

        # Write host_link_summary_d{d}.csv
        link_summary_file = os.path.join(output_dir, f"host_link_summary_d{d}.csv")
        with open(link_summary_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=["novel host", "known host"])
            writer.writeheader()
            for nh, kh in sorted(unique_host_links):
                writer.writerow({"novel host": nh, "known host": kh})

        summary.append({"d": d, "unique_query_subject_pairs": unique})

    with open(os.path.join(output_dir, "blast_query_subject_pair_counts.csv"), 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["d", "unique_query_subject_pairs"])
        writer.writeheader()
        writer.writerows(summary)


def run_all(
    target, known, novel, email,
    min_len=5000,
    segment_size=5000,
    d_range=(0, 12000, 2000),
    coverage_threshold=4000,
    identity_threshold=80.0,
    evalue_threshold=1e-3,
    output_dir="gem_output",
    force_overwrite=False,
    threads=1
):
    Entrez.email = email

    if os.path.exists(output_dir):
        if not force_overwrite:
            confirm = input(f"Output directory '{output_dir}' exists and will be removed. Continue? [y/N]: ")
            if confirm.lower() != 'y':
                print("Aborting.")
                return
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    filtered = os.path.join(output_dir, "known_filtered_5kb.fasta")
    blast_csv = os.path.join(output_dir, "target_gene_blast_result_10.csv")
    d_values = list(range(d_range[0], d_range[1] + 1, d_range[2]))
    thresholds = {
        "coverage": coverage_threshold,
        "identity": identity_threshold,
        "evalue": evalue_threshold
    }

    filter_sequences(known, filtered, min_len)
    build_db(target, os.path.join(output_dir, "target_gene_ref_db"))
    run_blast(filtered, os.path.join(output_dir, "target_gene_ref_db"), blast_csv, evalue_threshold, threads)

    hits = parse_blast_results(blast_csv, filtered)
    for d in d_values:
        coord_file = os.path.join(output_dir, f"target_gene_coordinates_d{d}.csv")
        write_expanded_coords(hits, d, segment_size, coord_file)

    extract_segments(filtered, d_values, output_dir)
    link_and_filter(novel, d_values, thresholds, output_dir, threads)

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="End-to-end pipeline for identifying target gene context and host linkage",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--target", required=True, help="FASTA file of target gene sequences")
    parser.add_argument("--known", required=True, help="FASTA file of known host sequences")
    parser.add_argument("--novel", required=True, help="FASTA file of novel host sequences to be linked")
    parser.add_argument("--email", required=True, help="Email for NCBI Entrez queries")
    parser.add_argument("--min-len", type=int, default=5000, help="Minimum length to keep known sequences")
    parser.add_argument("--segment-size", type=int, default=5000, help="Up/downstream window size in bp")
    parser.add_argument("--d-range", nargs=3, type=int, default=[0, 12000, 2000], help="Three space-separated integers: start, end, and step of expansion distances")
    parser.add_argument("--coverage-threshold", type=int, default=4000, help="Minimum total BLAST alignment length")
    parser.add_argument("--identity-threshold", type=float, default=80.0, help="Minimum BLAST identity percentage")
    parser.add_argument("--evalue-threshold", type=float, default=1e-3, help="Maximum acceptable BLAST e-value")
    parser.add_argument("--output-dir", default="gem_output", help="Directory to store all output files")
    parser.add_argument("--force", action="store_true", help="Force overwrite of output directory without confirmation")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads for BLAST")
    args = parser.parse_args()

    run_all(
        target=args.target,
        known=args.known,
        novel=args.novel,
        email=args.email,
        min_len=args.min_len,
        segment_size=args.segment_size,
        d_range=tuple(args.d_range),
        coverage_threshold=args.coverage_threshold,
        identity_threshold=args.identity_threshold,
        evalue_threshold=args.evalue_threshold,
        output_dir=args.output_dir,
        force_overwrite=args.force,
        threads=args.threads
    )

if __name__ == "__main__":
    main()
