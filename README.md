# GEM

**Version 1.1.2**

**Genetic Exchange Model (GEM)** is a cross-platform bioinformatics pipeline for analyzing genetic exchange between known and potential microbial hosts using comparative genomics. Version 1.1.2 improves host-link summarization and adds support for `--version` CLI display.

![Conda](https://img.shields.io/conda/vn/shihai_liu/gem?label=Install%20with%20conda)

---

## ✨ Features

* Filters input sequences by user-defined length
* Performs multi-threaded BLAST comparisons to detect homologous gene contexts
* Expands genetic regions upstream and downstream of target genes
* Infers genetic exchange by linking novel and known host sequences
* Automatically generates summary tables with predicted exchange events
* Tracks and summarizes host-host linkage pairs across multiple alignments
* Reports GEM version using `--version`

---

## 📦 Installation

### 🐧🍎 Option 1: Install with Conda (Recommended for Linux/macOS)

```bash
conda install -c shihai_liu -c conda-forge -c bioconda gem=1.1.2
```

### 🪟 Option 2: Windows Users

#### ✅ (A) Use WSL (Windows Subsystem for Linux)

```powershell
wsl --install
```
Then install GEM:

```bash
conda install -c shihai_liu -c conda-forge -c bioconda gem=1.1.2
```

#### ✅ (B) Native Windows (Advanced Users)

```bash
pip install git+https://github.com/LIUShi-Hai/GEM.git@v1.1.2
```

---

## 🚀 Usage

```bash
gem run-all --target target.fasta --known known.fasta --novel novel.fasta --email you@example.com --threads 4
```

```bash
gem --version
```

```bash
gem run-all --help
```

---

## 🔧 Key Parameters

* `--target`: Reference sequence (FASTA) of target gene (**required**)
* `--known`: FASTA file of known host sequences (**required**)
* `--novel`: FASTA of potential novel host sequences (**required**)
* `--email`: Your email (for NCBI Entrez)
* `--threads`: Number of BLAST threads (default: `1`)
* `--min-len`: Minimum sequence length (default: `5000`)
* `--segment-size`: Up/downstream context in bp (default: `5000`)
* `--d-range`: Expansion distances (default: `0 12000 2000`)
* `--coverage-threshold`: Total alignment length (default: `4000`)
* `--identity-threshold`: % identity threshold (default: `80.0`)
* `--evalue-threshold`: Max e-value (default: `1e-3`)

---

## ⏱ Run in Background with `nohup`

```bash
nohup yes | gem run-all ... > gem.log 2>&1 &
```

```bash
tail -f gem.log
```

---

## 🧪 Test Example

```bash
gem run-all --target test/target.fasta --known test/known.fasta --novel test/novel.fasta --email you@example.com --threads 2
```

---

## 🗂 Output Files

* `gem_output/`
  * `blast_query_subject_pair_counts.csv`
  * `Species_link_Genetic_Exchange_Prediction_d{d}.csv`
  * `host_link_summary_d{d}.csv`

---

## 🧩 Optional: Annotate Aligned CDS with Prokka

After running `gem run-all`, you can optionally annotate novel host genomes using [Prokka](https://github.com/tseemann/prokka) and extract coding sequences (CDS) overlapping the aligned regions.

### 📌 Purpose

After running `gem run-all`, use this tool to:

- Annotate novel genomes with [Prokka](https://github.com/tseemann/prokka)
- Identify overlapping CDS regions for gene context analysis
- Generate `Aligned_CDS_products_d{d}.csv` files for each expansion distance

### 🛠️ Requirements

Install dependencies in a new environment:

```bash
conda create -n gem_cds_env python=3.10 prokka biopython -c bioconda -c conda-forge
conda activate gem_cds_env
conda install -c bioconda -c conda-forge bcbio-gff
```

### ▶️ Usage

```bash
python annotate_aligned_cds.py --input-csv-dir ./gem_output --output-dir ./gem_output --novel-fasta ./test/novel.fasta --threads 4
```

**Arguments**:

* `--input-csv-dir`: Directory with `Species_link_Genetic_Exchange_Prediction_d*.csv` 
* `--output-dir`: Where results and annotations will be written
* `--novel-fasta`: FASTA file containing novel host genomes used in GEM
* `--threads`: Number of threads for Prokka

### 📂 Output

- `prokka_cds/` directory inside `--output-dir`
- Annotated `.gff` files in `prokka_cds/<qseqid>/`
- CDS product tables like:
  ```
  prokka_cds/Aligned_CDS_products_d0.csv
  prokka_cds/Aligned_CDS_products_d4000.csv
  ```

Each row includes:

* 'pair_num': Pair number from GEM
* 'qseqid': Novel genome contig ID
* 'sseqid': Known genome contig ID
* 'qstart/qend': Aligned region on query
* 'CDS_start/end': Coordinates of overlapping CDS
* 'strand': CDS strand (+ or -)
* 'gene': Gene name (if annotated)
* 'product': Functional description of CDS
* 'novel host': Species of the novel host
* 'known host': Species of the known host
  
### 📓 Logging

All runtime output is also saved to:
```bash
annotate_aligned_cds.log
```

### 💡 Notes

- This script is **optional** and intended for post-GEM analysis.
- It is safe to run using `nohup` for long jobs.
- It will automatically avoid re-annotating genomes that have already been processed.

---

## 📫 Contact

* 🧑‍🔬 Shihai Liu  
* 📧 [1330797686@qq.com](mailto:1330797686@qq.com)  
* 🔗 [GitHub: LIUShi-Hai/GEM](https://github.com/LIUShi-Hai/GEM)
