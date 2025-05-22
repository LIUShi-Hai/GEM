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

Installs all dependencies, including Biopython and NCBI BLAST+:

```bash
conda install -c shihai_liu -c conda-forge -c bioconda gem=1.1.2
```

### 🪟 Option 2: Windows Users

#### ✅ (A) Use WSL (Windows Subsystem for Linux)

1. Install WSL:

   ```powershell
   wsl --install
   ```
2. Open Ubuntu (or another Linux distro) and install GEM:

   ```bash
   conda install -c shihai_liu -c conda-forge -c bioconda gem=1.1.2
   ```

#### ✅ (B) Native Windows (Advanced Users)

1. Ensure system dependencies:

   * Python ≥ 3.10
   * Biopython
   * NCBI BLAST+ (`makeblastdb`, `blastn`) in PATH

2. Install from GitHub:

   ```bash
   pip install git+https://github.com/LIUShi-Hai/GEM.git@v1.1.2
   ```

---

## 🚀 Usage

Run the full GEM pipeline:

```bash
gem run-all --target target.fasta --known known.fasta --novel novel.fasta --email you@example.com --threads 4
```

Check the installed version:

```bash
gem --version
```

View all available parameters:

```bash
gem run-all --help
```

### 🔧 Key Parameters

* `--target`: Reference sequence (FASTA file) of the target gene (**required**)
* `--known`: FASTA file of known host sequences containing the target gene (**required**)
* `--novel`: FASTA file of potential novel host sequences lacking the target gene (**required**)
* `--email`: Your email address (required by NCBI Entrez)
* `--threads`: Number of threads for BLAST (default: `1`)
* `--min-len`: Minimum sequence length to retain (default: `5000`)
* `--segment-size`: Up/downstream extraction length in bp (default: `5000`)
* `--d-range`: Expansion distances: start end step (default: `0 12000 2000`)
* `--coverage-threshold`: Minimum alignment length for a pair (default: `4000`)
* `--identity-threshold`: Minimum % identity (default: `80.0`)
* `--evalue-threshold`: BLAST e-value cutoff (default: `1e-3`)

---

## ⏱ Run in Background with `nohup`

GEM may take time with large datasets. Run in background:

#### Option 1: Automatically confirm overwrite

```bash
nohup yes | gem run-all ... > gem.log 2>&1 &
```

#### Option 2: Manually delete previous output

```bash
rm -rf gem-output
nohup gem run-all ... > gem.log 2>&1 &
```

Monitor:

```bash
tail -f gem.log
```

---

## 🧪 Test Example

```bash
gem run-all --target test/target.fasta --known test/known.fasta --novel test/novel.fasta --email you@example.com --threads 2
```

### 🗂 Output Files

* `gem-output/`
  * `blast_query_subject_pair_counts.csv`: Number of valid query-subject pairs per `d`
  * `Species_link_Genetic_Exchange_Prediction_d{d}.csv`: Detailed alignments
  * `host_link_summary_d{d}.csv`: Summarized novel–known host linkages

---

## 📫 Contact

* 🧑‍🔬 Shihai Liu  
* 📧 [1330797686@qq.com](mailto:1330797686@qq.com)  
* 🔗 [GitHub: LIUShi-Hai/GEM](https://github.com/LIUShi-Hai/GEM)
