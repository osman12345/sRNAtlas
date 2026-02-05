<p align="center">
  <img src="assets/logo.svg" alt="sRNAtlas Logo" width="120" height="120">
</p>

<h1 align="center">sRNAtlas</h1>

<p align="center">
  <strong>Comprehensive Small RNA-seq Analysis Platform</strong>
</p>

<p align="center">
  <a href="#features">Features</a> •
  <a href="#installation">Installation</a> •
  <a href="#quick-start">Quick Start</a> •
  <a href="#documentation">Documentation</a> •
  <a href="#citation">Citation</a>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/python-3.9+-blue.svg" alt="Python">
  <img src="https://img.shields.io/badge/streamlit-1.28+-red.svg" alt="Streamlit">
  <img src="https://img.shields.io/badge/license-MIT-green.svg" alt="License">
  <img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg" alt="Platform">
  <img src="https://img.shields.io/badge/version-BETA-orange.svg" alt="Version">
</p>

<p align="center">
  <img src="https://img.shields.io/badge/tests-51%20passed-brightgreen.svg" alt="Tests">
  <img src="https://img.shields.io/badge/coverage-85%25-yellowgreen.svg" alt="Coverage">
</p>

---

## Overview

**sRNAtlas** is a powerful, user-friendly application for comprehensive small RNA sequencing (sRNA-seq) data analysis. Built with Streamlit, it provides an intuitive interface for researchers to process raw sequencing data through quality control, alignment, quantification, differential expression analysis, and functional enrichment—all without requiring command-line expertise.

### Why sRNAtlas?

- 🎯 **Purpose-built for small RNA**: Optimized parameters for miRNA, siRNA, piRNA analysis
- 🖥️ **No coding required**: Intuitive web interface for all analysis steps
- 🔬 **Complete pipeline**: From raw FASTQ to publication-ready figures
- 📊 **Interactive visualizations**: Explore your data with dynamic plots
- 💾 **Project management**: Save, load, and share analysis sessions
- 🧪 **Reproducible**: Consistent results with version-controlled parameters

---

## Features

### 📋 Analysis Pipeline

```
┌─────────────┐    ┌─────────────┐    ┌─────────────┐    ┌─────────────┐
│   Raw       │───▶│   Quality   │───▶│   Adapter   │───▶│  Reference  │
│   FASTQ     │    │   Control   │    │   Trimming  │    │  Database   │
└─────────────┘    └─────────────┘    └─────────────┘    └─────────────┘
                                                                │
┌─────────────┐    ┌─────────────┐    ┌─────────────┐           │
│  Functional │◀───│     DE      │◀───│    Read     │◀──────────┘
│  Enrichment │    │  Analysis   │    │  Counting   │
└─────────────┘    └─────────────┘    └─────────────┘
```

### Core Modules

| Module | Description | Key Features |
|--------|-------------|--------------|
| **📁 Project** | Organize your analysis | Sample management, metadata, save/load projects |
| **📊 Quality Control** | Assess read quality | Size distribution, quality scores, contamination check |
| **✂️ Trimming** | Remove adapters | Cutadapt integration, preset adapters, length filtering |
| **🗄️ Databases** | Reference management | miRBase, RNAcentral, custom FASTA, index building |
| **🔗 Alignment** | Map reads | Bowtie optimization for small RNA, multi-mapper handling |
| **🔬 Post-Align QC** | Alignment quality | Mapping stats, strand bias, 5' nucleotide analysis |
| **📈 Counting** | Quantification | Feature counting, count matrix generation |
| **🧬 DE Analysis** | Differential expression | pyDESeq2, volcano plots, heatmaps, PCA |
| **🔍 Novel miRNA** | Discovery | Identify unannotated small RNAs |
| **🧫 isomiR** | Variant analysis | Detect miRNA isoforms and modifications |
| **🎯 Targets** | Target prediction | miRNA target gene identification |
| **🧬 GO/Pathway** | Enrichment | Gene Ontology, KEGG pathway analysis |
| **⚡ Batch** | Automation | Full pipeline batch processing |
| **📋 Reports** | Export | HTML reports, figure export |

### Supported RNA Types

| RNA Type | Size Range | Characteristics |
|----------|-----------|-----------------|
| **miRNA** | 18-25 nt | Gene expression regulators, 5' U bias |
| **siRNA** | 20-24 nt | RNAi pathway, perfect complementarity |
| **piRNA** | 24-32 nt | Transposon silencing, germline |
| **tRF/tsRNA** | 14-40 nt | tRNA-derived fragments, stress response |
| **rsRF** | 15-40 nt | rRNA-derived fragments |
| **snoRNA** | 60-300 nt | Small nucleolar RNA, rRNA modification |
| **snRNA** | 100-300 nt | Small nuclear RNA, splicing |
| **Y RNA** | 80-120 nt | DNA replication, quality control |

**Plant-specific:** tasiRNA, phasiRNA, natsiRNA, hc-siRNA

### Advanced Features

- **🔬 Novel miRNA Discovery**: Identify unannotated small RNAs from unaligned reads
- **🧫 isomiR Analysis**: Detect 5'/3' variants, SNPs, and non-templated additions
- **📊 Multi-group Comparison**: ANOVA for >2 conditions with pairwise comparisons
- **🔥 Cluster Analysis**: Hierarchical clustering with interactive heatmaps
- **📈 Interactive Plots**: Zoom, pan, and export publication-ready figures

---

## Installation

### Prerequisites

- Python 3.9 or higher
- 8+ GB RAM recommended
- External tools: Bowtie, Samtools, Cutadapt

### Option 1: Conda (Recommended)

```bash
# Clone the repository
git clone https://github.com/osman12345/sRNAtlas.git
cd sRNAtlas

# Create conda environment
conda create -n srnatlas python=3.10 -y
conda activate srnatlas

# Install bioinformatics tools
conda install -c bioconda bowtie samtools -y

# Install Python dependencies
pip install -r requirements.txt
```

### Option 2: pip + Manual Tools

```bash
# Clone and setup
git clone https://github.com/osman12345/sRNAtlas.git
cd sRNAtlas

# Create virtual environment
python -m venv venv
source venv/bin/activate  # Linux/macOS
# or: venv\Scripts\activate  # Windows

# Install dependencies
pip install -r requirements.txt

# Install external tools separately
# See docs/INSTALLATION.md for platform-specific instructions
```

### Option 3: Docker

```bash
# Build the image
docker build -t srnatlas .

# Run the container
docker run -p 8501:8501 -v $(pwd)/data:/app/data srnatlas

# Or use docker-compose
docker-compose up -d
```

### Verify Installation

```bash
# Check Python packages
python -c "import streamlit; import pandas; import pysam; print('OK')"

# Check external tools
bowtie --version
samtools --version
cutadapt --version

# Run tests
python -m pytest
```

---

## Quick Start

### 1. Launch the Application

```bash
cd sRNAtlas
streamlit run app/main.py
```

Open your browser to `http://localhost:8501`

### 2. Create a Project

1. Click **Project** in the sidebar
2. Enter project name and select organism
3. Upload your FASTQ files

### 3. Run the Pipeline

| Step | Module | Action |
|------|--------|--------|
| 1 | **Quality Control** | Assess raw read quality |
| 2 | **Trimming** | Remove adapters (select preset) |
| 3 | **Databases** | Download miRBase + build index |
| 4 | **Alignment** | Map reads to reference |
| 5 | **Post-Align QC** | Verify alignment quality |
| 6 | **Counting** | Generate count matrix |
| 7 | **DE Analysis** | Compare conditions |
| 8 | **GO/Pathway** | Functional enrichment |

### 4. Export Results

- Download count matrices, DE results, and figures
- Generate HTML reports
- Save project for future analysis

📖 See the [Quick Start Guide](docs/QUICK_START.md) for a detailed walkthrough.

---

## Documentation

| Document | Description |
|----------|-------------|
| [Quick Start Guide](docs/QUICK_START.md) | Get running in 10 minutes |
| [User Guide](docs/USER_GUIDE.md) | Complete documentation |
| [Installation Guide](docs/INSTALLATION.md) | Platform-specific setup |

### In-App Help

The application includes comprehensive documentation accessible via the **Help** module:
- Quick start tutorials
- Module reference
- File format specifications
- FAQ and troubleshooting

---

## System Requirements

### Hardware

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| RAM | 8 GB | 16+ GB |
| Storage | 20 GB | 100+ GB |
| CPU | 4 cores | 8+ cores |

### Software

| Dependency | Version | Purpose |
|------------|---------|---------|
| Python | 3.9+ | Runtime |
| Bowtie | 1.3+ | Alignment (not Bowtie2) |
| Samtools | 1.17+ | BAM processing |
| Cutadapt | 4.0+ | Adapter trimming |

### Python Packages

```
streamlit>=1.28.0
streamlit-option-menu>=0.3.6
pandas>=2.0.0
numpy>=1.24.0
plotly>=5.18.0
pysam>=0.22.0
biopython>=1.81
scipy>=1.11.0
scikit-learn>=1.3.0
pydeseq2>=0.4.0
```

---

## Input/Output Formats

### Input Files

**FASTQ** (raw reads):
```
@SEQ_ID
TAGCTTATCAGACTGATGTTGA
+
IIIIIIIIIIIIIIIIIIIII
```

**Count Matrix** (CSV):
```csv
,Sample1,Sample2,Sample3,Sample4
hsa-miR-21-5p,1500,1450,2800,2750
hsa-miR-155-5p,200,180,650,620
```

**Sample Metadata** (CSV):
```csv
sample,condition,batch
Sample1,control,1
Sample2,control,1
Sample3,treatment,2
Sample4,treatment,2
```

### Output Files

| Output | Format | Description |
|--------|--------|-------------|
| Trimmed reads | FASTQ.gz | Adapter-free, filtered reads |
| Alignments | BAM/BAI | Sorted, indexed alignments |
| Count matrix | CSV | Read counts per feature |
| DE results | CSV | log2FC, p-value, FDR |
| Figures | PNG/HTML | Interactive visualizations |
| Reports | HTML | Comprehensive analysis summary |

---

## Configuration

### Alignment Parameters (Bowtie)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-v` | 1 | Mismatches allowed (0-3) |
| `-k` | 10 | Report up to k alignments |
| `--best` | On | Report best alignments first |
| `--strata` | On | Only report best stratum |

### Trimming Parameters (Cutadapt)

| Parameter | Default | Description |
|-----------|---------|-------------|
| Min length | 18 nt | Discard shorter reads |
| Max length | 35 nt | Discard longer reads |
| Quality cutoff | 20 | Trim low-quality bases |
| Error rate | 0.1 | Adapter matching tolerance |

### DE Analysis Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| FDR threshold | 0.05 | Significance cutoff |
| log2FC threshold | 0.585 | ~1.5-fold change |
| Min count | 10 | Filter low-count features |

---

## Testing

sRNAtlas includes a comprehensive test suite:

```bash
# Run all tests
python -m pytest

# Run with verbose output
python -m pytest -v

# Run specific test file
python -m pytest tests/test_file_handlers.py

# Run with coverage
python -m pytest --cov=utils --cov=modules
```

### Test Coverage

| Module | Tests | Status |
|--------|-------|--------|
| File handlers | 12 | ✅ Passing |
| Progress tracker | 5 | ✅ Passing |
| Cluster analysis | 6 | ✅ Passing |
| Novel miRNA | 8 | ✅ Passing |
| isomiR | 19 | ✅ Passing |
| **Total** | **51** | **All Passing** |

---

## Project Structure

```
sRNAtlas/
├── app/
│   └── main.py              # Application entry point
├── modules/
│   ├── qc_module.py         # Quality control
│   ├── trimming_module.py   # Adapter trimming
│   ├── alignment_module.py  # Read alignment
│   ├── counting_module.py   # Read counting
│   ├── de_module.py         # Differential expression
│   ├── novel_mirna_module.py # Novel miRNA discovery
│   ├── isomir_module.py     # isomiR analysis
│   └── ...
├── utils/
│   ├── file_handlers.py     # File I/O utilities
│   ├── plotting.py          # Visualization functions
│   ├── cluster_analysis.py  # Clustering utilities
│   └── ...
├── config/
│   └── settings.py          # Configuration
├── tests/
│   ├── conftest.py          # Test fixtures
│   └── test_*.py            # Test files
├── docs/
│   ├── QUICK_START.md
│   ├── USER_GUIDE.md
│   └── INSTALLATION.md
├── assets/
│   └── logo.svg             # Application logo
├── requirements.txt
└── README.md
```

---

## Troubleshooting

### Common Issues

| Problem | Cause | Solution |
|---------|-------|----------|
| 0 reads after trimming | Wrong adapter | Try different preset or check library kit |
| 0% alignment rate | Wrong reference | Verify organism, rebuild index |
| Memory error | Too many samples | Process in smaller batches |
| `bowtie: not found` | Not in PATH | Install via conda or add to PATH |
| `pysam` import error | Missing htslib | `conda install -c bioconda pysam` |

### Getting Help

1. Check the in-app **Help** module
2. Review [User Guide](docs/USER_GUIDE.md)
3. Search existing [GitHub Issues](https://github.com/YOUR_USERNAME/sRNAtlas/issues)
4. Open a new issue with:
   - Error message
   - Steps to reproduce
   - System information

---

## Contributing

Contributions are welcome! Please read our contributing guidelines before submitting pull requests.

### Development Setup

```bash
# Clone and setup
git clone https://github.com/osman12345/sRNAtlas.git
cd sRNAtlas

# Create dev environment
conda create -n srnatlas-dev python=3.10 -y
conda activate srnatlas-dev

# Install dev dependencies
pip install -r requirements.txt
pip install pytest pytest-cov black flake8

# Run tests
python -m pytest -v
```

### Code Style

- Follow PEP 8 guidelines
- Use type hints where appropriate
- Write docstrings for functions and classes
- Add tests for new features

---

## Roadmap

- [x] Core analysis pipeline
- [x] Novel miRNA discovery
- [x] isomiR analysis
- [x] Multi-group comparison
- [x] Cluster analysis
- [x] Unit test framework
- [x] Docker containerization
- [ ] Windows standalone installer
- [ ] Cloud deployment 
- [ ] API endpoints
- [ ] Batch job scheduler
- [ ] Enhanced visualization suite

---

## Citation

If you use sRNAtlas in your research, please cite:

```bibtex
@software{sRNAtlas,
  author = {Ayman Osman},
  title = {sRNAtlas: A Comprehensive Platform for Small RNA-seq Analysis},
  year = {2026},
  url = {https://github.com/osman12345/sRNAtlas},
  version = {BETA}
}
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

- [Streamlit](https://streamlit.io/) - Web application framework
- [Bowtie](http://bowtie-bio.sourceforge.net/) - Short read aligner
- [miRBase](https://www.mirbase.org/) - miRNA database
- [RNAcentral](https://rnacentral.org/) - ncRNA database
- [pyDESeq2](https://github.com/owkin/PyDESeq2) - Differential expression

---

<p align="center">
  <strong>Built with ❤️ for the small RNA research community</strong>
</p>

<p align="center">
  <a href="#srnatlas">Back to top ⬆️</a>
</p>
