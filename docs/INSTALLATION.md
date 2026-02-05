# sRNAtlas Installation Guide

Detailed setup instructions for all platforms.

---

## Quick Install

```bash
# Clone repository
git clone https://github.com/yourusername/sRNAtlas.git
cd sRNAtlas

# Create conda environment
conda create -n srnatlas python=3.10
conda activate srnatlas

# Install dependencies
pip install -r requirements.txt
conda install -c bioconda bowtie samtools

# Run
streamlit run app/main.py
```

---

## System Requirements

### Minimum
- Python 3.9+
- 4 GB RAM
- 10 GB disk space

### Recommended
- Python 3.10+
- 16 GB RAM
- 50 GB disk space
- SSD storage

---

## Platform-Specific Instructions

### Linux (Ubuntu/Debian)

```bash
# System dependencies
sudo apt-get update
sudo apt-get install -y build-essential zlib1g-dev libbz2-dev liblzma-dev

# Install Miniconda (if not installed)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Create environment
conda create -n srnatlas python=3.10
conda activate srnatlas

# Install bioinformatics tools
conda install -c bioconda bowtie samtools

# Install Python packages
pip install streamlit streamlit-option-menu pandas numpy plotly pysam cutadapt

# Clone and run
git clone https://github.com/yourusername/sRNAtlas.git
cd sRNAtlas
streamlit run app/main.py
```

### macOS

```bash
# Install Homebrew (if not installed)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install Miniconda
brew install --cask miniconda
conda init zsh  # or bash

# Create environment
conda create -n srnatlas python=3.10
conda activate srnatlas

# Install tools
conda install -c bioconda bowtie samtools
pip install cutadapt

# Install Python packages
pip install streamlit streamlit-option-menu pandas numpy plotly pysam

# Clone and run
git clone https://github.com/yourusername/sRNAtlas.git
cd sRNAtlas
streamlit run app/main.py
```

**Apple Silicon (M1/M2) Note:**
```bash
# If conda install fails, try:
CONDA_SUBDIR=osx-64 conda install -c bioconda bowtie samtools
```

### Windows (WSL2)

We recommend using Windows Subsystem for Linux:

```bash
# Enable WSL2 (PowerShell as Admin)
wsl --install

# Open Ubuntu terminal, then follow Linux instructions above
```

---

## Python Dependencies

### Core Requirements
```
streamlit>=1.28.0
streamlit-option-menu>=0.3.6
pandas>=2.0.0
numpy>=1.24.0
plotly>=5.18.0
pysam>=0.22.0
```

### Optional Dependencies
```
pydeseq2>=0.4.0      # Differential expression
scipy>=1.11.0        # Statistics
scikit-learn>=1.3.0  # Machine learning
matplotlib>=3.8.0    # Additional plotting
seaborn>=0.13.0      # Statistical visualization
```

---

## External Tools

### Bowtie (Required for Alignment)
```bash
# Conda (recommended)
conda install -c bioconda bowtie

# Verify
bowtie --version
```

**Note:** Use Bowtie v1, not Bowtie2. Bowtie v1 is optimized for short reads (<50 bp).

### Samtools (Required for BAM processing)
```bash
# Conda
conda install -c bioconda samtools

# Verify
samtools --version
```

### Cutadapt (Required for Trimming)
```bash
# pip (recommended)
pip install cutadapt

# Verify
cutadapt --version
```

---

## Running sRNAtlas

### Development Mode
```bash
cd sRNAtlas
streamlit run app/main.py
```

Access at: http://localhost:8501

### Production Mode
```bash
streamlit run app/main.py --server.port 80 --server.headless true
```

---

## Troubleshooting Installation

### "bowtie: command not found"
```bash
# Check if installed
which bowtie

# Or install via conda
conda install -c bioconda bowtie
```

### "No module named 'pysam'"
```bash
# Install with pip
pip install pysam

# Or use conda
conda install -c bioconda pysam
```

### Streamlit errors
```bash
# Clear cache
streamlit cache clear

# Reinstall streamlit
pip install --upgrade streamlit streamlit-option-menu
```

---

## Updating sRNAtlas

```bash
cd sRNAtlas
git pull origin main
pip install -r requirements.txt --upgrade
```

---

## Getting Help

- Check the Help module in the application
- Review documentation in docs/ folder
- Open an issue on GitHub

---

*sRNAtlas Installation Guide v1.3.0*
