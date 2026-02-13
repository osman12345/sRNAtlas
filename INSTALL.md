# Installation Guide for sRNA-seq Analysis WebTool

This guide provides step-by-step installation instructions for all platforms.

**No R required!** This tool uses Python-only packages (pyDESeq2, g:Profiler) for all analyses.

---

## Table of Contents

1. [Quick Start (Recommended)](#quick-start-recommended)
2. [macOS Installation](#macos-installation)
3. [Linux Installation](#linux-installation)
4. [Windows Installation](#windows-installation)
5. [Docker Installation](#docker-installation)
6. [Server Deployment](#server-deployment)
7. [Troubleshooting](#troubleshooting)

---

## Quick Start (Recommended)

The easiest way to install on any platform is using **Conda**:

```bash
# 1. Navigate to the tool directory
cd sRNA_WebTool

# 2. Run the installation script
chmod +x install_conda.sh
./install_conda.sh

# 3. Activate and run
conda activate srna-webtool
streamlit run app/main.py
```

Open http://localhost:8501 in your browser.

---

## macOS Installation

### Option 1: Using Conda (Recommended)

```bash
# Install Miniconda if not already installed
# Apple Silicon (M1/M2/M3):
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
bash Miniconda3-latest-MacOSX-arm64.sh

# Intel Mac:
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh

# Then run the install script
cd sRNA_WebTool
./install_conda.sh
```

### Option 2: Using Homebrew + venv

```bash
# Run the macOS install script
chmod +x install_macos.sh
./install_macos.sh
```

### Option 3: Manual Installation

```bash
# Install Python via Homebrew
brew install python@3.11

# Create virtual environment
python3 -m venv venv
source venv/bin/activate

# Install packages
pip install -r requirements.txt

# Run
streamlit run app/main.py
```

---

## Linux Installation

### Ubuntu/Debian

```bash
# Run the Linux install script
chmod +x install_linux.sh
./install_linux.sh
```

Or manually:

```bash
# Install system dependencies
sudo apt-get update
sudo apt-get install -y python3 python3-pip python3-venv

# Create virtual environment
python3 -m venv venv
source venv/bin/activate

# Install packages
pip install -r requirements.txt

# Run
streamlit run app/main.py
```

### CentOS/RHEL/Fedora

```bash
# Use the conda installation method
./install_conda.sh
```

### Using Conda (All Distributions)

```bash
# Install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Run install script
./install_conda.sh
```

---

## Windows Installation

### Option 1: Using Conda (Recommended)

1. Download and install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
2. Open "Anaconda Prompt" from Start Menu
3. Navigate to the tool directory:
   ```cmd
   cd path\to\sRNA_WebTool
   ```
4. Run:
   ```cmd
   conda create -n srna-webtool python=3.11 -y
   conda activate srna-webtool
   pip install -r requirements.txt
   streamlit run app/main.py
   ```

### Option 2: Using WSL2 (Full Linux Compatibility)

```powershell
# In PowerShell as Administrator
wsl --install -d Ubuntu
```

After restart, open Ubuntu and follow the [Linux Installation](#linux-installation) instructions.

---

## Docker Installation

Docker provides a consistent environment across all platforms.

### Prerequisites

- Install [Docker Desktop](https://www.docker.com/products/docker-desktop/)

### Quick Start

```bash
# Build and run with docker-compose
docker-compose up -d

# Or build manually
docker build -t srna-webtool .
docker run -p 8501:8501 -v $(pwd)/data:/app/data srna-webtool
```

Open http://localhost:8501 in your browser.

### With Data Volumes

```bash
# Create directories for data persistence
mkdir -p data results

# Run with mounted volumes
docker run -p 8501:8501 \
    -v $(pwd)/data:/app/data \
    -v $(pwd)/results:/app/results \
    srna-webtool
```

---

## Server Deployment

### Debian/Ubuntu Server

See [DEBIAN_SERVER_DEPLOY.md](DEBIAN_SERVER_DEPLOY.md) for detailed instructions including:
- Systemd service setup
- Nginx reverse proxy
- SSL/HTTPS with Let's Encrypt
- Production configuration

Quick start:

```bash
# SSH into your server
ssh username@your-server-ip

# Upload and install
scp -r sRNA_WebTool username@your-server-ip:~/
ssh username@your-server-ip
cd sRNA_WebTool
./install_server.sh

# Run
conda activate srna-webtool
streamlit run app/main.py --server.port 8501 --server.headless true --server.address 0.0.0.0
```

### Streamlit Cloud

See [STREAMLIT_CLOUD_DEPLOY.md](STREAMLIT_CLOUD_DEPLOY.md) for free cloud hosting.

---

## Core vs Full Pipeline

### Core Installation (Default)

Works everywhere, including Streamlit Cloud:
- ✅ Differential Expression (pyDESeq2)
- ✅ GO/KEGG Enrichment (g:Profiler)
- ✅ Interactive Visualizations
- ✅ Report Generation

### Full Pipeline (Optional)

Requires additional tools (local installation only):
- ✅ All core features
- ✅ FASTQ Quality Control
- ✅ Bowtie Alignment (optimized for small RNA)
- ✅ Read Counting from BAM

**Why Bowtie instead of Bowtie2?** Bowtie is preferred for small RNA-seq because:
- Designed for short reads (<50 bp)
- Exact mismatch control with `-v` mode
- No gap alignment (appropriate for miRNA/siRNA)
- Faster for very short sequences

To add full pipeline support:

```bash
# Conda (recommended)
conda install -c bioconda bowtie samtools pysam cutadapt miranda -y

# macOS (Homebrew)
brew install bowtie samtools
pip install pysam cutadapt

# Linux (apt)
sudo apt-get install bowtie samtools
pip install pysam cutadapt
```

### Optional: miRanda for Animal Target Prediction

miRanda is used for animal miRNA target prediction. If not installed, the tool falls back to local seed matching.

```bash
# Conda (recommended)
conda install -c bioconda miranda -y

# From source
wget http://cbio.mskcc.org/microrna_data/miRanda-aug2010.tar.gz
tar xzf miRanda-aug2010.tar.gz
cd miRanda-3.3a && ./configure && make && sudo make install
```

---

## Troubleshooting

### Common Issues

#### "ModuleNotFoundError: No module named 'streamlit'"

```bash
# Make sure environment is activated
conda activate srna-webtool
# or
source venv/bin/activate

# Reinstall
pip install streamlit
```

#### "pyDESeq2 installation fails"

```bash
pip install --upgrade pip setuptools wheel
pip install pydeseq2 --no-cache-dir
```

#### "pysam installation fails" (macOS)

```bash
brew install htslib
export HTSLIB_LIBRARY_DIR=/opt/homebrew/lib  # Apple Silicon
export HTSLIB_INCLUDE_DIR=/opt/homebrew/include
pip install pysam
```

#### "Streamlit won't start / Port in use"

```bash
# Check what's using the port
lsof -i :8501  # macOS/Linux
netstat -ano | findstr 8501  # Windows

# Use a different port
streamlit run app/main.py --server.port 8502
```

#### Memory issues during analysis

- Reduce number of threads in Settings
- Process fewer samples at a time
- Use a machine with more RAM
- Consider using cloud computing

---

## Verifying Installation

Run this to verify all packages are installed:

```bash
python -c "
import streamlit; print(f'✓ Streamlit: {streamlit.__version__}')
import pandas; print(f'✓ Pandas: {pandas.__version__}')
import numpy; print(f'✓ NumPy: {numpy.__version__}')
import scipy; print(f'✓ SciPy: {scipy.__version__}')
import plotly; print(f'✓ Plotly: {plotly.__version__}')
import sklearn; print(f'✓ scikit-learn: {sklearn.__version__}')
import Bio; print(f'✓ Biopython: {Bio.__version__}')
import pysam; print(f'✓ pysam: {pysam.__version__}')
import cutadapt; print(f'✓ cutadapt: {cutadapt.__version__}')
import pydeseq2; print(f'✓ pyDESeq2: {pydeseq2.__version__}')
import gprofiler; print('✓ g:Profiler: OK')
import loguru; print('✓ Loguru: OK')
import yaml; print('✓ PyYAML: OK')
try:
    import kaleido; print('✓ Kaleido: OK')
except: print('⚠ Kaleido: not installed (optional for image export)')
print('\\n✅ All packages installed!')
"
```

### Verify External Tools (Optional)

For full pipeline functionality with FASTQ processing:

```bash
# Check alignment tools
bowtie --version
samtools --version

# Check adapter trimming
cutadapt --version

# Check miRanda (optional, for animal target prediction)
miranda -h 2>&1 | head -1

# Check pysam (Python BAM handling)
python -c "import pysam; print(f'pysam: {pysam.__version__}')"
```

---

## Quick Reference

| Platform | Recommended Method | Command |
|----------|-------------------|---------|
| macOS | Conda | `./install_conda.sh` |
| Linux | Conda | `./install_conda.sh` |
| Windows | Conda | See Windows section |
| Docker | Docker Compose | `docker-compose up -d` |
| Server | install_server.sh | `./install_server.sh` |
| Cloud | Streamlit Cloud | See deployment guide |

---

## Getting Help

If you encounter issues not covered here:

1. Check the error message carefully
2. Search the error online
3. Check GitHub Issues (if available)
4. Try using Docker for a clean environment

---

## System Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| Python | 3.9 | 3.11 |
| RAM | 4GB | 8GB+ |
| Storage | 5GB | 20GB+ |
| CPU | 2 cores | 4+ cores |
