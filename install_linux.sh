#!/bin/bash
# =============================================================================
# sRNA-seq Analysis WebTool - Linux Installation Script
# For Ubuntu/Debian, CentOS/RHEL/Fedora
# No R required - uses Python-only packages
# =============================================================================

set -e  # Exit on error

echo "=========================================="
echo "sRNA-seq WebTool - Linux Installation"
echo "=========================================="
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Detect distribution
if [ -f /etc/os-release ]; then
    . /etc/os-release
    DISTRO=$ID
    VERSION=$VERSION_ID
else
    DISTRO="unknown"
fi

echo "Detected distribution: $DISTRO $VERSION"

# Install system dependencies
echo ""
echo "Installing system dependencies..."

if [[ "$DISTRO" == "ubuntu" || "$DISTRO" == "debian" ]]; then
    sudo apt-get update
    sudo apt-get install -y \
        python3 python3-pip python3-venv \
        build-essential \
        libhts-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
        libbz2-dev liblzma-dev zlib1g-dev

elif [[ "$DISTRO" == "fedora" ]]; then
    sudo dnf install -y \
        python3 python3-pip python3-devel \
        gcc gcc-c++ make \
        libcurl-devel openssl-devel libxml2-devel \
        bzip2-devel xz-devel zlib-devel

elif [[ "$DISTRO" == "centos" || "$DISTRO" == "rhel" || "$DISTRO" == "rocky" || "$DISTRO" == "almalinux" ]]; then
    sudo yum install -y epel-release || true
    sudo yum install -y \
        python3 python3-pip python3-devel \
        gcc gcc-c++ make \
        libcurl-devel openssl-devel libxml2-devel \
        bzip2-devel xz-devel zlib-devel

else
    echo -e "${YELLOW}Warning: Unknown distribution. Please install Python 3.9+ manually.${NC}"
fi

echo -e "${GREEN}✓ System dependencies installed${NC}"

# Create virtual environment
echo ""
echo "Creating Python virtual environment..."
python3 -m venv venv
source venv/bin/activate

# Upgrade pip
pip install --upgrade pip wheel setuptools

# Install Python packages
echo ""
echo "Installing Python packages..."

pip install \
    "streamlit>=1.28.0" \
    "streamlit-option-menu>=0.3.6" \
    "pandas>=2.0.0" \
    "numpy>=1.24.0" \
    "scipy>=1.10.0" \
    "plotly>=5.15.0" \
    "matplotlib>=3.7.0" \
    "seaborn>=0.12.0" \
    "statsmodels>=0.14.0" \
    "scikit-learn>=1.3.0" \
    "biopython>=1.81" \
    "pydeseq2>=0.4.0" \
    "anndata>=0.10.0" \
    "gprofiler-official>=1.0.0" \
    "openpyxl>=3.1.0" \
    "xlsxwriter>=3.1.0" \
    "tqdm>=4.65.0" \
    "loguru>=0.7.0" \
    "pyyaml>=6.0" \
    "jinja2>=3.1.0"

echo -e "${GREEN}✓ Core Python packages installed${NC}"

# Ask about full pipeline tools
echo ""
echo "==========================================="
echo "Optional: Full Pipeline Support"
echo "==========================================="
echo ""
echo "The full pipeline (FASTQ → Alignment → Counting) requires:"
echo "  - Bowtie (alignment - optimized for small RNA)"
echo "  - Samtools (BAM processing)"
echo "  - pysam (Python BAM interface)"
echo ""
echo "Without these, you can still use:"
echo "  - DE Analysis (from count matrix)"
echo "  - GO/KEGG Enrichment"
echo "  - All visualizations and reports"
echo ""

read -p "Install full pipeline tools? (y/n): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo ""
    echo "Installing bioinformatics tools..."

    if [[ "$DISTRO" == "ubuntu" || "$DISTRO" == "debian" ]]; then
        sudo apt-get install -y bowtie samtools
    elif [[ "$DISTRO" == "fedora" ]]; then
        sudo dnf install -y bowtie samtools
    else
        echo -e "${YELLOW}Installing via conda (apt/yum packages not available)...${NC}"
        if command -v conda &> /dev/null; then
            conda install -c bioconda bowtie samtools -y
        else
            echo -e "${RED}Conda not found. Please install Bowtie and Samtools manually.${NC}"
        fi
    fi

    echo "Installing pysam..."
    pip install pysam || echo -e "${YELLOW}Warning: pysam installation failed.${NC}"

    echo -e "${GREEN}✓ Full pipeline tools installed${NC}"
    FULL_PIPELINE=true
else
    echo -e "${YELLOW}Skipping full pipeline tools${NC}"
    FULL_PIPELINE=false
fi

# Verify installation
echo ""
echo "==========================================="
echo "Verifying Installation"
echo "==========================================="
echo ""

python -c "import streamlit; print(f'  ✓ Streamlit {streamlit.__version__}')"
python -c "import pydeseq2; print(f'  ✓ pyDESeq2 {pydeseq2.__version__}')"
python -c "import pandas; print(f'  ✓ Pandas {pandas.__version__}')"
python -c "import numpy; print(f'  ✓ NumPy {numpy.__version__}')"
python -c "import plotly; print(f'  ✓ Plotly {plotly.__version__}')"
python -c "import gprofiler; print('  ✓ g:Profiler (GO/KEGG analysis)')"

if [ "$FULL_PIPELINE" = true ]; then
    echo ""
    command -v bowtie &> /dev/null && echo "  ✓ Bowtie $(bowtie --version 2>&1 | head -1 | awk '{print $3}')"
    command -v samtools &> /dev/null && echo "  ✓ Samtools $(samtools --version | head -1 | awk '{print $2}')"
    python -c "import pysam; print(f'  ✓ pysam {pysam.__version__}')" 2>/dev/null || true
fi

# Success message
echo ""
echo -e "${GREEN}==========================================="
echo "Installation Complete!"
echo "==========================================${NC}"
echo ""
echo "To run the tool:"
echo ""
echo "  1. Activate the virtual environment:"
echo "     ${YELLOW}source venv/bin/activate${NC}"
echo ""
echo "  2. Start the tool:"
echo "     ${YELLOW}streamlit run app/main.py${NC}"
echo ""
echo "  3. Open http://localhost:8501 in your browser"
echo ""

if [ "$FULL_PIPELINE" = true ]; then
    echo "Full pipeline available: FASTQ → QC → Alignment → Counting → DE → Enrichment"
else
    echo "Available modules: DE Analysis, GO/KEGG Enrichment, Visualizations, Reports"
    echo ""
    echo "To add full pipeline later (Ubuntu/Debian):"
    echo "  sudo apt-get install bowtie samtools"
    echo "  source venv/bin/activate"
    echo "  pip install pysam"
fi
echo ""
