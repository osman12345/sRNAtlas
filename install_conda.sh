#!/bin/bash
# =============================================================================
# sRNA-seq Analysis WebTool - Conda Installation Script
# Works on macOS, Linux, and Windows (WSL/Git Bash)
# No R required - uses Python-only packages
# =============================================================================

set -e  # Exit on error

echo "=========================================="
echo "sRNA-seq Analysis WebTool - Conda Setup"
echo "=========================================="
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo -e "${RED}Error: Conda not found!${NC}"
    echo ""
    echo "Please install Miniconda first:"
    echo ""
    echo "  macOS (Apple Silicon):"
    echo "    curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh"
    echo "    bash Miniconda3-latest-MacOSX-arm64.sh"
    echo ""
    echo "  macOS (Intel):"
    echo "    curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
    echo "    bash Miniconda3-latest-MacOSX-x86_64.sh"
    echo ""
    echo "  Linux:"
    echo "    curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    echo "    bash Miniconda3-latest-Linux-x86_64.sh"
    echo ""
    echo "  Windows: Download from https://docs.conda.io/en/latest/miniconda.html"
    echo ""
    exit 1
fi

echo -e "${GREEN}✓ Conda found: $(conda --version)${NC}"

# Environment name
ENV_NAME="srna-webtool"

# Check if environment already exists
if conda env list | grep -q "^${ENV_NAME} "; then
    echo -e "${YELLOW}Environment '${ENV_NAME}' already exists.${NC}"
    read -p "Remove and recreate? (y/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing environment..."
        conda env remove -n ${ENV_NAME} -y
    else
        echo "Using existing environment."
        source $(conda info --base)/etc/profile.d/conda.sh 2>/dev/null || true
        conda activate ${ENV_NAME}
        echo -e "${GREEN}✓ Environment activated${NC}"
        echo ""
        echo "Run the app with:"
        echo "  streamlit run app/main.py"
        exit 0
    fi
fi

# Create environment
echo ""
echo "Creating conda environment with Python 3.11..."
conda create -n ${ENV_NAME} python=3.11 -y

# Activate environment
source $(conda info --base)/etc/profile.d/conda.sh 2>/dev/null || true
conda activate ${ENV_NAME}

echo -e "${GREEN}✓ Environment created and activated${NC}"

# Upgrade pip
echo ""
echo "Upgrading pip..."
pip install --upgrade pip setuptools wheel

# Install core Python packages
echo ""
echo "Installing Python packages (this may take a few minutes)..."
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

echo -e "${GREEN}✓ Core packages installed${NC}"

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
    echo "Installing bioinformatics tools via conda..."
    conda install -c bioconda -c conda-forge bowtie samtools pysam -y
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
    if command -v bowtie &> /dev/null; then
        echo "  ✓ Bowtie $(bowtie --version 2>&1 | head -1 | awk '{print $3}')"
    fi
    if command -v samtools &> /dev/null; then
        echo "  ✓ Samtools $(samtools --version | head -1 | awk '{print $2}')"
    fi
    python -c "import pysam; print(f'  ✓ pysam {pysam.__version__}')" 2>/dev/null || true
fi

# Success message
echo ""
echo -e "${GREEN}==========================================="
echo "Installation Complete!"
echo "==========================================${NC}"
echo ""
echo "To use the tool:"
echo ""
echo "  1. Activate the environment:"
echo "     ${YELLOW}conda activate ${ENV_NAME}${NC}"
echo ""
echo "  2. Run the application:"
echo "     ${YELLOW}streamlit run app/main.py${NC}"
echo ""
echo "  3. Open in browser:"
echo "     http://localhost:8501"
echo ""

if [ "$FULL_PIPELINE" = true ]; then
    echo "Full pipeline available: FASTQ → QC → Alignment → Counting → DE → Enrichment"
else
    echo "Available modules: DE Analysis, GO/KEGG Enrichment, Visualizations, Reports"
    echo ""
    echo "To add full pipeline later:"
    echo "  conda activate ${ENV_NAME}"
    echo "  conda install -c bioconda bowtie samtools pysam -y"
fi
echo ""
