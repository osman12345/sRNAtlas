#!/bin/bash
# =============================================================================
# sRNA-seq Analysis WebTool - macOS Installation Script
# No R required - uses Python-only packages
# =============================================================================

set -e  # Exit on error

echo "=========================================="
echo "sRNA-seq WebTool - macOS Installation"
echo "=========================================="
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Detect architecture
ARCH=$(uname -m)
echo "Detected architecture: $ARCH"

# Check for Homebrew
if ! command -v brew &> /dev/null; then
    echo ""
    echo "Homebrew not found. Installing..."
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

    # Add to PATH for Apple Silicon
    if [[ "$ARCH" == "arm64" ]]; then
        echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.zprofile
        eval "$(/opt/homebrew/bin/brew shellenv)"
    fi
fi

echo -e "${GREEN}✓ Homebrew available${NC}"

# Install Python if needed
if ! command -v python3 &> /dev/null; then
    echo "Installing Python 3.11..."
    brew install python@3.11
fi
echo -e "${GREEN}✓ Python available: $(python3 --version)${NC}"

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

# Core packages
pip install \
    "streamlit>=1.28.0" \
    "streamlit-option-menu>=0.3.6" \
    "pandas>=2.0.0" \
    "numpy>=1.24.0" \
    "scipy>=1.10.0" \
    "plotly>=5.15.0" \
    "kaleido==0.2.1" \
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
    echo "Installing bioinformatics tools via Homebrew..."
    brew install bowtie samtools htslib 2>/dev/null || true

    echo "Installing pysam..."
    if [[ "$ARCH" == "arm64" ]]; then
        export HTSLIB_LIBRARY_DIR=/opt/homebrew/lib
        export HTSLIB_INCLUDE_DIR=/opt/homebrew/include
    else
        export HTSLIB_LIBRARY_DIR=/usr/local/lib
        export HTSLIB_INCLUDE_DIR=/usr/local/include
    fi
    pip install pysam || echo -e "${YELLOW}Warning: pysam installation failed. BAM processing will be limited.${NC}"

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
    command -v bowtie &> /dev/null && echo "  ✓ Bowtie installed"
    command -v samtools &> /dev/null && echo "  ✓ Samtools installed"
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
    echo "To add full pipeline later:"
    echo "  brew install bowtie samtools htslib"
    echo "  source venv/bin/activate"
    echo "  pip install pysam"
fi
echo ""
