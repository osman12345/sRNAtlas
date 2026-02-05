#!/bin/bash
# =============================================================================
# sRNA-seq Analysis WebTool - Server Installation Script (No R Required)
# For Debian/Ubuntu servers
# =============================================================================

set -e  # Exit on error

echo "=========================================="
echo "sRNA-seq Analysis WebTool - Server Setup"
echo "=========================================="
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo -e "${RED}Conda not found!${NC}"
    echo ""
    echo "Please install Miniconda first:"
    echo "  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    echo "  bash Miniconda3-latest-Linux-x86_64.sh -b"
    echo "  ~/miniconda3/bin/conda init bash"
    echo "  source ~/.bashrc"
    echo ""
    exit 1
fi

echo -e "${GREEN}✓ Conda found${NC}"

# Environment name
ENV_NAME="srna-webtool"

# Check if environment already exists
if conda env list | grep -q "^${ENV_NAME} "; then
    echo -e "${YELLOW}Environment '${ENV_NAME}' already exists.${NC}"
    read -p "Do you want to remove and recreate it? (y/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing environment..."
        conda env remove -n ${ENV_NAME} -y
    else
        echo "Activating existing environment..."
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate ${ENV_NAME}
        echo -e "${GREEN}✓ Environment activated${NC}"
        exit 0
    fi
fi

echo ""
echo "Creating conda environment with Python 3.11..."
conda create -n ${ENV_NAME} python=3.11 -y

# Activate environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate ${ENV_NAME}

echo -e "${GREEN}✓ Environment created and activated${NC}"

# Upgrade pip
echo ""
echo "Upgrading pip..."
pip install --upgrade pip

# Install Python packages
echo ""
echo "Installing Python packages..."
pip install \
    streamlit>=1.28.0 \
    streamlit-option-menu>=0.3.6 \
    pandas>=2.0.0 \
    numpy>=1.24.0 \
    scipy>=1.10.0 \
    plotly>=5.15.0 \
    matplotlib>=3.7.0 \
    seaborn>=0.12.0 \
    statsmodels>=0.14.0 \
    scikit-learn>=1.3.0 \
    biopython>=1.81 \
    pydeseq2>=0.4.0 \
    anndata>=0.10.0 \
    gprofiler-official>=1.0.0 \
    openpyxl>=3.1.0 \
    xlsxwriter>=3.1.0 \
    tqdm>=4.65.0 \
    loguru>=0.7.0 \
    pyyaml>=6.0 \
    jinja2>=3.1.0

echo -e "${GREEN}✓ Python packages installed${NC}"

# Optional: Install bioinformatics tools for full pipeline
echo ""
read -p "Install Bowtie and Samtools for full pipeline? (y/n): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Installing bioinformatics tools..."
    conda install -c bioconda bowtie samtools pysam -y
    echo -e "${GREEN}✓ Bowtie and Samtools installed${NC}"
fi

# Verify installation
echo ""
echo "Verifying installation..."
python -c "import streamlit; print(f'  Streamlit: {streamlit.__version__}')"
python -c "import pydeseq2; print(f'  pyDESeq2: {pydeseq2.__version__}')"
python -c "import pandas; print(f'  Pandas: {pandas.__version__}')"
python -c "import plotly; print(f'  Plotly: {plotly.__version__}')"

echo ""
echo -e "${GREEN}=========================================="
echo "Installation Complete!"
echo "==========================================${NC}"
echo ""
echo "To use the tool:"
echo ""
echo "  1. Activate the environment:"
echo "     conda activate ${ENV_NAME}"
echo ""
echo "  2. Run the application:"
echo "     streamlit run app/main.py --server.port 8501 --server.headless true --server.address 0.0.0.0"
echo ""
echo "  3. Access in browser:"
echo "     http://YOUR_SERVER_IP:8501"
echo ""
echo "For production deployment, see DEBIAN_SERVER_DEPLOY.md"
echo ""
