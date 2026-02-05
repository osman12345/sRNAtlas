"""
Adapter Trimming Module for sRNAtlas
Integrates Cutadapt for adapter removal and quality trimming
"""
import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import subprocess
import tempfile
import shutil
import gzip
import re

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from config import config


# Common adapter sequences for small RNA libraries organized by platform
ADAPTER_PRESETS = {
    # ============ ILLUMINA PLATFORM - SMALL RNA ============
    "Illumina TruSeq Small RNA (RA3)": {
        "3_adapter": "TGGAATTCTCGGGTGCCAAGG",
        "5_adapter": None,
        "platform": "Illumina",
        "description": "TruSeq Small RNA 3' adapter RA3 (most common for small RNA)"
    },
    "Illumina TruSeq Small RNA (RA3+RA5)": {
        "3_adapter": "TGGAATTCTCGGGTGCCAAGG",
        "5_adapter": "GTTCAGAGTTCTACAGTCCGACGATC",
        "platform": "Illumina",
        "description": "TruSeq Small RNA with both 3' (RA3) and 5' (RA5) adapters"
    },
    "Illumina Small RNA v1.5": {
        "3_adapter": "ATCTCGTATGCCGTCTTCTGCTTG",
        "5_adapter": None,
        "platform": "Illumina",
        "description": "Older Illumina Small RNA Kit v1.5"
    },

    # ============ ILLUMINA PLATFORM - UNIVERSAL/TRUSEQ ============
    "Illumina Universal": {
        "3_adapter": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        "5_adapter": None,
        "platform": "Illumina",
        "description": "Illumina Universal adapter (full sequence)"
    },
    "Illumina TruSeq DNA/RNA Read 1": {
        "3_adapter": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        "5_adapter": None,
        "platform": "Illumina",
        "description": "TruSeq DNA/RNA Read 1 adapter"
    },
    "Illumina TruSeq DNA/RNA Read 2": {
        "3_adapter": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        "5_adapter": None,
        "platform": "Illumina",
        "description": "TruSeq DNA/RNA Read 2 adapter"
    },
    "Illumina TruSeq Index (i7)": {
        "3_adapter": "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        "5_adapter": None,
        "platform": "Illumina",
        "description": "TruSeq Index adapter with i7 barcode region"
    },
    "Illumina TruSeq Index (i5)": {
        "3_adapter": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        "5_adapter": None,
        "platform": "Illumina",
        "description": "TruSeq Index adapter with i5 barcode region"
    },

    # ============ ILLUMINA PLATFORM - NEXTERA ============
    "Illumina Nextera": {
        "3_adapter": "CTGTCTCTTATACACATCT",
        "5_adapter": None,
        "platform": "Illumina",
        "description": "Nextera Transposase adapter"
    },
    "Illumina Nextera Read 1": {
        "3_adapter": "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
        "5_adapter": None,
        "platform": "Illumina",
        "description": "Nextera DNA Read 1 adapter"
    },
    "Illumina Nextera Read 2": {
        "3_adapter": "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
        "5_adapter": None,
        "platform": "Illumina",
        "description": "Nextera DNA Read 2 adapter"
    },

    # ============ ILLUMINA PLATFORM - THIRD PARTY KITS ============
    "NEBNext Small RNA": {
        "3_adapter": "AGATCGGAAGAGCACACGTCT",
        "5_adapter": None,
        "platform": "Illumina",
        "description": "NEBNext small RNA library prep for Illumina"
    },
    "NEBNext Multiplex Small RNA": {
        "3_adapter": "AGATCGGAAGAGCACACGTCT",
        "5_adapter": "GTTCAGAGTTCTACAGTCCGACGATC",
        "platform": "Illumina",
        "description": "NEBNext Multiplex small RNA with 5' adapter"
    },
    "QIAseq miRNA": {
        "3_adapter": "AACTGTAGGCACCATCAAT",
        "5_adapter": None,
        "platform": "Illumina",
        "description": "QIAseq miRNA library prep for Illumina"
    },
    "QIAseq miRNA UDI": {
        "3_adapter": "AACTGTAGGCACCATCAATC",
        "5_adapter": None,
        "platform": "Illumina",
        "description": "QIAseq miRNA with Unique Dual Indexes"
    },
    "Lexogen Small RNA": {
        "3_adapter": "TGGAATTCTCGGGTGCCAAGGAACTCC",
        "5_adapter": None,
        "platform": "Illumina",
        "description": "Lexogen Small RNA-Seq Library Prep"
    },
    "Clontech/Takara SMARTer smRNA": {
        "3_adapter": "AAAAAAAAAA",
        "5_adapter": None,
        "platform": "Illumina",
        "description": "Takara SMARTer smRNA-Seq (uses poly-A tail)"
    },
    "Bioo Scientific NEXTflex": {
        "3_adapter": "TGGAATTCTCGGGTGCCAAGG",
        "5_adapter": None,
        "platform": "Illumina",
        "description": "Bioo Scientific NEXTflex Small RNA-Seq v3"
    },
    "IDT for Illumina": {
        "3_adapter": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        "5_adapter": None,
        "platform": "Illumina",
        "description": "IDT for Illumina DNA/RNA UD Indexes"
    },
    "PerkinElmer NEXTFLEX": {
        "3_adapter": "TGGAATTCTCGGGTGCCAAGG",
        "5_adapter": "GTTCAGAGTTCTACAGTCCGACGATC",
        "platform": "Illumina",
        "description": "PerkinElmer NEXTFLEX Small RNA-Seq Kit v4"
    },
    "Diagenode CATS Small RNA": {
        "3_adapter": "AAAAAAAAAA",
        "5_adapter": None,
        "platform": "Illumina",
        "description": "Diagenode CATS Small RNA-Seq Kit (poly-A)"
    },

    # ============ ION TORRENT PLATFORM ============
    "Ion Torrent Small RNA": {
        "3_adapter": "ATCACCGACTGCCCATAGAGAGG",
        "5_adapter": None,
        "platform": "Ion Torrent",
        "description": "Ion Torrent small RNA library (Ion S5/Proton/Genexus)"
    },
    "Ion Torrent A-adapter": {
        "3_adapter": "CCATCTCATCCCTGCGTGTCTCCGAC",
        "5_adapter": None,
        "platform": "Ion Torrent",
        "description": "Ion Torrent A-adapter (alternative)"
    },
    "Ion Torrent P1-adapter": {
        "3_adapter": "CCTCTCTATGGGCAGTCGGTGAT",
        "5_adapter": None,
        "platform": "Ion Torrent",
        "description": "Ion Torrent P1-adapter"
    },

    # ============ BGI/MGI PLATFORM ============
    "BGI/MGI Small RNA": {
        "3_adapter": "AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA",
        "5_adapter": None,
        "platform": "BGI/MGI",
        "description": "BGI DNBSEQ/MGI small RNA adapter"
    },
    "BGI/MGI Universal": {
        "3_adapter": "AAGTCGGAGGCCAAGCGGTC",
        "5_adapter": None,
        "platform": "BGI/MGI",
        "description": "BGI/MGI universal adapter (short version)"
    },
    "MGIEasy Small RNA": {
        "3_adapter": "AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA",
        "5_adapter": "GAACGACATGGCTACGATCCGACTT",
        "platform": "BGI/MGI",
        "description": "MGIEasy Small RNA Library Prep"
    },

    # ============ OXFORD NANOPORE ============
    "Nanopore cDNA-PCR": {
        "3_adapter": "ACTTGCCTGTCGCTCTATCTTC",
        "5_adapter": "TTTCTGTTGGTGCTGATATTGC",
        "platform": "Nanopore",
        "description": "Oxford Nanopore cDNA-PCR sequencing kit"
    },
    "Nanopore Direct RNA": {
        "3_adapter": "GGCGTCTGCTTGGGTGTTTAACC",
        "5_adapter": None,
        "platform": "Nanopore",
        "description": "Oxford Nanopore Direct RNA (usually no trimming needed)"
    },

    # ============ PACBIO ============
    "PacBio SMRTbell": {
        "3_adapter": "ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT",
        "5_adapter": None,
        "platform": "PacBio",
        "description": "PacBio SMRTbell adapter (HiFi/Sequel)"
    },

    # ============ LEGACY PLATFORMS ============
    "454/Roche A-adapter": {
        "3_adapter": "GCCTCCCTCGCGCCATCAG",
        "5_adapter": None,
        "platform": "454",
        "description": "454/Roche A-adapter (legacy)"
    },
    "454/Roche B-adapter": {
        "3_adapter": "GCCTTGCCAGCCCGCTCAG",
        "5_adapter": None,
        "platform": "454",
        "description": "454/Roche B-adapter (legacy)"
    },
    "SOLiD Small RNA": {
        "3_adapter": "CGCCTTGGCCGTACAGCAG",
        "5_adapter": None,
        "platform": "SOLiD",
        "description": "ABI SOLiD small RNA adapter (legacy)"
    },

    # ============ SPECIAL ============
    "Poly-A Tail": {
        "3_adapter": "AAAAAAAAAAAAAAAAAAAAAAAAA",
        "5_adapter": None,
        "platform": "Special",
        "description": "Poly-A tail trimming (25 A's)"
    },
    "Poly-T Head": {
        "3_adapter": None,
        "5_adapter": "TTTTTTTTTTTTTTTTTTTTTTTTT",
        "platform": "Special",
        "description": "Poly-T head trimming (25 T's)"
    },

    # ============ CUSTOM ============
    "Custom": {
        "3_adapter": "",
        "5_adapter": None,
        "platform": "Custom",
        "description": "Enter your own adapter sequence"
    }
}

# Group adapters by platform for easier selection
ADAPTER_PLATFORMS = {
    "Illumina - Small RNA": [
        "Illumina TruSeq Small RNA (RA3)",
        "Illumina TruSeq Small RNA (RA3+RA5)",
        "Illumina Small RNA v1.5",
        "NEBNext Small RNA",
        "NEBNext Multiplex Small RNA",
        "QIAseq miRNA",
        "QIAseq miRNA UDI",
        "Lexogen Small RNA",
        "Clontech/Takara SMARTer smRNA",
        "Bioo Scientific NEXTflex",
        "PerkinElmer NEXTFLEX",
        "Diagenode CATS Small RNA"
    ],
    "Illumina - Universal/TruSeq": [
        "Illumina Universal",
        "Illumina TruSeq DNA/RNA Read 1",
        "Illumina TruSeq DNA/RNA Read 2",
        "Illumina TruSeq Index (i7)",
        "Illumina TruSeq Index (i5)",
        "IDT for Illumina"
    ],
    "Illumina - Nextera": [
        "Illumina Nextera",
        "Illumina Nextera Read 1",
        "Illumina Nextera Read 2"
    ],
    "Ion Torrent": [
        "Ion Torrent Small RNA",
        "Ion Torrent A-adapter",
        "Ion Torrent P1-adapter"
    ],
    "BGI/MGI": [
        "BGI/MGI Small RNA",
        "BGI/MGI Universal",
        "MGIEasy Small RNA"
    ],
    "Nanopore": [
        "Nanopore cDNA-PCR",
        "Nanopore Direct RNA"
    ],
    "PacBio": [
        "PacBio SMRTbell"
    ],
    "Legacy": [
        "454/Roche A-adapter",
        "454/Roche B-adapter",
        "SOLiD Small RNA"
    ],
    "Special": [
        "Poly-A Tail",
        "Poly-T Head"
    ],
    "Custom": [
        "Custom"
    ]
}


def check_cutadapt() -> Tuple[bool, str]:
    """Check if cutadapt is installed and get version"""
    try:
        result = subprocess.run(
            ['cutadapt', '--version'],
            capture_output=True,
            text=True
        )
        if result.returncode == 0:
            version = result.stdout.strip()
            return True, version
        return False, "Not found"
    except FileNotFoundError:
        return False, "Not installed"


def run_cutadapt(
    input_file: Path,
    output_file: Path,
    adapter_3: str,
    adapter_5: Optional[str] = None,
    min_length: int = 18,
    max_length: int = 35,
    quality_cutoff: int = 20,
    trim_n: bool = True,
    discard_untrimmed: bool = False,
    cores: int = 4,
    error_rate: float = 0.1
) -> Dict:
    """
    Run cutadapt on a single FASTQ file

    Args:
        input_file: Input FASTQ file path
        output_file: Output trimmed FASTQ file path
        adapter_3: 3' adapter sequence
        adapter_5: Optional 5' adapter sequence
        min_length: Minimum read length after trimming
        max_length: Maximum read length after trimming
        quality_cutoff: Quality score cutoff for trimming
        trim_n: Whether to trim N bases from ends
        discard_untrimmed: Whether to discard reads without adapters
        cores: Number of CPU cores to use
        error_rate: Maximum error rate for adapter matching

    Returns:
        Dictionary with trimming statistics
    """
    # Build cutadapt command
    cmd = [
        'cutadapt',
        '-a', adapter_3,  # 3' adapter
        '-m', str(min_length),  # Minimum length
        '-M', str(max_length),  # Maximum length
        '-q', str(quality_cutoff),  # Quality cutoff
        '-e', str(error_rate),  # Error rate
        '-j', str(cores),  # Cores
        '-o', str(output_file),  # Output
    ]

    # Optional parameters
    if adapter_5:
        cmd.extend(['-g', adapter_5])  # 5' adapter

    if trim_n:
        cmd.append('--trim-n')

    if discard_untrimmed:
        cmd.append('--discard-untrimmed')

    # Add input file
    cmd.append(str(input_file))

    try:
        # Verify input file exists and has content
        if not input_file.exists():
            return {
                'status': 'error',
                'error': f'Input file not found: {input_file}'
            }

        file_size = input_file.stat().st_size
        if file_size == 0:
            return {
                'status': 'error',
                'error': f'Input file is empty: {input_file}'
            }

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600  # 1 hour timeout
        )

        if result.returncode != 0:
            return {
                'status': 'error',
                'error': result.stderr or result.stdout or 'Unknown error'
            }

        # Parse cutadapt output (cutadapt reports to stderr)
        # Try stderr first, then stdout (some versions differ)
        output_text = result.stderr if result.stderr else result.stdout
        stats = parse_cutadapt_output(output_text)

        # If no reads parsed, store raw output for debugging
        if stats['total_reads'] == 0:
            stats['debug_stderr'] = result.stderr[:1000] if result.stderr else 'No stderr output'
            stats['debug_stdout'] = result.stdout[:1000] if result.stdout else 'No stdout output'
            stats['input_file_size'] = file_size
            stats['debug_cmd'] = ' '.join(cmd)

        stats['status'] = 'success'
        stats['output_file'] = str(output_file)

        return stats

    except subprocess.TimeoutExpired:
        return {
            'status': 'error',
            'error': 'Trimming timed out (>1 hour)'
        }
    except Exception as e:
        return {
            'status': 'error',
            'error': str(e)
        }


def parse_cutadapt_output(stderr: str) -> Dict:
    """Parse cutadapt stderr output for statistics"""
    stats = {
        'total_reads': 0,
        'reads_with_adapter': 0,
        'reads_written': 0,
        'reads_too_short': 0,
        'reads_too_long': 0,
        'total_bp_input': 0,
        'total_bp_output': 0,
        'adapter_trimmed_bp': 0,
        'quality_trimmed_bp': 0
    }

    # Parse statistics from cutadapt output
    patterns = {
        'total_reads': r'Total reads processed:\s+([\d,]+)',
        'reads_with_adapter': r'Reads with adapters:\s+([\d,]+)',
        'reads_written': r'Reads written \(passing filters\):\s+([\d,]+)',
        'reads_too_short': r'Reads that were too short:\s+([\d,]+)',
        'reads_too_long': r'Reads that were too long:\s+([\d,]+)',
        'total_bp_input': r'Total basepairs processed:\s+([\d,]+)',
        'total_bp_output': r'Total written \(filtered\):\s+([\d,]+)',
    }

    for key, pattern in patterns.items():
        match = re.search(pattern, stderr)
        if match:
            stats[key] = int(match.group(1).replace(',', ''))

    # Calculate percentages
    if stats['total_reads'] > 0:
        stats['adapter_rate'] = 100 * stats['reads_with_adapter'] / stats['total_reads']
        stats['pass_rate'] = 100 * stats['reads_written'] / stats['total_reads']
    else:
        stats['adapter_rate'] = 0
        stats['pass_rate'] = 0

    return stats


def batch_trim_files(
    input_files: List[Path],
    output_dir: Path,
    adapter_3: str,
    adapter_5: Optional[str] = None,
    min_length: int = 18,
    max_length: int = 35,
    quality_cutoff: int = 20,
    cores: int = 4,
    progress_callback=None
) -> Tuple[List[Path], pd.DataFrame]:
    """
    Batch trim multiple FASTQ files

    Args:
        input_files: List of input FASTQ file paths
        output_dir: Output directory for trimmed files
        adapter_3: 3' adapter sequence
        adapter_5: Optional 5' adapter sequence
        min_length: Minimum read length
        max_length: Maximum read length
        quality_cutoff: Quality cutoff
        cores: Number of CPU cores
        progress_callback: Optional progress callback function

    Returns:
        Tuple of (output file paths, statistics dataframe)
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    output_files = []
    all_stats = []

    for i, input_file in enumerate(input_files):
        if progress_callback:
            progress_callback(i / len(input_files), f"Trimming {input_file.name}...")

        # Generate output filename
        stem = input_file.stem
        if stem.endswith('.fastq'):
            stem = stem[:-6]
        elif stem.endswith('.fq'):
            stem = stem[:-3]

        output_file = output_dir / f"{stem}_trimmed.fastq.gz"

        # Run cutadapt
        stats = run_cutadapt(
            input_file=input_file,
            output_file=output_file,
            adapter_3=adapter_3,
            adapter_5=adapter_5,
            min_length=min_length,
            max_length=max_length,
            quality_cutoff=quality_cutoff,
            cores=cores
        )

        stats['sample'] = input_file.stem

        # Check if output file was created and has content
        if stats['status'] == 'success':
            if output_file.exists():
                output_size = output_file.stat().st_size
                stats['output_size'] = output_size
                if output_size > 0:
                    output_files.append(output_file)
                else:
                    stats['warning'] = 'Output file is empty'
            else:
                stats['warning'] = 'Output file not created'

        all_stats.append(stats)

    # Create statistics dataframe
    stats_df = pd.DataFrame(all_stats)

    return output_files, stats_df


def render_trimming_page():
    """Render the adapter trimming page"""
    st.header("✂️ Adapter Trimming")

    st.markdown("""
    Remove adapters and low-quality bases from raw FASTQ files using **Cutadapt**.
    This is typically the first step in small RNA-seq analysis.
    """)

    # Check cutadapt installation
    cutadapt_ok, version = check_cutadapt()

    if cutadapt_ok:
        st.success(f"✅ Cutadapt {version} installed")
    else:
        st.error("""
        ❌ Cutadapt not found. Please install it:
        ```bash
        pip install cutadapt
        # or
        conda install -c bioconda cutadapt
        ```
        """)
        return

    tab1, tab2, tab3 = st.tabs([
        "📁 Input Files",
        "⚙️ Settings",
        "📊 Results"
    ])

    with tab1:
        render_trimming_input()

    with tab2:
        render_trimming_settings()

    with tab3:
        render_trimming_results()


def render_trimming_input():
    """Render input file section"""
    st.subheader("Input FASTQ Files")

    # Check if files are available from Project module
    project_files = st.session_state.get('project_fastq_files', [])
    project_names = st.session_state.get('project_fastq_names', [])
    has_project_files = project_files and len(project_files) > 0

    if has_project_files:
        st.success(f"✅ **{len(project_files)} files available from Project**")

        # Let user choose source
        input_source = st.radio(
            "Input Source",
            options=["Use files from Project", "Upload new files"],
            index=0,
            horizontal=True
        )

        if input_source == "Use files from Project":
            # Show project files
            file_info = []
            for path, name in zip(project_files, project_names):
                if isinstance(path, Path) and path.exists():
                    size_mb = path.stat().st_size / (1024 * 1024)
                    file_info.append({
                        'Filename': name,
                        'Size (MB)': f"{size_mb:.2f}",
                        'Status': '✅ Ready'
                    })

            if file_info:
                st.dataframe(pd.DataFrame(file_info), width="stretch")

            # Store project files as input
            st.session_state.trim_input_files = project_files
            st.session_state.trim_use_project_files = True
            return
        else:
            st.info("Upload new FASTQ files below (these will be used instead of Project files)")

    else:
        st.info("📁 Upload your FASTQ files here, or load them first in the **Project** module.")

    # File uploader for direct upload
    uploaded_files = st.file_uploader(
        "Upload FASTQ files",
        type=None,  # Accept all types - filter manually for .fastq.gz support
        accept_multiple_files=True,
        help="Upload raw FASTQ files (.fastq, .fq, .fastq.gz, .fq.gz)"
    )

    # Filter to valid FASTQ files
    if uploaded_files:
        valid_extensions = ('.fastq', '.fq', '.fastq.gz', '.fq.gz')
        valid_files = []
        for f in uploaded_files:
            fname = f.name.lower()
            if fname.endswith(valid_extensions):
                valid_files.append(f)
            else:
                st.warning(f"⚠️ Skipping '{f.name}' - not a valid FASTQ file")
        uploaded_files = valid_files if valid_files else None

    if uploaded_files:
        st.success(f"✅ Uploaded {len(uploaded_files)} files")

        # Show file list
        file_info = []
        for f in uploaded_files:
            size_mb = f.size / (1024 * 1024)
            file_info.append({
                'Filename': f.name,
                'Size (MB)': f"{size_mb:.2f}"
            })

        st.dataframe(pd.DataFrame(file_info), width="stretch")

        # Store in session state
        st.session_state.trim_input_files = uploaded_files
        st.session_state.trim_use_project_files = False


def render_trimming_settings():
    """Render trimming settings section"""
    st.subheader("Trimming Settings")

    # Adapter selection - organized by platform
    st.markdown("### Adapter Sequences")

    # Platform selection first
    platform = st.selectbox(
        "Sequencing Platform",
        options=list(ADAPTER_PLATFORMS.keys()),
        index=0,  # Default to Illumina
        help="Select your sequencing platform to see available adapter presets"
    )

    # Get adapters for selected platform
    platform_adapters = ADAPTER_PLATFORMS[platform]

    # Adapter preset selection
    preset = st.selectbox(
        "Adapter Preset",
        options=platform_adapters,
        help="Select your library prep kit or choose Custom"
    )

    preset_info = ADAPTER_PRESETS[preset]

    # Show platform and description
    col_info1, col_info2 = st.columns([1, 3])
    with col_info1:
        st.markdown(f"**Platform:** {preset_info.get('platform', 'Unknown')}")
    with col_info2:
        st.caption(preset_info['description'])

    # Adapter sequence inputs
    col1, col2 = st.columns(2)

    with col1:
        # 3' Adapter
        default_3 = preset_info.get('3_adapter', '') or ''
        if preset == "Custom":
            adapter_3 = st.text_input(
                "3' Adapter Sequence",
                value="",
                help="Enter your 3' adapter sequence"
            )
        else:
            adapter_3 = st.text_input(
                "3' Adapter Sequence",
                value=default_3,
                help="3' adapter sequence (can be modified)"
            )

    with col2:
        # 5' Adapter - check if preset has one
        preset_5_adapter = preset_info.get('5_adapter')
        has_preset_5 = preset_5_adapter is not None and preset_5_adapter != ""

        if has_preset_5:
            st.info(f"ℹ️ This kit includes a 5' adapter")

        use_5_adapter = st.checkbox(
            "Use 5' adapter",
            value=has_preset_5,
            help="Enable if your library prep uses a 5' adapter"
        )

        if use_5_adapter:
            default_5 = preset_5_adapter if has_preset_5 else ""
            adapter_5 = st.text_input(
                "5' Adapter Sequence",
                value=default_5,
                help="Enter 5' adapter sequence"
            )
        else:
            adapter_5 = None

    st.session_state.trim_adapter_3 = adapter_3
    st.session_state.trim_adapter_5 = adapter_5

    # Platform-specific notes
    if platform == "Nanopore":
        st.info("💡 **Nanopore Note:** Direct RNA sequencing often doesn't require adapter trimming. "
                "For cDNA-PCR libraries, use the cDNA-PCR preset.")
    elif platform == "PacBio":
        st.info("💡 **PacBio Note:** HiFi/CCS reads are typically pre-processed. "
                "Use this for raw subreads if needed.")
    elif platform == "Ion Torrent":
        st.info("💡 **Ion Torrent Note:** The small RNA kit adapter differs from standard Ion adapters.")
    elif platform == "BGI/MGI":
        st.info("💡 **BGI/MGI Note:** DNBSEQ and MGISeq platforms use different adapters than Illumina.")

    # Length filtering
    st.markdown("### Length Filtering")

    col1, col2 = st.columns(2)

    with col1:
        min_length = st.slider(
            "Minimum Length",
            10, 50, 18,
            help="Reads shorter than this will be discarded (18 nt typical for miRNA)"
        )

    with col2:
        max_length = st.slider(
            "Maximum Length",
            25, 200, 35,
            help="Reads longer than this will be discarded (35 nt for miRNA, increase for piRNA/snoRNA)"
        )

    st.session_state.trim_min_length = min_length
    st.session_state.trim_max_length = max_length

    # Quality settings
    st.markdown("### Quality Settings")

    col1, col2 = st.columns(2)

    with col1:
        quality_cutoff = st.slider(
            "Quality Cutoff (Phred)",
            0, 30, 20,
            help="Trim low-quality bases from 3' end"
        )

    with col2:
        error_rate = st.slider(
            "Adapter Error Rate",
            0.0, 0.3, 0.1, 0.01,
            help="Maximum allowed error rate for adapter matching"
        )

    st.session_state.trim_quality = quality_cutoff
    st.session_state.trim_error_rate = error_rate

    # Advanced options
    with st.expander("Advanced Options"):
        col1, col2 = st.columns(2)

        with col1:
            trim_n = st.checkbox("Trim N bases", value=True)
            discard_untrimmed = st.checkbox(
                "Discard reads without adapter",
                value=False,
                help="Only keep reads where adapter was found"
            )

        with col2:
            cores = st.slider("CPU Cores", 1, 16, 4)

    st.session_state.trim_n = trim_n
    st.session_state.trim_discard = discard_untrimmed
    st.session_state.trim_cores = cores

    # Run trimming
    st.divider()

    input_files = st.session_state.get('trim_input_files')

    if not input_files:
        st.warning("Please upload FASTQ files in the Input tab first.")
        return

    if not adapter_3:
        st.warning("Please specify a 3' adapter sequence.")
        return

    if st.button("✂️ Run Trimming", type="primary", width="stretch"):
        # Create temp directory for output
        temp_dir = Path(tempfile.mkdtemp())
        output_dir = temp_dir / "output"
        output_dir.mkdir()

        # Handle input files - could be Path objects (from Project) or UploadedFile objects
        input_paths = []
        use_project = st.session_state.get('trim_use_project_files', False)

        if use_project:
            # Files are already Path objects from Project module
            for file_path in input_files:
                if isinstance(file_path, Path) and file_path.exists():
                    input_paths.append(file_path)
        else:
            # Files are UploadedFile objects - save to temp directory
            input_dir = temp_dir / "input"
            input_dir.mkdir()

            for uf in input_files:
                file_path = input_dir / uf.name
                with open(file_path, 'wb') as f:
                    f.write(uf.getbuffer())
                input_paths.append(file_path)

        if not input_paths:
            st.error("No valid input files found")
            return

        # Progress tracking
        progress_bar = st.progress(0)
        status_text = st.empty()

        def update_progress(progress, message):
            progress_bar.progress(progress)
            status_text.text(message)

        with st.spinner("Running Cutadapt..."):
            output_files, stats_df = batch_trim_files(
                input_files=input_paths,
                output_dir=output_dir,
                adapter_3=adapter_3,
                adapter_5=adapter_5,
                min_length=min_length,
                max_length=max_length,
                quality_cutoff=quality_cutoff,
                cores=cores,
                progress_callback=update_progress
            )

        progress_bar.progress(1.0)
        status_text.text("Trimming complete!")

        # Store results
        st.session_state.trim_output_files = output_files
        st.session_state.trim_stats = stats_df
        st.session_state.trim_output_dir = output_dir

        if len(output_files) > 0:
            st.success(f"Successfully trimmed {len(output_files)} files!")
            st.info("Go to the Results tab to view statistics and download files.")
        else:
            st.error("Trimming failed. Check the error messages above.")


def render_trimming_results():
    """Render trimming results section"""
    st.subheader("Trimming Results")

    stats_df = st.session_state.get('trim_stats')
    output_files = st.session_state.get('trim_output_files')

    if stats_df is None:
        st.info("No trimming results yet. Run trimming first.")
        return

    # Summary metrics
    col1, col2, col3, col4 = st.columns(4)

    successful = (stats_df['status'] == 'success').sum()
    total_reads = stats_df['total_reads'].sum()
    reads_written = stats_df['reads_written'].sum()
    avg_adapter_rate = stats_df['adapter_rate'].mean()

    with col1:
        st.metric("Files Processed", f"{successful}/{len(stats_df)}")

    with col2:
        st.metric("Total Reads", f"{total_reads:,}")

    with col3:
        st.metric("Reads Passed", f"{reads_written:,}")

    with col4:
        st.metric("Avg Adapter Rate", f"{avg_adapter_rate:.1f}%")

    # Check for issues (0 reads processed)
    if total_reads == 0:
        st.error("⚠️ **No reads were processed!** This usually means:")
        st.markdown("""
        1. **Files are empty or corrupted** - Check your input files
        2. **File format issue** - Ensure files are valid FASTQ format
        3. **Cutadapt can't read the files** - Check debug info below
        """)

        # Show debug info if available
        with st.expander("🔍 Debug Information"):
            for _, row in stats_df.iterrows():
                st.markdown(f"**{row['sample']}:**")
                if 'input_file_size' in row and pd.notna(row.get('input_file_size')):
                    st.write(f"- Input file size: {row['input_file_size']:,} bytes")
                if 'output_size' in row and pd.notna(row.get('output_size')):
                    st.write(f"- Output file size: {row['output_size']:,} bytes")
                if 'debug_cmd' in row and row.get('debug_cmd'):
                    st.write("**Command:**")
                    st.code(row['debug_cmd'], language='bash')
                if 'debug_stderr' in row and row.get('debug_stderr'):
                    st.write("**Cutadapt output (stderr):**")
                    st.code(row['debug_stderr'], language='text')
                if 'debug_stdout' in row and row.get('debug_stdout'):
                    st.write("**Cutadapt output (stdout):**")
                    st.code(row['debug_stdout'], language='text')
                st.divider()

    # Detailed statistics table
    st.subheader("Per-Sample Statistics")

    display_cols = ['sample', 'total_reads', 'reads_with_adapter', 'reads_written',
                    'reads_too_short', 'adapter_rate', 'pass_rate', 'status']
    display_cols = [c for c in display_cols if c in stats_df.columns]

    st.dataframe(
        stats_df[display_cols].style.format({
            'adapter_rate': '{:.1f}%',
            'pass_rate': '{:.1f}%',
            'total_reads': '{:,}',
            'reads_with_adapter': '{:,}',
            'reads_written': '{:,}',
            'reads_too_short': '{:,}'
        }),
        width="stretch"
    )

    # Visualization
    st.subheader("Visualization")

    import plotly.express as px
    import plotly.graph_objects as go

    col1, col2 = st.columns(2)

    with col1:
        # Read retention plot
        fig = go.Figure()

        fig.add_trace(go.Bar(
            name='Passed',
            x=stats_df['sample'],
            y=stats_df['reads_written'],
            marker_color='#00A087'
        ))

        fig.add_trace(go.Bar(
            name='Too Short',
            x=stats_df['sample'],
            y=stats_df['reads_too_short'],
            marker_color='#E64B35'
        ))

        fig.update_layout(
            title="Read Retention",
            barmode='stack',
            xaxis_title="Sample",
            yaxis_title="Reads"
        )

        st.plotly_chart(fig, width="stretch")

    with col2:
        # Adapter detection rate
        fig = px.bar(
            stats_df,
            x='sample',
            y='adapter_rate',
            title="Adapter Detection Rate",
            color='adapter_rate',
            color_continuous_scale='RdYlGn'
        )
        fig.update_layout(
            xaxis_title="Sample",
            yaxis_title="Adapter Rate (%)"
        )
        st.plotly_chart(fig, width="stretch")

    # Download section
    st.divider()
    st.subheader("Download Results")

    col1, col2 = st.columns(2)

    with col1:
        # Download statistics
        csv = stats_df.to_csv(index=False)
        st.download_button(
            "📥 Download Statistics (CSV)",
            csv,
            "trimming_statistics.csv",
            "text/csv"
        )

    with col2:
        # Download trimmed files
        if output_files:
            import io
            import zipfile

            buffer = io.BytesIO()
            with zipfile.ZipFile(buffer, 'w', zipfile.ZIP_DEFLATED) as zf:
                for of in output_files:
                    if of.exists():
                        zf.write(of, of.name)

            st.download_button(
                "📥 Download Trimmed Files (ZIP)",
                buffer.getvalue(),
                "trimmed_fastq.zip",
                "application/zip"
            )

    # Option to use trimmed files in pipeline
    st.divider()
    if st.button("➡️ Use Trimmed Files for QC/Alignment"):
        st.session_state.fastq_files = output_files
        st.success("Trimmed files are now ready for QC and alignment!")
