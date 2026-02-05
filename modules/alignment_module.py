"""
Alignment Module for sRNAtlas
Handles Bowtie alignment optimized for small RNA sequences
"""
import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
import subprocess
import tempfile
import os
from typing import List, Dict, Optional
from datetime import datetime

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from config import config


def check_bowtie_installed() -> bool:
    """Check if bowtie is installed and available"""
    try:
        result = subprocess.run(['bowtie', '--version'], capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False


def check_samtools_installed() -> bool:
    """Check if samtools is installed and available"""
    try:
        result = subprocess.run(['samtools', '--version'], capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False


def validate_bam_file(bam_file: Path) -> Dict:
    """
    Validate a BAM file is properly formatted and not corrupted

    Args:
        bam_file: Path to BAM file

    Returns:
        Dictionary with validation status
    """
    result = {
        'valid': False,
        'error': None,
        'read_count': 0
    }

    try:
        # Check file exists and has content
        if not bam_file.exists():
            result['error'] = "BAM file does not exist"
            return result

        if bam_file.stat().st_size == 0:
            result['error'] = "BAM file is empty"
            return result

        # Check BAM magic bytes (should start with gzip header for BGZF)
        with open(bam_file, 'rb') as f:
            header = f.read(4)
            if header[:2] != b'\x1f\x8b':
                result['error'] = "Invalid BAM file (not BGZF compressed)"
                return result

        # Use samtools quickcheck if available
        quickcheck = subprocess.run(
            ['samtools', 'quickcheck', str(bam_file)],
            capture_output=True,
            text=True
        )

        if quickcheck.returncode != 0:
            result['error'] = f"BAM file failed quickcheck: {quickcheck.stderr}"
            return result

        # Count reads to verify file is readable
        count_result = subprocess.run(
            ['samtools', 'view', '-c', str(bam_file)],
            capture_output=True,
            text=True
        )

        if count_result.returncode == 0:
            result['read_count'] = int(count_result.stdout.strip())
            result['valid'] = True
        else:
            result['error'] = f"Failed to read BAM: {count_result.stderr}"

    except Exception as e:
        result['error'] = str(e)

    return result


def build_bowtie_index(fasta_file: Path, output_prefix: Path) -> Dict:
    """
    Build Bowtie index from FASTA file

    Args:
        fasta_file: Path to reference FASTA
        output_prefix: Prefix for index files

    Returns:
        Dictionary with build status and info
    """
    cmd = [
        'bowtie-build',
        '--threads', str(config.alignment.threads),
        str(fasta_file),
        str(output_prefix)
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

        if result.returncode == 0:
            return {
                'status': 'success',
                'index_prefix': str(output_prefix),
                'message': 'Index built successfully'
            }
        else:
            return {
                'status': 'error',
                'error': result.stderr
            }

    except subprocess.TimeoutExpired:
        return {
            'status': 'error',
            'error': 'Index building timed out after 1 hour'
        }
    except Exception as e:
        return {
            'status': 'error',
            'error': str(e)
        }


def run_alignment(
    fastq_file: Path,
    index_prefix: Path,
    output_dir: Path,
    settings: Optional[Dict] = None
) -> Dict:
    """
    Run Bowtie alignment for a single sample

    Args:
        fastq_file: Path to FASTQ file
        index_prefix: Bowtie index prefix
        output_dir: Output directory for BAM files
        settings: Alignment settings (uses config if None)

    Returns:
        Dictionary with alignment results
    """
    if settings is None:
        settings = config.alignment

    # Extract sample name
    sample_name = fastq_file.stem
    for ext in ['.fastq', '.fq', '.gz']:
        sample_name = sample_name.replace(ext, '')

    output_dir.mkdir(parents=True, exist_ok=True)

    sam_file = output_dir / f"{sample_name}.sam"
    bam_file = output_dir / f"{sample_name}_sorted.bam"
    log_file = output_dir / f"{sample_name}_alignment.log"
    unaligned_file = output_dir / f"{sample_name}_unaligned.fq"

    # Build bowtie command
    # Bowtie1 syntax: bowtie [options] <index_prefix> <reads> [output]
    # All options MUST come before positional arguments

    # Start with command and options
    cmd = ['bowtie']

    # Threading
    cmd.extend(['-p', str(settings.threads)])

    # Mismatch settings
    cmd.extend(['-v', str(settings.mismatches)])  # Allow n mismatches

    # Multi-mapping settings (must be before positional args)
    if hasattr(settings, 'report_all') and settings.report_all:
        cmd.append('-a')
    else:
        cmd.extend(['-k', str(settings.max_alignments)])

    if hasattr(settings, 'suppress_multi') and settings.suppress_multi > 0:
        cmd.extend(['-m', str(settings.suppress_multi)])

    if hasattr(settings, 'best_mode') and settings.best_mode:
        cmd.append('--best')

    if hasattr(settings, 'strata') and settings.strata:
        cmd.append('--strata')

    # Output unaligned reads
    cmd.extend(['--un', str(unaligned_file)])

    # SAM output flag - with -S, bowtie writes SAM to stdout
    cmd.append('-S')

    # Positional arguments: <index_prefix> <reads>
    # Note: With -S flag, SAM output goes to stdout, NOT to a file argument
    cmd.append(str(index_prefix))  # Index prefix
    cmd.append(str(fastq_file))    # Input FASTQ
    # Do NOT add sam_file here - Bowtie1 with -S writes SAM to stdout

    result = {
        'sample': sample_name,
        'status': 'pending',
        'bam_file': None,
        'stats': None,
        'error': None
    }

    try:
        # Run alignment
        # With -S flag, bowtie writes SAM to stdout and stats to stderr
        # Redirect stdout to SAM file, capture stderr for stats
        with open(sam_file, 'w') as sam_out:
            align_result = subprocess.run(
                cmd,
                stdout=sam_out,  # SAM output goes here
                stderr=subprocess.PIPE,
                text=True,
                timeout=7200  # 2 hour timeout
            )

        # Write alignment stats (stderr) to log file
        with open(log_file, 'w') as log:
            log.write(f"Command: {' '.join(cmd)}\n\n")
            log.write("Bowtie Alignment Statistics:\n")
            log.write(align_result.stderr)

        if align_result.returncode != 0:
            result['status'] = 'error'
            result['error'] = align_result.stderr
            return result

        # Parse alignment statistics from stderr
        stats = parse_bowtie_log(align_result.stderr)

        # Convert SAM to BAM (two-step process to avoid truncation issues)
        unsorted_bam = output_dir / f"{sample_name}_unsorted.bam"

        # Step 1: Convert SAM to BAM
        view_result = subprocess.run(
            ['samtools', 'view', '-bS', '-o', str(unsorted_bam), str(sam_file)],
            capture_output=True,
            text=True
        )
        if view_result.returncode != 0:
            result['status'] = 'error'
            result['error'] = f"SAM to BAM conversion failed: {view_result.stderr}"
            return result

        # Step 2: Sort BAM
        sort_result = subprocess.run(
            ['samtools', 'sort', '-@', str(settings.threads), '-o', str(bam_file), str(unsorted_bam)],
            capture_output=True,
            text=True
        )
        if sort_result.returncode != 0:
            result['status'] = 'error'
            result['error'] = f"BAM sorting failed: {sort_result.stderr}"
            return result

        # Remove unsorted BAM
        if unsorted_bam.exists():
            unsorted_bam.unlink()

        # Index BAM
        index_result = subprocess.run(
            ['samtools', 'index', str(bam_file)],
            capture_output=True,
            text=True
        )
        if index_result.returncode != 0:
            result['status'] = 'error'
            result['error'] = f"BAM indexing failed: {index_result.stderr}"
            return result

        # Validate BAM file before proceeding
        validation = validate_bam_file(bam_file)
        if not validation['valid']:
            result['status'] = 'error'
            result['error'] = f"BAM validation failed: {validation['error']}"
            return result

        # Remove SAM file to save space
        if sam_file.exists():
            sam_file.unlink()

        # Get BAM statistics
        flagstat_result = subprocess.run(
            ['samtools', 'flagstat', str(bam_file)],
            capture_output=True,
            text=True
        )

        result['status'] = 'success'
        result['bam_file'] = str(bam_file)
        result['stats'] = stats
        result['flagstat'] = flagstat_result.stdout
        result['bam_read_count'] = validation['read_count']

    except subprocess.TimeoutExpired:
        result['status'] = 'error'
        result['error'] = 'Alignment timed out'
    except Exception as e:
        result['status'] = 'error'
        result['error'] = str(e)

    return result


def parse_bowtie_log(log_text: str) -> Dict:
    """
    Parse Bowtie stderr output for statistics

    Bowtie output format (varies by version):
    # reads processed: 1000000
    # reads with at least one alignment: 850000 (85.00%)
    # reads that failed to align: 100000 (10.00%)
    # reads with alignments suppressed due to -m: 50000 (5.00%)
    Reported 1234567 alignments

    Args:
        log_text: Bowtie stderr output

    Returns:
        Dictionary with alignment statistics
    """
    stats = {
        'total_reads': 0,
        'aligned_reads': 0,
        'not_aligned': 0,
        'suppressed': 0,
        'alignment_rate': 0.0,
        'total_alignments': 0
    }

    for line in log_text.split('\n'):
        line = line.strip()

        if 'reads processed:' in line:
            try:
                stats['total_reads'] = int(line.split(':')[1].strip())
            except (ValueError, IndexError):
                pass

        # Match both "at least one alignment" and "at least one reported alignment"
        elif 'reads with at least one' in line and 'alignment' in line:
            try:
                num_str = line.split(':')[1].split('(')[0].strip()
                stats['aligned_reads'] = int(num_str)
            except (ValueError, IndexError):
                pass

        elif 'reads that failed to align:' in line:
            try:
                num_str = line.split(':')[1].split('(')[0].strip()
                stats['not_aligned'] = int(num_str)
            except (ValueError, IndexError):
                pass

        elif 'alignments suppressed due to -m:' in line:
            try:
                num_str = line.split(':')[1].split('(')[0].strip()
                stats['suppressed'] = int(num_str)
            except (ValueError, IndexError):
                pass

        elif line.startswith('Reported') and 'alignments' in line:
            try:
                # "Reported 4151903 alignments"
                stats['total_alignments'] = int(line.split()[1])
            except (ValueError, IndexError):
                pass

    # Calculate alignment rate
    if stats['total_reads'] > 0:
        stats['alignment_rate'] = (stats['aligned_reads'] / stats['total_reads']) * 100

    return stats


def run_batch_alignment(
    fastq_files: List[Path],
    index_prefix: Path,
    output_dir: Path,
    progress_callback=None
) -> Dict:
    """
    Run alignment for multiple samples

    Args:
        fastq_files: List of FASTQ files
        index_prefix: Bowtie index prefix
        output_dir: Output directory
        progress_callback: Optional callback for progress updates

    Returns:
        Dictionary with all alignment results
    """
    results = {}
    total = len(fastq_files)

    for i, fastq in enumerate(fastq_files):
        if progress_callback:
            progress_callback(i / total, f"Aligning {fastq.name}...")

        result = run_alignment(fastq, index_prefix, output_dir)
        results[result['sample']] = result

    return results


def render_alignment_page():
    """Render the alignment page"""
    st.header("🔗 Alignment")

    st.info("""
    **Bowtie Aligner** - Optimized for small RNA sequences (18-35 nt)

    Bowtie is preferred over Bowtie2 for small RNA because:
    - Designed for short reads (<50 bp)
    - Exact mismatch control with `-v` mode
    - No gap alignment (appropriate for miRNA/siRNA)
    - Faster for very short sequences
    """)

    # Check dependencies
    col1, col2 = st.columns(2)

    with col1:
        bowtie_ok = check_bowtie_installed()
        if bowtie_ok:
            st.success("✅ Bowtie installed")
        else:
            st.error("❌ Bowtie not found. Please install it first.")

    with col2:
        samtools_ok = check_samtools_installed()
        if samtools_ok:
            st.success("✅ Samtools installed")
        else:
            st.error("❌ Samtools not found. Please install it first.")

    if not (bowtie_ok and samtools_ok):
        st.warning("Please install the required tools to use the alignment module.")
        st.code("""
# Install with conda
conda install -c bioconda bowtie samtools

# Or with apt (Ubuntu/Debian)
sudo apt-get install bowtie samtools
        """)
        return

    tab1, tab2, tab3 = st.tabs(["📁 Setup", "▶️ Run Alignment", "📊 Results"])

    with tab1:
        render_alignment_setup()

    with tab2:
        render_run_alignment()

    with tab3:
        render_alignment_results()


def render_alignment_setup():
    """Render alignment setup section"""
    st.subheader("Reference Setup")

    # Refresh button to rescan databases
    col_refresh, col_spacer = st.columns([1, 4])
    with col_refresh:
        if st.button("🔄 Refresh", help="Rescan reference directories"):
            from utils.reference_manager import get_reference_manager
            ref_manager = get_reference_manager()
            ref_manager.refresh()
            st.rerun()

    ref_option = st.radio(
        "Reference Source",
        options=["Available References", "Use existing index", "Build from FASTA", "Go to Database Downloads"]
    )

    if ref_option == "Go to Database Downloads":
        st.info("👉 Go to **🗄️ Databases** in the sidebar to download miRBase or RNAcentral references.")
        return

    if ref_option == "Available References":
        # Get reference manager
        from utils.reference_manager import get_reference_manager
        ref_manager = get_reference_manager()
        ref_manager._auto_detect_references()  # Rescan for new files
        references = ref_manager.list_references()

        if not references:
            st.warning("""
            No reference databases found.

            **To add references:**
            1. Go to **🗄️ Databases** in the sidebar
            2. Download from miRBase or RNAcentral
            3. Return here and click 'Refresh'
            """)
            return

        # Show available references with index status
        st.markdown("**Available Reference Databases:**")

        # Create a nice display
        ref_data = []
        for ref in references:
            ref_path = ref_manager.get_reference_path(ref['id'])
            has_index = False
            if ref_path:
                index_prefix = ref_path.with_suffix('')
                has_index = (index_prefix.parent / f"{index_prefix.name}.1.ebwt").exists()

            ref_data.append({
                'id': ref['id'],
                'name': ref.get('name', ref['id']),
                'species': ref.get('species', 'Unknown'),
                'source': ref.get('source', 'Custom'),
                'sequences': ref.get('sequences', 0),
                'has_index': has_index,
                'path': str(ref_path) if ref_path else None
            })

        # Create selectbox with descriptive options
        ref_options = {
            r['id']: f"{'✅' if r['has_index'] else '⚠️'} {r['name']} ({r['species']}, {r['sequences']:,} seqs)"
            for r in ref_data
        }

        selected_id = st.selectbox(
            "Select Reference",
            options=list(ref_options.keys()),
            format_func=lambda x: ref_options[x],
            key="alignment_ref_selector",
            help="✅ = Index ready, ⚠️ = Index needs to be built"
        )

        if selected_id:
            selected_ref = next((r for r in ref_data if r['id'] == selected_id), None)

            if selected_ref:
                # Show details
                with st.expander("📋 Reference Details", expanded=False):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.markdown(f"**Name:** {selected_ref['name']}")
                        st.markdown(f"**Species:** {selected_ref['species']}")
                        st.markdown(f"**Source:** {selected_ref['source']}")
                    with col2:
                        st.markdown(f"**Sequences:** {selected_ref['sequences']:,}")
                        st.markdown(f"**Index Ready:** {'Yes ✅' if selected_ref['has_index'] else 'No ⚠️'}")
                    if selected_ref['path']:
                        st.caption(f"Path: {selected_ref['path']}")

                # Store in session state
                st.session_state.selected_reference = selected_ref
                st.session_state.reference_fasta = selected_ref['path']

                if selected_ref['has_index']:
                    index_prefix = Path(selected_ref['path']).with_suffix('')
                    st.success(f"✅ Ready for alignment! Index exists for {selected_ref['name']}")
                    st.session_state.reference_index = str(index_prefix)
                else:
                    st.warning(f"⚠️ Bowtie index not found. Please build the index first.")

                    if st.button("🔨 Build Bowtie Index", type="primary", key="build_idx_btn"):
                        with st.spinner(f"Building index for {selected_ref['name']}... This may take a while."):
                            ref_path = Path(selected_ref['path'])
                            index_prefix = ref_path.with_suffix('')
                            result = build_bowtie_index(ref_path, index_prefix)

                            if result['status'] == 'success':
                                st.success("✅ Bowtie index built successfully!")
                                st.session_state.reference_index = str(index_prefix)
                                st.rerun()
                            else:
                                st.error(f"Error building index: {result['error']}")

    elif ref_option == "Use existing index":
        index_path = st.text_input(
            "Bowtie Index Prefix",
            placeholder="/path/to/index/reference",
            help="Path to existing Bowtie index (without .ebwt extension)"
        )

        if index_path:
            # Check if index exists (Bowtie uses .ebwt extension)
            if Path(f"{index_path}.1.ebwt").exists():
                st.success("✅ Bowtie index found")
                st.session_state.reference_index = index_path
            elif Path(f"{index_path}.1.ebwtl").exists():
                st.success("✅ Bowtie large index found")
                st.session_state.reference_index = index_path
            else:
                st.error("Index files not found at this location")

    elif ref_option == "Build from FASTA":
        fasta_file = st.file_uploader(
            "Upload Reference FASTA",
            type=['fasta', 'fa', 'fna', 'gz']
        )

        if fasta_file:
            output_name = st.text_input(
                "Index Name",
                value=fasta_file.name.replace('.fasta', '').replace('.fa', '').replace('.gz', '')
            )

            if st.button("🔨 Build Index", type="primary"):
                with st.spinner("Building Bowtie index... This may take a while."):
                    # Save FASTA file
                    temp_dir = Path(tempfile.mkdtemp())
                    fasta_path = temp_dir / fasta_file.name
                    with open(fasta_path, 'wb') as f:
                        f.write(fasta_file.getbuffer())

                    # Build index
                    index_prefix = temp_dir / output_name
                    result = build_bowtie_index(fasta_path, index_prefix)

                    if result['status'] == 'success':
                        st.success("Index built successfully!")
                        st.session_state.reference_index = str(index_prefix)
                    else:
                        st.error(f"Error building index: {result['error']}")

    # Alignment Settings
    st.divider()
    st.subheader("⚙️ Alignment Settings (Bowtie)")

    col1, col2 = st.columns(2)

    with col1:
        threads = st.slider("Threads", 1, 32, config.alignment.threads)

        mismatches = st.selectbox(
            "Mismatches Allowed (-v)",
            options=[0, 1, 2, 3],
            index=1,
            help="Number of mismatches allowed in the entire read. For small RNA, 0-1 is recommended."
        )

        max_alignments = st.slider(
            "Max Alignments (-k)",
            1, 50, config.alignment.max_alignments,
            help="Report up to k alignments per read"
        )

    with col2:
        best_mode = st.checkbox(
            "Best mode (--best)",
            value=config.alignment.best_mode,
            help="Report alignments in best-to-worst order"
        )

        strata = st.checkbox(
            "Strata mode (--strata)",
            value=config.alignment.strata,
            help="Only report alignments in the best stratum (requires --best)"
        )

        suppress_multi = st.number_input(
            "Suppress if >m alignments (-m)",
            min_value=0, max_value=100, value=0,
            help="Suppress reads with more than m alignments (0=disabled)"
        )

    # Save settings
    if st.button("💾 Save Settings"):
        config.alignment.threads = threads
        config.alignment.mismatches = mismatches
        config.alignment.max_alignments = max_alignments
        config.alignment.best_mode = best_mode
        config.alignment.strata = strata
        config.alignment.suppress_multi = suppress_multi

        st.session_state.alignment_settings = config.alignment
        st.success("Settings saved!")


def render_run_alignment():
    """Render alignment execution section"""
    st.subheader("Run Alignment")

    # Check prerequisites
    if not st.session_state.get('reference_index'):
        st.warning("Please set up a reference index in the Setup tab first.")
        return

    st.success(f"Reference: {st.session_state.reference_index}")

    # Check if trimmed files are available from previous step
    trim_output_files = st.session_state.get('trim_output_files')
    trim_output_dir = st.session_state.get('trim_output_dir')
    trimmed_files_available = trim_output_files and len(trim_output_files) > 0

    # Check if project files (untrimmed) are available
    project_files = st.session_state.get('project_fastq_files', [])
    project_files_available = project_files and len(project_files) > 0

    # Build input source options
    input_options = []
    if trimmed_files_available:
        input_options.append("Use trimmed files")
    if project_files_available:
        input_options.append("Use untrimmed files (from Project)")
    input_options.append("Upload new files")

    if trimmed_files_available:
        st.success(f"✅ {len(trim_output_files)} trimmed files available from Trimming step")
    if project_files_available:
        st.info(f"📁 {len(project_files)} untrimmed files available from Project")

    input_source = st.radio(
        "Input Source",
        options=input_options,
        index=0,
        help="Use trimmed files, untrimmed project files, or upload new FASTQ files"
    )

    if input_source == "Use trimmed files":
        st.info(f"Using trimmed files from: {trim_output_dir or 'trimming output'}")

        # Show file list
        file_info = [{'Filename': Path(f).name, 'Source': 'Trimming Step'} for f in trim_output_files]
        st.dataframe(pd.DataFrame(file_info), width="stretch")

        st.session_state.alignment_input_files = [Path(f) for f in trim_output_files]
        st.session_state.alignment_input_source = 'trimming'

    elif input_source == "Use untrimmed files (from Project)":
        st.warning("⚠️ Using untrimmed files - adapters may affect alignment quality")

        # Show file list
        project_names = st.session_state.get('project_fastq_names', [])
        file_info = [{'Filename': Path(f).name, 'Original': project_names[i] if i < len(project_names) else 'N/A', 'Source': 'Project'}
                     for i, f in enumerate(project_files)]
        st.dataframe(pd.DataFrame(file_info), width="stretch")

        st.session_state.alignment_input_files = [Path(f) for f in project_files]
        st.session_state.alignment_input_source = 'project'

    else:
        st.session_state.alignment_input_source = 'upload'

    # Upload FASTQ files (if not using trimmed)
    uploaded_files = None
    if st.session_state.get('alignment_input_source') == 'upload':
        uploaded_files = st.file_uploader(
            "Upload FASTQ files",
            type=['fastq', 'fq', 'gz'],
            accept_multiple_files=True
        )

    if uploaded_files:
        st.info(f"Selected {len(uploaded_files)} files for alignment")

        # Show file list
        with st.expander("View files"):
            for f in uploaded_files:
                st.text(f.name)

    # Determine if we have files to align
    input_source = st.session_state.get('alignment_input_source')
    has_files = uploaded_files or input_source in ['trimming', 'project']

    # Output directory (always show when we have files)
    if has_files:
        output_dir = st.text_input(
            "Output Directory",
            value="./output/alignment",
            help="Directory to save BAM files",
            key="alignment_output_dir_input"
        )

    if has_files:
        if st.button("🚀 Start Alignment", type="primary"):
            # Create progress tracking
            progress_bar = st.progress(0)
            status_text = st.empty()

            # Get input files based on source
            fastq_files = []

            if input_source == 'trimming':
                # Use files from trimming step
                fastq_files = st.session_state.get('alignment_input_files', [])
                status_text.text("Using trimmed files from previous step...")
            elif input_source == 'project':
                # Use untrimmed files from project
                fastq_files = st.session_state.get('alignment_input_files', [])
                status_text.text("Using untrimmed files from project...")
            else:
                # Save uploaded files temporarily
                temp_dir = Path(tempfile.mkdtemp())
                for uf in uploaded_files:
                    file_path = temp_dir / uf.name
                    with open(file_path, 'wb') as f:
                        f.write(uf.getbuffer())
                    fastq_files.append(file_path)

            if not fastq_files:
                st.error("No input files found!")
                return

            # Run alignment
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

            results = {}
            bam_files = []
            for i, fastq in enumerate(fastq_files):
                progress = (i + 1) / len(fastq_files)
                progress_bar.progress(progress)
                status_text.text(f"Aligning {fastq.name}... ({i+1}/{len(fastq_files)})")

                result = run_alignment(
                    fastq,
                    Path(st.session_state.reference_index),
                    output_path
                )
                results[result['sample']] = result
                if result.get('bam_file'):
                    bam_files.append(result['bam_file'])

            # Save results for next step (Counting)
            st.session_state.alignment_results = results
            st.session_state.alignment_bam_files = bam_files
            st.session_state.alignment_output_dir = str(output_path)

            # Show summary
            n_success = sum(1 for r in results.values() if r['status'] == 'success')
            n_failed = len(results) - n_success

            if n_failed == 0:
                st.success(f"✅ All {n_success} samples aligned successfully!")
            else:
                st.warning(f"Completed: {n_success} success, {n_failed} failed")

            st.balloons()


def render_alignment_results():
    """Render alignment results section"""
    st.subheader("Alignment Results")

    if not st.session_state.get('alignment_results'):
        st.info("No alignment results yet. Run alignment in the 'Run Alignment' tab.")
        return

    results = st.session_state.alignment_results

    # Summary table
    summary_data = []
    for sample, result in results.items():
        if result['status'] == 'success' and result['stats']:
            stats = result['stats']
            summary_data.append({
                'Sample': sample,
                'Total Reads': stats['total_reads'],
                'Aligned': stats['aligned_reads'],
                'Not Aligned': stats['not_aligned'],
                'Suppressed (-m)': stats.get('suppressed', 0),
                'Alignment Rate': f"{stats['alignment_rate']:.1f}%",
                'Status': '✅'
            })
        else:
            summary_data.append({
                'Sample': sample,
                'Total Reads': 'N/A',
                'Aligned': 'N/A',
                'Not Aligned': 'N/A',
                'Suppressed (-m)': 'N/A',
                'Alignment Rate': 'N/A',
                'Status': '❌'
            })

    summary_df = pd.DataFrame(summary_data)
    st.dataframe(summary_df, width="stretch")

    # Check for low alignment rates and show diagnostic info
    low_alignment = [r for r in results.values()
                     if r['status'] == 'success' and r['stats'] and r['stats']['alignment_rate'] < 10]

    if low_alignment:
        with st.expander("⚠️ Low Alignment Rate - Troubleshooting", expanded=True):
            st.warning("One or more samples have very low alignment rates (<10%). Common causes:")

            st.markdown("""
            **1. Wrong Reference Database**
            - Check that your reference matches your organism
            - Verify the reference index was built correctly
            - Current reference: `{}`

            **2. Adapters Not Trimmed**
            - Raw reads contain adapter sequences that prevent alignment
            - Go to Trimming module and verify trimming was successful
            - Check that the correct adapter preset was selected

            **3. Organism Mismatch**
            - Your reads may be from a different organism than the reference
            - Try using a broader reference (e.g., RNAcentral instead of miRBase)

            **4. Reference Too Restrictive**
            - miRBase only contains known miRNAs
            - Try using custom reference or RNAcentral for broader coverage

            **5. Alignment Parameters Too Strict**
            - Try increasing allowed mismatches (currently: {})
            - Disable the -m suppression parameter
            """.format(
                st.session_state.get('reference_index', 'Not set'),
                config.alignment.mismatches
            ))

            # Show a few sample reads for debugging
            st.markdown("**Diagnostic: Check a sample read from your FASTQ**")
            st.info("Compare your read sequences with your reference to verify they should match.")

    # Visualization
    successful = [r for r in results.values() if r['status'] == 'success' and r['stats']]

    if successful:
        import plotly.express as px

        # Alignment rate bar chart
        rates_df = pd.DataFrame([
            {'Sample': r['sample'], 'Alignment Rate': r['stats']['alignment_rate']}
            for r in successful
        ])

        fig = px.bar(
            rates_df,
            x='Sample',
            y='Alignment Rate',
            title='Alignment Rates by Sample',
            color='Alignment Rate',
            color_continuous_scale='RdYlGn'
        )
        fig.update_layout(xaxis_tickangle=-45)
        st.plotly_chart(fig, width="stretch")

        # Stacked bar for read distribution
        dist_data = []
        for r in successful:
            stats = r['stats']
            dist_data.append({
                'Sample': r['sample'],
                'Category': 'Aligned',
                'Reads': stats['aligned_reads']
            })
            dist_data.append({
                'Sample': r['sample'],
                'Category': 'Not Aligned',
                'Reads': stats['not_aligned']
            })
            if stats.get('suppressed', 0) > 0:
                dist_data.append({
                    'Sample': r['sample'],
                    'Category': 'Suppressed',
                    'Reads': stats['suppressed']
                })

        dist_df = pd.DataFrame(dist_data)

        fig2 = px.bar(
            dist_df,
            x='Sample',
            y='Reads',
            color='Category',
            title='Read Distribution by Sample',
            barmode='stack'
        )
        fig2.update_layout(xaxis_tickangle=-45)
        st.plotly_chart(fig2, width="stretch")

    # Download results
    st.divider()
    csv = summary_df.to_csv(index=False)
    st.download_button(
        "📥 Download Alignment Summary",
        csv,
        "alignment_summary.csv",
        "text/csv"
    )
