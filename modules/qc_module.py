"""
Quality Control Module for sRNAtlas
Provides comprehensive QC analysis for small RNA-seq data
"""
import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
import subprocess
import tempfile
import json
import os
from typing import List, Dict, Optional, Tuple
from collections import Counter, defaultdict

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from config import config
from utils.plotting import (
    plot_read_length_distribution,
    plot_quality_distribution,
    plot_rna_type_distribution
)
from utils.caching import hash_file, hash_dataframe
from utils.qc_scorecard import (
    QCStatus,
    QCMetric,
    SampleScorecard,
    QC_THRESHOLDS,
    evaluate_metric,
    classify_size_distribution,
    classify_plant_sirna,
    generate_sample_scorecard,
    scorecard_to_dataframe,
    get_status_emoji,
    get_status_color
)


def analyze_fastq_basic(fastq_file: Path, max_reads: int = 100000) -> Dict:
    """
    Analyze FASTQ file with basic Python (no external tools needed)
    This is a wrapper that handles caching via file hash.

    Args:
        fastq_file: Path to FASTQ file
        max_reads: Maximum reads to sample for analysis

    Returns:
        Dictionary with QC metrics
    """
    file_hash = hash_file(fastq_file)
    return _analyze_fastq_basic_cached(file_hash, str(fastq_file), max_reads)


@st.cache_data(show_spinner=False, ttl=3600)
def _analyze_fastq_basic_cached(_file_hash: str, fastq_file_path: str, max_reads: int = 100000) -> Dict:
    """Cached implementation of FASTQ analysis"""
    import gzip

    fastq_file = Path(fastq_file_path)
    lengths = []
    qualities = []
    gc_contents = []
    read_count = 0
    total_bases = 0

    # Open file (handle gzip)
    opener = gzip.open if str(fastq_file).endswith('.gz') else open
    mode = 'rt' if str(fastq_file).endswith('.gz') else 'r'

    try:
        with opener(fastq_file, mode) as f:
            while read_count < max_reads:
                # Read 4 lines (one FASTQ record)
                header = f.readline()
                if not header:
                    break

                sequence = f.readline().strip()
                plus = f.readline()
                quality = f.readline().strip()

                if not sequence:
                    break

                read_count += 1
                seq_len = len(sequence)
                lengths.append(seq_len)
                total_bases += seq_len

                # Calculate GC content
                gc = (sequence.upper().count('G') + sequence.upper().count('C')) / seq_len * 100
                gc_contents.append(gc)

                # Calculate mean quality score for this read
                if quality:
                    qual_scores = [ord(c) - 33 for c in quality]  # Phred+33
                    qualities.append(np.mean(qual_scores))

        # Calculate statistics
        lengths = np.array(lengths)
        qualities = np.array(qualities)
        gc_contents = np.array(gc_contents)

        # Length distribution
        length_counts = Counter(lengths)

        return {
            'status': 'success',
            'total_reads': read_count,
            'total_bases': total_bases,
            'mean_length': float(np.mean(lengths)),
            'median_length': float(np.median(lengths)),
            'min_length': int(np.min(lengths)),
            'max_length': int(np.max(lengths)),
            'std_length': float(np.std(lengths)),
            'mean_quality': float(np.mean(qualities)) if len(qualities) > 0 else 0,
            'mean_gc': float(np.mean(gc_contents)),
            'length_distribution': dict(length_counts),
            'lengths': lengths.tolist()[:10000],  # Keep sample for plotting
            'qualities': qualities.tolist()[:10000],
        }

    except Exception as e:
        return {
            'status': 'error',
            'error': str(e)
        }


def run_fastqc(fastq_files: List[Path], output_dir: Path, threads: int = 4) -> Dict:
    """
    Run FastQC on FASTQ files

    Args:
        fastq_files: List of FASTQ file paths
        output_dir: Directory to save FastQC outputs
        threads: Number of threads

    Returns:
        Dictionary with FastQC results
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    results = {}

    for fastq in fastq_files:
        sample_name = fastq.stem.replace('.fastq', '').replace('.fq', '').replace('.gz', '')

        cmd = [
            'fastqc',
            '--outdir', str(output_dir),
            '--threads', str(threads),
            '--extract',
            str(fastq)
        ]

        try:
            subprocess.run(cmd, check=True, capture_output=True)

            # Parse FastQC data
            fastqc_dir = output_dir / f"{sample_name}_fastqc"
            fastqc_data = parse_fastqc_data(fastqc_dir / 'fastqc_data.txt')

            results[sample_name] = {
                'status': 'success',
                'data': fastqc_data
            }

        except subprocess.CalledProcessError as e:
            results[sample_name] = {
                'status': 'error',
                'error': str(e)
            }
        except FileNotFoundError:
            results[sample_name] = {
                'status': 'error',
                'error': 'FastQC not found. Please install FastQC.'
            }

    return results


def parse_fastqc_data(data_file: Path) -> Dict:
    """
    Parse FastQC data file

    Args:
        data_file: Path to fastqc_data.txt

    Returns:
        Dictionary with parsed QC metrics
    """
    if not data_file.exists():
        return {}

    data = {
        'basic_stats': {},
        'per_base_quality': [],
        'per_sequence_quality': [],
        'per_base_content': [],
        'gc_content': [],
        'sequence_length': [],
        'duplication': [],
        'overrepresented': [],
        'adapter_content': [],
        'module_status': {}
    }

    current_module = None
    current_data = []

    with open(data_file, 'r') as f:
        for line in f:
            line = line.strip()

            if line.startswith('>>'):
                # Module header
                if current_module and current_data:
                    data[current_module] = current_data
                    current_data = []

                parts = line[2:].split('\t')
                module_name = parts[0].lower().replace(' ', '_')
                status = parts[1] if len(parts) > 1 else 'unknown'

                current_module = module_name
                data['module_status'][module_name] = status

            elif line.startswith('#'):
                continue

            elif line and current_module:
                current_data.append(line.split('\t'))

        # Save last module
        if current_module and current_data:
            data[current_module] = current_data

    # Parse basic statistics
    if 'basic_statistics' in data:
        for row in data['basic_statistics']:
            if len(row) >= 2:
                data['basic_stats'][row[0]] = row[1]

    return data


def analyze_size_distribution(
    lengths: List[int],
    size_ranges: Optional[Dict[str, Tuple[int, int]]] = None
) -> Dict:
    """
    Analyze read length distribution for small RNA types

    Args:
        lengths: List of read lengths
        size_ranges: Expected size ranges for RNA types

    Returns:
        Dictionary with size distribution analysis
    """
    if size_ranges is None:
        size_ranges = config.qc.size_ranges

    length_counts = Counter(lengths)
    total_reads = len(lengths)

    # Categorize reads by size
    categorized = {}
    for rna_type, (min_len, max_len) in size_ranges.items():
        count = sum(c for l, c in length_counts.items() if min_len <= l <= max_len)
        pct = 100 * count / total_reads if total_reads > 0 else 0
        categorized[rna_type] = {
            'count': count,
            'percentage': pct,
            'size_range': f"{min_len}-{max_len} nt"
        }

    # Peak detection
    sorted_lengths = sorted(length_counts.items())
    peaks = []
    for i, (length, count) in enumerate(sorted_lengths):
        if i > 0 and i < len(sorted_lengths) - 1:
            prev_count = sorted_lengths[i-1][1]
            next_count = sorted_lengths[i+1][1]
            if count > prev_count and count > next_count:
                peaks.append({
                    'length': length,
                    'count': count,
                    'percentage': 100 * count / total_reads
                })

    # Sort peaks by count
    peaks.sort(key=lambda x: x['count'], reverse=True)

    return {
        'length_counts': dict(length_counts),
        'total_reads': total_reads,
        'mean_length': np.mean(lengths),
        'median_length': np.median(lengths),
        'min_length': min(lengths),
        'max_length': max(lengths),
        'categorized': categorized,
        'peaks': peaks[:5]  # Top 5 peaks
    }


def render_rna_type_reference():
    """Render comprehensive small RNA type reference guide"""
    st.markdown("""
    ### Small RNA Types Reference

    Small RNAs are non-coding RNA molecules typically <200 nucleotides. Below is a comprehensive guide to the major types detected in small RNA-seq experiments.
    """)

    # Core small RNA types table
    st.markdown("#### Core Small RNA Types")
    core_types = pd.DataFrame([
        {'Type': 'miRNA', 'Full Name': 'MicroRNA', 'Size (nt)': '18-25',
         'Description': 'Gene expression regulators via mRNA targeting. Show 5\' U bias.'},
        {'Type': 'siRNA', 'Full Name': 'Small interfering RNA', 'Size (nt)': '20-24',
         'Description': 'RNAi pathway, requires perfect complementarity to targets.'},
        {'Type': 'piRNA', 'Full Name': 'PIWI-interacting RNA', 'Size (nt)': '24-32',
         'Description': 'Transposon silencing in germline. Show 5\' U or 10A bias.'},
        {'Type': 'tRF/tsRNA', 'Full Name': 'tRNA-derived fragment', 'Size (nt)': '14-40',
         'Description': 'Stress-induced fragments from tRNAs. Gene regulation roles.'},
        {'Type': 'rsRF', 'Full Name': 'rRNA-derived fragment', 'Size (nt)': '15-40',
         'Description': 'Fragments from ribosomal RNA processing.'},
    ])
    st.dataframe(core_types, width="stretch", hide_index=True)

    # Other small RNA types
    st.markdown("#### Other Small RNA Types")
    other_types = pd.DataFrame([
        {'Type': 'snoRNA', 'Full Name': 'Small nucleolar RNA', 'Size (nt)': '60-300',
         'Description': 'Guide chemical modifications of rRNA and snRNA.'},
        {'Type': 'snRNA', 'Full Name': 'Small nuclear RNA', 'Size (nt)': '100-300',
         'Description': 'Core splicing machinery components (U1, U2, U4, U5, U6).'},
        {'Type': 'Y RNA', 'Full Name': 'Y RNA', 'Size (nt)': '80-120',
         'Description': 'DNA replication initiation, RNA quality control.'},
        {'Type': 'vtRNA', 'Full Name': 'Vault RNA', 'Size (nt)': '80-150',
         'Description': 'Associated with vault ribonucleoprotein complex.'},
        {'Type': '7SL/SRP', 'Full Name': 'Signal recognition particle RNA', 'Size (nt)': '~300',
         'Description': 'Protein targeting to endoplasmic reticulum.'},
        {'Type': '7SK', 'Full Name': '7SK RNA', 'Size (nt)': '~330',
         'Description': 'Transcription regulation, P-TEFb sequestration.'},
    ])
    st.dataframe(other_types, width="stretch", hide_index=True)

    # Plant-specific types
    st.markdown("#### Plant-Specific Small RNAs")
    plant_types = pd.DataFrame([
        {'Type': 'tasiRNA', 'Full Name': 'Trans-acting siRNA', 'Size (nt)': '21',
         'Description': 'Phased production from TAS loci, target developmental genes.'},
        {'Type': 'phasiRNA', 'Full Name': 'Phased siRNA', 'Size (nt)': '21 or 24',
         'Description': 'Phased pattern, enriched in reproductive tissues.'},
        {'Type': 'natsiRNA', 'Full Name': 'Natural antisense siRNA', 'Size (nt)': '21-24',
         'Description': 'Derived from overlapping sense/antisense transcripts.'},
        {'Type': 'hc-siRNA', 'Full Name': 'Heterochromatic siRNA', 'Size (nt)': '24',
         'Description': 'DNA methylation and transposon silencing via RdDM pathway.'},
    ])
    st.dataframe(plant_types, width="stretch", hide_index=True)

    # Size interpretation guide
    st.markdown("""
    #### Interpreting Size Distribution Peaks

    | Peak Position | Likely RNA Type | Interpretation |
    |---------------|-----------------|----------------|
    | 21-23 nt | miRNA | Good miRNA library enrichment |
    | 24 nt | siRNA/hc-siRNA | Plant siRNA or heterochromatic siRNA |
    | 26-32 nt | piRNA | PIWI pathway active (animal germline) |
    | 28-35 nt | tRF | tRNA fragments, may indicate stress |
    | Broad distribution | Degradation | Check RNA quality, may need size selection |
    """)


def check_rrna_contamination(
    annotations: pd.DataFrame,
    count_matrix: pd.DataFrame
) -> Dict:
    """
    Check for rRNA contamination in the data (cached wrapper)

    Args:
        annotations: DataFrame with RNA annotations
        count_matrix: Count matrix

    Returns:
        Dictionary with rRNA contamination analysis
    """
    ann_hash = hash_dataframe(annotations)
    count_hash = hash_dataframe(count_matrix)
    return _check_rrna_contamination_cached(ann_hash, count_hash, annotations, count_matrix)


@st.cache_data(show_spinner=False, ttl=3600)
def _check_rrna_contamination_cached(
    _ann_hash: str,
    _count_hash: str,
    annotations: pd.DataFrame,
    count_matrix: pd.DataFrame
) -> Dict:
    """Cached implementation of rRNA contamination check"""
    if 'RNA_category' not in annotations.columns:
        return {'error': 'No RNA_category column in annotations'}

    # Get rRNA reads
    rrna_mask = annotations['RNA_category'] == 'rRNA'
    rrna_ids = annotations[rrna_mask].index

    # Calculate rRNA reads per sample
    sample_cols = [c for c in count_matrix.columns if c not in annotations.columns]

    rrna_reads = {}
    total_reads = {}

    for sample in sample_cols:
        sample_total = count_matrix[sample].sum()
        rrna_total = count_matrix.loc[count_matrix.index.isin(rrna_ids), sample].sum()

        total_reads[sample] = sample_total
        rrna_reads[sample] = rrna_total

    # Calculate percentages
    rrna_pct = {s: 100 * rrna_reads[s] / total_reads[s] if total_reads[s] > 0 else 0
                for s in sample_cols}

    return {
        'rrna_reads': rrna_reads,
        'total_reads': total_reads,
        'rrna_percentage': rrna_pct,
        'mean_rrna_pct': np.mean(list(rrna_pct.values())),
        'samples_high_rrna': [s for s, p in rrna_pct.items() if p > 50]
    }


def render_qc_page():
    """Render the QC analysis page - single page with all results"""
    st.header("📊 Quality Control")

    # Check if we have data
    if not st.session_state.get('project_name'):
        st.warning("Please set up a project first in the Project Setup page.")
        return

    # === INPUT SECTION ===
    render_qc_input_section()

    # === RESULTS SECTION ===
    if st.session_state.get('qc_results'):
        st.divider()
        render_qc_results_section()


def render_qc_input_section():
    """Render input file selection and run analysis button"""
    st.subheader("📁 Select Files for QC Analysis")

    # Check for files from Project module
    project_files = st.session_state.get('project_fastq_files', [])
    project_names = st.session_state.get('project_fastq_names', [])
    has_project_files = project_files and len(project_files) > 0

    # Input source selection
    if has_project_files:
        st.success(f"✅ **{len(project_files)} files available from Project**")
        input_source = st.radio(
            "Input Source",
            options=["Use files from Project", "Upload new files"],
            index=0,
            horizontal=True,
            key="qc_input_source_main"
        )
    else:
        input_source = "Upload new files"
        st.info("📁 Upload FASTQ files here, or load them first in the **Project** module.")

    # Show files or upload widget
    files_to_analyze = []
    file_names = []

    if input_source == "Use files from Project" and has_project_files:
        # Show project files in a compact list
        with st.expander(f"📄 Files to analyze ({len(project_files)})", expanded=False):
            for name in project_names:
                st.caption(f"• {name}")

        files_to_analyze = [p for p in project_files if isinstance(p, Path) and p.exists()]
        file_names = project_names

    else:
        # File uploader
        uploaded_files = st.file_uploader(
            "Upload FASTQ files",
            type=None,
            accept_multiple_files=True,
            help="Upload FASTQ files (.fastq, .fq, .fastq.gz, .fq.gz)",
            key="qc_fastq_upload_main"
        )

        if uploaded_files:
            valid_extensions = ('.fastq', '.fq', '.fastq.gz', '.fq.gz')
            valid_files = [f for f in uploaded_files if f.name.lower().endswith(valid_extensions)]

            if valid_files:
                st.success(f"Selected {len(valid_files)} files")
                files_to_analyze = valid_files
                file_names = [f.name for f in valid_files]

    # Run Analysis button
    if files_to_analyze:
        st.divider()

        col1, col2 = st.columns([3, 1])
        with col1:
            st.markdown("**Ready to analyze.** This will run: Read Statistics, Size Distribution, and Contamination Check.")
        with col2:
            run_qc = st.button("🔍 Run QC Analysis", type="primary", key="run_qc_main")

        if run_qc:
            run_comprehensive_qc(files_to_analyze, file_names, input_source == "Use files from Project")


def run_comprehensive_qc(files, file_names, files_are_paths=False):
    """Run all QC analyses on the provided files"""
    progress_bar = st.progress(0)
    status_text = st.empty()

    all_results = {}
    all_lengths = []
    contamination_results = {}

    # Save uploaded files if needed
    if not files_are_paths:
        temp_dir = Path(tempfile.mkdtemp())
        saved_files = []
        for uf in files:
            file_path = temp_dir / uf.name
            with open(file_path, 'wb') as f:
                f.write(uf.getbuffer())
            saved_files.append(file_path)
        files = saved_files

    total_steps = len(files) * 2  # Basic analysis + contamination check
    current_step = 0

    # Step 1: Basic FASTQ analysis
    for i, (file_path, name) in enumerate(zip(files, file_names)):
        status_text.text(f"Analyzing {name}... (basic stats)")
        current_step += 1
        progress_bar.progress(current_step / total_steps)

        result = analyze_fastq_basic(file_path)
        result['filename'] = name
        all_results[name] = result

        if result['status'] == 'success':
            all_lengths.extend(result.get('lengths', []))

    # Step 2: Contamination check
    for i, (file_path, name) in enumerate(zip(files, file_names)):
        status_text.text(f"Checking contamination in {name}...")
        current_step += 1
        progress_bar.progress(current_step / total_steps)

        contam_result = analyze_fastq_contamination(file_path)
        contamination_results[name] = contam_result

    progress_bar.progress(1.0)
    status_text.text("✅ Analysis complete!")

    # Store all results
    st.session_state.qc_results = {
        'files': file_names,
        'file_paths': files,
        'analysis': all_results,
        'all_lengths': all_lengths,
        'contamination': contamination_results
    }
    st.session_state.fastq_qc_done = True

    st.success(f"✅ Analyzed {len(files)} files!")
    st.rerun()


def render_qc_scorecard():
    """Render QC Scorecard with traffic-light status indicators"""
    qc_results = st.session_state.get('qc_results')
    if not qc_results:
        return

    analysis = qc_results.get('analysis', {})
    contamination = qc_results.get('contamination', {})

    # Check if we're analyzing plant data
    is_plant = st.session_state.get('is_plant_data', False)

    # Generate scorecards for all samples
    scorecards = []

    for sample_name, result in analysis.items():
        if result.get('status') != 'success':
            continue

        # Get contamination results for this sample
        contam_result = contamination.get(sample_name, {})

        # Build QC results dict for scorecard
        qc_data = {
            'total_reads': result.get('total_reads', 0),
            'mean_quality': result.get('mean_quality', 0),
            'lengths': result.get('lengths', []),
            'length_distribution': result.get('length_distribution', {})
        }

        # Generate the scorecard
        # Use 'or {}' pattern to handle both missing keys and None values
        alignment_data = (st.session_state.get('alignment_results') or {}).get(sample_name)

        scorecard = generate_sample_scorecard(
            sample_name=sample_name,
            qc_results=qc_data,
            trim_results=None,  # Not available at this stage
            alignment_results=alignment_data,
            contamination_results=contam_result if contam_result.get('status') == 'success' else None,
            is_plant=is_plant
        )
        scorecards.append(scorecard)

    if not scorecards:
        return

    # === Render Scorecard UI ===
    st.subheader("🎯 QC Scorecard Summary")

    # Overall status summary
    status_counts = {
        QCStatus.OK: sum(1 for sc in scorecards if sc.overall_status == QCStatus.OK),
        QCStatus.WARNING: sum(1 for sc in scorecards if sc.overall_status == QCStatus.WARNING),
        QCStatus.CRITICAL: sum(1 for sc in scorecards if sc.overall_status == QCStatus.CRITICAL)
    }

    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Total Samples", len(scorecards))
    with col2:
        st.metric("✅ Passed", status_counts[QCStatus.OK],
                 delta_color="normal" if status_counts[QCStatus.OK] == len(scorecards) else "off")
    with col3:
        st.metric("⚠️ Warnings", status_counts[QCStatus.WARNING],
                 delta_color="off" if status_counts[QCStatus.WARNING] > 0 else "normal")
    with col4:
        st.metric("❌ Critical", status_counts[QCStatus.CRITICAL],
                 delta_color="inverse" if status_counts[QCStatus.CRITICAL] > 0 else "normal")

    # Detailed scorecard table
    st.markdown("#### Sample Status Overview")

    # Create styled table
    scorecard_data = []
    for sc in scorecards:
        row = {
            'Sample': sc.sample_name,
            'Status': f"{get_status_emoji(sc.overall_status)} {sc.overall_status.value.upper()}"
        }

        # Add key metrics with status
        for metric in sc.metrics:
            metric_emoji = get_status_emoji(metric.status)

            if metric.name == 'total_reads':
                row['Reads'] = f"{metric_emoji} {metric.value:,.0f}"
            elif metric.name == 'mean_quality':
                row['Quality'] = f"{metric_emoji} {metric.value:.1f}"
            elif metric.name == 'adapter_pct':
                row['Adapters'] = f"{metric_emoji} {metric.value:.1f}%"
            elif metric.name == 'size_pattern':
                row['Pattern'] = f"{metric_emoji} {metric.unit}"
            elif metric.name == 'alignment_rate':
                row['Aligned'] = f"{metric_emoji} {metric.value:.1f}%"

        scorecard_data.append(row)

    if scorecard_data:
        scorecard_df = pd.DataFrame(scorecard_data)
        st.dataframe(scorecard_df, width="stretch", hide_index=True)

    # Expandable detailed view per sample
    st.markdown("#### Detailed Sample Scorecards")

    for sc in scorecards:
        status_emoji = get_status_emoji(sc.overall_status)
        status_color = get_status_color(sc.overall_status)

        with st.expander(f"{status_emoji} **{sc.sample_name}** - {sc.summary}"):
            # Metrics grid
            st.markdown("**Metrics:**")

            # Display metrics in columns
            metric_cols = st.columns(3)

            for i, metric in enumerate(sc.metrics):
                col_idx = i % 3
                with metric_cols[col_idx]:
                    metric_emoji = get_status_emoji(metric.status)

                    if metric.unit in ['%', 'reads', 'Phred']:
                        # Numeric metric
                        value_str = f"{metric.value:,.1f} {metric.unit}" if metric.unit else f"{metric.value:,.1f}"
                    else:
                        # Pattern/categorical metric
                        value_str = metric.unit

                    st.markdown(f"""
                    **{metric_emoji} {metric.name.replace('_', ' ').title()}**
                    {value_str}
                    *{metric.description}*
                    """)

            # Thresholds info
            if any(m.threshold_warning for m in sc.metrics):
                with st.popover("📊 Threshold Details"):
                    st.markdown("**QC Thresholds Used:**")
                    for metric in sc.metrics:
                        if metric.threshold_warning or metric.threshold_critical:
                            direction = QC_THRESHOLDS.get(metric.name, {}).get('direction', 'min')
                            if direction == 'min':
                                st.markdown(f"- **{metric.name}**: ⚠️ < {metric.threshold_warning}, ❌ < {metric.threshold_critical}")
                            else:
                                st.markdown(f"- **{metric.name}**: ⚠️ > {metric.threshold_warning}, ❌ > {metric.threshold_critical}")

            # Recommendations
            if sc.recommendations:
                st.markdown("**Recommendations:**")
                for rec in sc.recommendations:
                    if "No issues" in rec:
                        st.success(f"✓ {rec}")
                    else:
                        st.warning(f"→ {rec}")

    # Option to view/customize thresholds
    with st.expander("⚙️ QC Threshold Settings"):
        st.markdown("""
        **Current QC Thresholds:**

        These thresholds determine the traffic-light status for each metric.
        """)

        threshold_data = []
        for metric_name, thresholds in QC_THRESHOLDS.items():
            threshold_data.append({
                'Metric': metric_name.replace('_', ' ').title(),
                'Warning': thresholds.get('warning', '-'),
                'Critical': thresholds.get('critical', '-'),
                'Direction': '↑ Higher is better' if thresholds.get('direction') == 'min' else '↓ Lower is better'
            })

        st.dataframe(pd.DataFrame(threshold_data), width="stretch", hide_index=True)

        # Plant data toggle
        st.checkbox(
            "Plant Data Mode (enables 21nt/24nt classification)",
            key="is_plant_data",
            help="Enable plant-specific siRNA classification (21nt miRNA/ta-siRNA vs 24nt hc-siRNA)"
        )

    # Download scorecard
    if scorecards:
        scorecard_df_full = scorecard_to_dataframe(scorecards)
        csv = scorecard_df_full.to_csv(index=False)
        st.download_button(
            "📥 Download Full Scorecard",
            csv,
            "qc_scorecard.csv",
            "text/csv"
        )

    st.divider()


def render_qc_results_section():
    """Render all QC results on a single page"""
    qc_results = st.session_state.get('qc_results')
    if not qc_results:
        return

    analysis = qc_results.get('analysis', {})
    contamination = qc_results.get('contamination', {})
    all_lengths = qc_results.get('all_lengths', [])

    # === SECTION 0: QC Scorecard (Traffic Light Summary) ===
    render_qc_scorecard()

    # === SECTION 0.5: Multi-Sample QC Overlay (collapsible) ===
    with st.expander("🔍 Multi-Sample QC Overlay & Outlier Detection", expanded=False):
        render_multi_sample_qc_overlay()

    # === SECTION 1: Summary Statistics ===
    st.subheader("📈 Read Statistics Summary")

    # Summary table
    summary_data = []
    for name, result in analysis.items():
        if result.get('status') == 'success':
            summary_data.append({
                'Sample': name,
                'Total Reads': f"{result.get('total_reads', 0):,}",
                'Mean Length': f"{result.get('mean_length', 0):.1f}",
                'Mean Quality': f"{result.get('mean_quality', 0):.1f}",
                'GC %': f"{result.get('mean_gc', 0):.1f}",
                'Status': '✅'
            })
        else:
            summary_data.append({
                'Sample': name,
                'Total Reads': 'Error',
                'Mean Length': '-',
                'Mean Quality': '-',
                'GC %': '-',
                'Status': '❌'
            })

    if summary_data:
        st.dataframe(pd.DataFrame(summary_data), width="stretch")

    # Quick metrics
    successful = [r for r in analysis.values() if r.get('status') == 'success']
    if successful:
        col1, col2, col3, col4 = st.columns(4)
        total_reads = sum(r.get('total_reads', 0) for r in successful)
        avg_quality = np.mean([r.get('mean_quality', 0) for r in successful])
        avg_gc = np.mean([r.get('mean_gc', 0) for r in successful])
        avg_length = np.mean([r.get('mean_length', 0) for r in successful])

        with col1:
            st.metric("Total Reads", f"{total_reads:,}")
        with col2:
            st.metric("Avg Quality", f"{avg_quality:.1f}")
        with col3:
            st.metric("Avg GC%", f"{avg_gc:.1f}%")
        with col4:
            st.metric("Avg Length", f"{avg_length:.1f} nt")

    # === SECTION 2: Size Distribution ===
    st.divider()
    st.subheader("🔬 Size Distribution")

    if all_lengths:
        import plotly.graph_objects as go

        # Create histogram
        length_counts = Counter(all_lengths)
        lengths = sorted(length_counts.keys())
        counts = [length_counts[l] for l in lengths]

        fig = go.Figure()

        # Add histogram bars
        fig.add_trace(go.Bar(
            x=lengths,
            y=counts,
            name='Read counts',
            marker_color='steelblue'
        ))

        # Add RNA type shaded regions
        rna_ranges = [
            ('miRNA', 18, 25, 'rgba(255, 99, 132, 0.2)'),
            ('siRNA', 20, 24, 'rgba(54, 162, 235, 0.2)'),
            ('piRNA', 24, 32, 'rgba(255, 206, 86, 0.2)'),
            ('tRF', 28, 36, 'rgba(75, 192, 192, 0.2)'),
        ]

        for rna_type, start, end, color in rna_ranges:
            fig.add_vrect(x0=start, x1=end, fillcolor=color, opacity=0.3, line_width=0)

        # Add legend for RNA types
        for rna_type, start, end, color in rna_ranges:
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode='markers',
                marker=dict(size=10, color=color.replace('0.2', '0.6')),
                name=f"{rna_type} ({start}-{end} nt)",
                showlegend=True
            ))

        fig.update_layout(
            title='Read Length Distribution with RNA Type Annotations',
            xaxis_title='Read Length (nt)',
            yaxis_title='Count',
            xaxis=dict(range=[10, max(60, max(lengths) + 5)]),
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
            height=400
        )

        st.plotly_chart(fig, width="stretch")

        # Peak detection
        if lengths:
            peak_length = max(length_counts, key=length_counts.get)
            st.info(f"📊 **Peak read length:** {peak_length} nt ({length_counts[peak_length]:,} reads)")

            # Interpretation
            if 20 <= peak_length <= 24:
                st.success("✅ Peak consistent with **miRNA** (20-24 nt)")
            elif 26 <= peak_length <= 32:
                st.info("ℹ️ Peak consistent with **piRNA** (26-32 nt)")
            elif peak_length > 40:
                st.warning("⚠️ Peak length suggests reads may not be trimmed or are not small RNA")

    # === SECTION 3: Contamination Check ===
    st.divider()
    st.subheader("⚠️ Contamination Check")

    if contamination:
        # Adapter content
        st.markdown("**Adapter Content:**")
        adapter_data = []
        for name, result in contamination.items():
            if result.get('status') == 'success':
                row = {'Sample': name}
                for adapter_name, pct in result.get('adapter_content', {}).items():
                    row[adapter_name] = f"{pct:.1f}%"
                adapter_data.append(row)

        if adapter_data:
            adapter_df = pd.DataFrame(adapter_data)
            st.dataframe(adapter_df, width="stretch")

            # Warnings for high adapter content
            high_adapter_samples = []
            for name, result in contamination.items():
                if result.get('adapters_detected'):
                    high_adapter_samples.append((name, result['adapters_detected']))

            if high_adapter_samples:
                st.warning("⚠️ **High adapter content detected:**")
                for sample, adapters in high_adapter_samples:
                    st.markdown(f"- **{sample}**: {', '.join(adapters)}")
                st.markdown("*Consider trimming adapters before alignment.*")
            else:
                st.success("✅ No significant adapter contamination detected")

        # Poly-A/T content
        st.markdown("**Poly-A/T Content:**")
        poly_data = []
        for name, result in contamination.items():
            if result.get('status') == 'success':
                poly_data.append({
                    'Sample': name,
                    'Poly-A %': f"{result.get('poly_a_pct', 0):.2f}%",
                    'Poly-T %': f"{result.get('poly_t_pct', 0):.2f}%",
                    'N Content %': f"{result.get('n_content_pct', 0):.2f}%"
                })

        if poly_data:
            st.dataframe(pd.DataFrame(poly_data), width="stretch")

    else:
        st.info("No contamination data available")

    # === Download Section ===
    st.divider()
    col1, col2 = st.columns(2)

    with col1:
        if summary_data:
            csv = pd.DataFrame(summary_data).to_csv(index=False)
            st.download_button(
                "📥 Download QC Summary",
                csv,
                "qc_summary.csv",
                "text/csv"
            )

    with col2:
        if st.button("🔄 Re-run QC Analysis"):
            st.session_state.qc_results = None
            st.rerun()


def render_qc_upload():
    """Legacy function - kept for compatibility"""
    st.subheader("Upload Data for QC Analysis")

    data_type = st.radio(
        "What data do you want to analyze?",
        options=["FASTQ files", "Count matrix with annotations"]
    )

    if data_type == "FASTQ files":
        # Check for files from Project module
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
                horizontal=True,
                key="qc_input_source"
            )

            if input_source == "Use files from Project":
                # Show project files
                for name in project_names:
                    st.caption(f"📄 {name}")

                if st.button("🔍 Run QC Analysis", type="primary", key="qc_project"):
                    # Files are already saved - use directly
                    saved_files = [p for p in project_files if isinstance(p, Path) and p.exists()]

                    if not saved_files:
                        st.error("No valid files found. Please re-upload files in the Project module.")
                        return

                    progress_bar = st.progress(0)
                    status_text = st.empty()

                    all_results = {}
                    all_lengths = []

                    for i, (file_path, name) in enumerate(zip(saved_files, project_names)):
                        status_text.text(f"Analyzing {name}...")
                        progress_bar.progress((i + 1) / len(saved_files))

                        result = analyze_fastq_basic(file_path)
                        result['filename'] = name
                        all_results[name] = result

                        if result['status'] == 'success':
                            all_lengths.extend(result.get('lengths', []))

                    progress_bar.progress(1.0)
                    status_text.text("Analysis complete!")

                    # Store results
                    st.session_state.qc_results = {
                        'files': project_names,
                        'file_paths': saved_files,
                        'analysis': all_results,
                        'all_lengths': all_lengths
                    }
                    st.session_state.fastq_qc_done = True

                    st.success(f"✅ Analyzed {len(saved_files)} files!")
                    st.rerun()
                return

        # Direct file upload
        if not has_project_files:
            st.info("📁 Upload FASTQ files here, or load them first in the **Project** module.")

        uploaded_files = st.file_uploader(
            "Upload FASTQ files",
            type=None,  # Accept all types - filter manually for .fastq.gz support
            accept_multiple_files=True,
            help="Upload FASTQ files (.fastq, .fq, .fastq.gz, .fq.gz)",
            key="qc_fastq_upload"
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
            st.success(f"Uploaded {len(uploaded_files)} files")

            if st.button("🔍 Run QC Analysis", type="primary", key="qc_upload"):
                with st.spinner("Running quality control analysis..."):
                    # Save files temporarily
                    temp_dir = Path(tempfile.mkdtemp())
                    saved_files = []

                    for uf in uploaded_files:
                        file_path = temp_dir / uf.name
                        with open(file_path, 'wb') as f:
                            f.write(uf.getbuffer())
                        saved_files.append(file_path)

                    # Run analysis
                    st.session_state.qc_results = {
                        'files': [f.name for f in saved_files],
                        'file_paths': saved_files,
                        'temp_dir': str(temp_dir)
                    }

                    st.success("QC analysis complete!")
                    st.rerun()

    else:  # Count matrix
        col1, col2 = st.columns(2)

        with col1:
            count_file = st.file_uploader(
                "Upload count matrix (CSV)",
                type=['csv', 'tsv', 'txt']
            )

        with col2:
            annot_file = st.file_uploader(
                "Upload annotations (CSV) - optional",
                type=['csv', 'tsv', 'txt']
            )

        if count_file:
            try:
                # Detect delimiter
                content = count_file.getvalue().decode('utf-8')
                sep = '\t' if '\t' in content.split('\n')[0] else ','

                count_file.seek(0)
                counts = pd.read_csv(count_file, sep=sep, index_col=0)

                st.success(f"Loaded count matrix: {counts.shape[0]} features × {counts.shape[1]} samples")

                # Store in session state
                st.session_state.count_matrix = counts

                if annot_file:
                    annot_file.seek(0)
                    annotations = pd.read_csv(annot_file, sep=sep, index_col=0)
                    st.session_state.annotations = annotations
                    st.success(f"Loaded annotations: {len(annotations)} entries")

            except Exception as e:
                st.error(f"Error loading file: {e}")


def render_read_stats():
    """Render read statistics section"""
    st.subheader("📈 Read Statistics")

    # Check for FASTQ QC results first
    qc_results = st.session_state.get('qc_results')
    if qc_results and 'analysis' in qc_results:
        analysis = qc_results['analysis']

        st.success(f"✅ FASTQ QC results for {len(analysis)} files")

        # Summary table
        summary_data = []
        for name, result in analysis.items():
            if result.get('status') == 'success':
                summary_data.append({
                    'Sample': name,
                    'Total Reads': f"{result['total_reads']:,}",
                    'Mean Length': f"{result['mean_length']:.1f}",
                    'Length Range': f"{result['min_length']}-{result['max_length']}",
                    'Mean Quality': f"{result['mean_quality']:.1f}",
                    'GC %': f"{result['mean_gc']:.1f}%"
                })
            else:
                summary_data.append({
                    'Sample': name,
                    'Total Reads': 'Error',
                    'Mean Length': '-',
                    'Length Range': '-',
                    'Mean Quality': '-',
                    'GC %': result.get('error', 'Unknown error')
                })

        if summary_data:
            st.dataframe(pd.DataFrame(summary_data), width="stretch")

        # Aggregated metrics
        successful = [r for r in analysis.values() if r.get('status') == 'success']
        if successful:
            col1, col2, col3, col4 = st.columns(4)

            total_reads = sum(r['total_reads'] for r in successful)
            avg_length = np.mean([r['mean_length'] for r in successful])
            avg_quality = np.mean([r['mean_quality'] for r in successful])
            avg_gc = np.mean([r['mean_gc'] for r in successful])

            with col1:
                st.metric("Total Reads", f"{total_reads:,}")
            with col2:
                st.metric("Avg Read Length", f"{avg_length:.1f} nt")
            with col3:
                st.metric("Avg Quality Score", f"{avg_quality:.1f}")
            with col4:
                st.metric("Avg GC Content", f"{avg_gc:.1f}%")

            # Length distribution plot
            st.subheader("Read Length Distribution")

            all_lengths = qc_results.get('all_lengths', [])
            if all_lengths:
                import plotly.express as px
                fig = px.histogram(
                    x=all_lengths,
                    nbins=50,
                    title="Read Length Distribution (All Samples)",
                    labels={'x': 'Read Length (nt)', 'y': 'Count'}
                )
                fig.update_layout(showlegend=False)
                st.plotly_chart(fig, width="stretch")

        return

    # Fall back to count matrix analysis
    if st.session_state.get('count_matrix') is not None:
        counts = st.session_state.count_matrix

        # Identify sample columns (numeric only)
        sample_cols = counts.select_dtypes(include=[np.number]).columns.tolist()

        if not sample_cols:
            st.warning("No numeric columns found in count matrix")
            return

        # Calculate statistics
        stats_df = pd.DataFrame({
            'Total Reads': counts[sample_cols].sum(),
            'Detected Features': (counts[sample_cols] > 0).sum(),
            'Mean Counts': counts[sample_cols].mean(),
            'Median Counts': counts[sample_cols].median(),
            'Max Counts': counts[sample_cols].max()
        })

        # Display summary metrics
        col1, col2, col3, col4 = st.columns(4)

        with col1:
            st.metric("Total Samples", len(sample_cols))

        with col2:
            st.metric("Total Features", len(counts))

        with col3:
            st.metric("Mean Reads/Sample", f"{stats_df['Total Reads'].mean():,.0f}")

        with col4:
            st.metric("Features Detected (>0 in any)", (counts[sample_cols].sum(axis=1) > 0).sum())

        # Show statistics table
        st.dataframe(stats_df, width="stretch")

        # Library size distribution
        st.subheader("Library Size Distribution")

        import plotly.express as px
        fig = px.bar(
            x=sample_cols,
            y=counts[sample_cols].sum().values,
            labels={'x': 'Sample', 'y': 'Total Reads'},
            title='Library Sizes'
        )
        fig.update_layout(xaxis_tickangle=-45)
        st.plotly_chart(fig, width="stretch")

        # Download statistics
        csv = stats_df.to_csv()
        st.download_button(
            "📥 Download Statistics",
            csv,
            "qc_statistics.csv",
            "text/csv"
        )

    else:
        st.info("Please upload data in the 'Upload Data' tab first.")


def render_size_distribution():
    """Render size distribution analysis"""
    st.subheader("🔬 Read Size Distribution")

    # Check for FASTQ QC results first
    qc_results = st.session_state.get('qc_results')
    if qc_results and 'analysis' in qc_results:
        all_lengths = qc_results.get('all_lengths', [])

        if all_lengths:
            st.success(f"✅ Size distribution from FASTQ QC ({len(all_lengths):,} reads sampled)")

            lengths = np.array(all_lengths)

            # Metrics
            col1, col2, col3, col4 = st.columns(4)

            with col1:
                st.metric("Total Sampled", f"{len(lengths):,}")
            with col2:
                st.metric("Mean Length", f"{np.mean(lengths):.1f} nt")
            with col3:
                st.metric("Median Length", f"{np.median(lengths):.0f} nt")
            with col4:
                st.metric("Range", f"{np.min(lengths)}-{np.max(lengths)} nt")

            # Size distribution plot with RNA type annotations
            import plotly.graph_objects as go

            fig = go.Figure()

            # Histogram
            fig.add_trace(go.Histogram(
                x=lengths,
                nbinsx=50,
                name="Read Lengths",
                marker_color='#3498db'
            ))

            # Add RNA type reference regions (non-overlapping display)
            rna_ranges = [
                ('miRNA', 18, 25, '#1f77b4'),
                ('siRNA', 20, 24, '#ff7f0e'),
                ('piRNA', 24, 32, '#2ca02c'),
                ('tRF', 28, 40, '#d62728'),
                ('snoRNA', 60, 90, '#9467bd'),
                ('snRNA', 100, 150, '#8c564b'),
            ]

            # Add shaded regions without text annotations
            for rna_type, start, end, color in rna_ranges:
                fig.add_vrect(
                    x0=start, x1=end,
                    fillcolor=color, opacity=0.15,
                    line_width=1,
                    line_color=color,
                )

            # Add a legend using invisible scatter traces
            for rna_type, start, end, color in rna_ranges:
                fig.add_trace(go.Scatter(
                    x=[None], y=[None],
                    mode='markers',
                    marker=dict(size=10, color=color),
                    name=f"{rna_type} ({start}-{end} nt)",
                    showlegend=True
                ))

            fig.update_layout(
                title="Read Length Distribution with RNA Type Ranges",
                xaxis_title="Read Length (nt)",
                yaxis_title="Count",
                showlegend=True,
                legend=dict(
                    orientation="h",
                    yanchor="bottom",
                    y=1.02,
                    xanchor="center",
                    x=0.5,
                    font=dict(size=10)
                ),
                height=450
            )

            st.plotly_chart(fig, width="stretch")

            # Categorize reads by RNA type
            st.subheader("Reads by Expected RNA Type")

            categories = {
                'miRNA (18-25 nt)': np.sum((lengths >= 18) & (lengths <= 25)),
                'siRNA (20-24 nt)': np.sum((lengths >= 20) & (lengths <= 24)),
                'piRNA (24-32 nt)': np.sum((lengths >= 24) & (lengths <= 32)),
                'tRF/tsRNA (14-40 nt)': np.sum((lengths >= 14) & (lengths <= 40)),
                'rsRF (15-40 nt)': np.sum((lengths >= 15) & (lengths <= 40)),
                'Short (<14 nt)': np.sum(lengths < 14),
                'Long (>50 nt)': np.sum(lengths > 50),
            }

            cat_df = pd.DataFrame([
                {'RNA Type': k, 'Count': v, 'Percentage': f"{100*v/len(lengths):.1f}%"}
                for k, v in categories.items()
            ])
            cat_df = cat_df.sort_values('Count', ascending=False)

            st.dataframe(cat_df, width="stretch", hide_index=True)

            # RNA Type Reference Guide
            with st.expander("📖 Small RNA Type Reference Guide"):
                render_rna_type_reference()

            # Peak detection
            st.subheader("Detected Peaks")
            length_counts = Counter(lengths)
            sorted_counts = sorted(length_counts.items(), key=lambda x: x[1], reverse=True)
            top_peaks = sorted_counts[:5]

            peaks_df = pd.DataFrame([
                {'Length (nt)': length, 'Count': count, 'Percentage': f"{100*count/len(lengths):.2f}%"}
                for length, count in top_peaks
            ])
            st.dataframe(peaks_df, width="stretch", hide_index=True)

            return

    # Fall back to annotations-based analysis
    if st.session_state.get('annotations') is not None:
        annotations = st.session_state.annotations

        if 'sequence_length' in annotations.columns:
            lengths = annotations['sequence_length'].dropna().astype(int).tolist()

            # Size distribution analysis
            size_analysis = analyze_size_distribution(lengths)

            # Metrics
            col1, col2, col3, col4 = st.columns(4)

            with col1:
                st.metric("Total Sequences", f"{size_analysis['total_reads']:,}")

            with col2:
                st.metric("Mean Length", f"{size_analysis['mean_length']:.1f} nt")

            with col3:
                st.metric("Median Length", f"{size_analysis['median_length']:.0f} nt")

            with col4:
                st.metric("Range", f"{size_analysis['min_length']}-{size_analysis['max_length']} nt")

            # Size distribution plot
            fig = plot_read_length_distribution(
                lengths,
                title="Sequence Length Distribution",
                expected_ranges=config.qc.size_ranges
            )
            st.plotly_chart(fig, width="stretch")

            # Categorized by RNA type
            st.subheader("Sequences by Expected RNA Type Size")

            cat_df = pd.DataFrame(size_analysis['categorized']).T
            cat_df = cat_df.sort_values('count', ascending=False)

            st.dataframe(cat_df, width="stretch")

            # Peak detection
            if size_analysis['peaks']:
                st.subheader("Detected Peaks")
                peaks_df = pd.DataFrame(size_analysis['peaks'])
                st.dataframe(peaks_df, width="stretch")

        else:
            st.warning("No 'sequence_length' column found in annotations")

    else:
        st.info("Run FASTQ QC analysis first, or upload annotations with sequence length information.")


def analyze_fastq_contamination(fastq_file: Path, max_reads: int = 50000) -> Dict:
    """
    Analyze FASTQ file for potential contamination indicators (cached wrapper)

    Args:
        fastq_file: Path to FASTQ file
        max_reads: Maximum reads to sample

    Returns:
        Dictionary with contamination metrics
    """
    file_hash = hash_file(fastq_file)
    return _analyze_fastq_contamination_cached(file_hash, str(fastq_file), max_reads)


@st.cache_data(show_spinner=False, ttl=3600)
def _analyze_fastq_contamination_cached(_file_hash: str, fastq_file_path: str, max_reads: int = 50000) -> Dict:
    """Cached implementation of contamination analysis"""
    import gzip
    from collections import Counter

    fastq_file = Path(fastq_file_path)

    # Common adapter sequences (first 12bp) for multiple platforms
    adapters = {
        # Illumina adapters
        'Illumina Universal': 'AGATCGGAAGAG',
        'Illumina Small RNA': 'TGGAATTCTCGG',
        'Nextera': 'CTGTCTCTTATA',
        # Ion Torrent adapters
        'Ion Torrent A': 'CCATCTCATCCC',
        'Ion Torrent P1': 'CCTCTCTATGGG',
        'Ion Torrent smRNA': 'ATCACCGACTGC',
        # BGI/MGI adapters
        'BGI/MGI': 'AGTCGGAGGCCA',
        # Nanopore adapters
        'Nanopore cDNA': 'ACTTGCCTGTCG',
        # QIAseq adapters
        'QIAseq': 'AACTGTAGGCAC',
    }

    adapter_counts = {name: 0 for name in adapters}
    poly_a_count = 0
    poly_t_count = 0
    n_content_total = 0
    total_bases = 0
    read_count = 0

    # K-mer analysis (6-mers for overrepresentation)
    kmer_size = 6
    kmer_counts = Counter()

    # Open file
    opener = gzip.open if str(fastq_file).endswith('.gz') else open
    mode = 'rt' if str(fastq_file).endswith('.gz') else 'r'

    try:
        with opener(fastq_file, mode) as f:
            while read_count < max_reads:
                header = f.readline()
                if not header:
                    break

                sequence = f.readline().strip().upper()
                plus = f.readline()
                quality = f.readline()

                if not sequence:
                    break

                read_count += 1
                seq_len = len(sequence)
                total_bases += seq_len

                # Check for adapters
                for name, adapter in adapters.items():
                    if adapter in sequence:
                        adapter_counts[name] += 1

                # Check for poly-A/T (>= 10 consecutive)
                if 'AAAAAAAAAA' in sequence:
                    poly_a_count += 1
                if 'TTTTTTTTTT' in sequence:
                    poly_t_count += 1

                # Count N bases
                n_content_total += sequence.count('N')

                # K-mer counting (sample every 10th read for speed)
                if read_count % 10 == 0:
                    for i in range(len(sequence) - kmer_size + 1):
                        kmer = sequence[i:i+kmer_size]
                        if 'N' not in kmer:
                            kmer_counts[kmer] += 1

        # Calculate percentages
        adapter_pct = {name: 100 * count / read_count if read_count > 0 else 0
                      for name, count in adapter_counts.items()}

        # Find overrepresented k-mers (>1% of sampled k-mers)
        total_kmers = sum(kmer_counts.values())
        overrep_kmers = []
        if total_kmers > 0:
            for kmer, count in kmer_counts.most_common(20):
                pct = 100 * count / total_kmers
                if pct > 0.5:  # More than 0.5% is notable for 6-mers
                    overrep_kmers.append({
                        'kmer': kmer,
                        'count': count,
                        'percentage': pct
                    })

        return {
            'status': 'success',
            'reads_analyzed': read_count,
            'adapter_content': adapter_pct,
            'adapters_detected': [name for name, pct in adapter_pct.items() if pct > 1],
            'poly_a_reads': poly_a_count,
            'poly_a_pct': 100 * poly_a_count / read_count if read_count > 0 else 0,
            'poly_t_reads': poly_t_count,
            'poly_t_pct': 100 * poly_t_count / read_count if read_count > 0 else 0,
            'n_content_pct': 100 * n_content_total / total_bases if total_bases > 0 else 0,
            'overrepresented_kmers': overrep_kmers[:10]
        }

    except Exception as e:
        return {
            'status': 'error',
            'error': str(e)
        }


def render_contamination_check():
    """Render contamination check section"""
    st.subheader("⚠️ Contamination Check")

    # Check for FASTQ QC results first
    qc_results = st.session_state.get('qc_results')
    if qc_results and 'file_paths' in qc_results:
        file_paths = qc_results.get('file_paths', [])

        if file_paths:
            st.info("🔬 **FASTQ-level contamination analysis** - checks for adapters, poly-A/T, and overrepresented sequences")

            # Check if contamination analysis already done
            if 'contamination_results' not in st.session_state:
                if st.button("🔍 Run Contamination Analysis", type="primary"):
                    progress = st.progress(0)
                    status = st.empty()

                    contamination_results = {}
                    valid_files = [p for p in file_paths if isinstance(p, Path) and p.exists()]

                    for i, file_path in enumerate(valid_files):
                        status.text(f"Analyzing {file_path.name}...")
                        progress.progress((i + 1) / len(valid_files))

                        result = analyze_fastq_contamination(file_path)
                        result['filename'] = file_path.name
                        contamination_results[file_path.name] = result

                    progress.progress(1.0)
                    status.text("Analysis complete!")

                    st.session_state.contamination_results = contamination_results
                    st.rerun()
            else:
                # Display contamination results
                contamination_results = st.session_state.contamination_results

                st.success(f"✅ Contamination analysis for {len(contamination_results)} files")

                # Adapter content summary
                st.subheader("🧬 Adapter Content")

                adapter_data = []
                for name, result in contamination_results.items():
                    if result.get('status') == 'success':
                        row = {'Sample': name}
                        row.update(result['adapter_content'])
                        adapter_data.append(row)

                if adapter_data:
                    adapter_df = pd.DataFrame(adapter_data)

                    # Highlight high adapter content
                    def highlight_high(val):
                        if isinstance(val, (int, float)) and val > 5:
                            return 'background-color: #ffcccc'
                        return ''

                    # Use map (newer pandas) or applymap (older pandas)
                    try:
                        styled_df = adapter_df.style.map(highlight_high, subset=adapter_df.columns[1:])
                    except AttributeError:
                        styled_df = adapter_df.style.applymap(highlight_high, subset=adapter_df.columns[1:])
                    st.dataframe(styled_df, width="stretch")

                    # Warning for high adapter content
                    high_adapter_samples = []
                    for name, result in contamination_results.items():
                        if result.get('status') == 'success':
                            if result.get('adapters_detected'):
                                high_adapter_samples.append((name, result['adapters_detected']))

                    if high_adapter_samples:
                        st.warning("⚠️ **High adapter content detected:**")
                        for sample, adapters in high_adapter_samples:
                            st.markdown(f"- **{sample}**: {', '.join(adapters)}")
                        st.markdown("*Consider trimming adapters before alignment.*")

                # Poly-A/T and N content
                st.subheader("🔤 Sequence Quality Indicators")

                quality_data = []
                for name, result in contamination_results.items():
                    if result.get('status') == 'success':
                        quality_data.append({
                            'Sample': name,
                            'Poly-A %': f"{result['poly_a_pct']:.2f}%",
                            'Poly-T %': f"{result['poly_t_pct']:.2f}%",
                            'N Content %': f"{result['n_content_pct']:.3f}%",
                            'Reads Analyzed': f"{result['reads_analyzed']:,}"
                        })

                if quality_data:
                    st.dataframe(pd.DataFrame(quality_data), width="stretch", hide_index=True)

                # Overrepresented k-mers
                st.subheader("🔁 Overrepresented Sequences (6-mers)")

                for name, result in contamination_results.items():
                    if result.get('status') == 'success' and result.get('overrepresented_kmers'):
                        with st.expander(f"📄 {name}"):
                            kmer_df = pd.DataFrame(result['overrepresented_kmers'])
                            kmer_df['percentage'] = kmer_df['percentage'].apply(lambda x: f"{x:.3f}%")
                            st.dataframe(kmer_df, width="stretch", hide_index=True)

                # Re-run option
                if st.button("🔄 Re-run Analysis"):
                    del st.session_state.contamination_results
                    st.rerun()

            return

    # Fall back to count matrix analysis
    if st.session_state.get('count_matrix') is not None and st.session_state.get('annotations') is not None:
        counts = st.session_state.count_matrix
        annotations = st.session_state.annotations

        # Check rRNA contamination
        rrna_analysis = check_rrna_contamination(annotations, counts)

        if 'error' not in rrna_analysis:
            # Summary metrics
            col1, col2 = st.columns(2)

            with col1:
                st.metric(
                    "Mean rRNA %",
                    f"{rrna_analysis['mean_rrna_pct']:.1f}%",
                    delta=None if rrna_analysis['mean_rrna_pct'] < 20 else "High!",
                    delta_color="off" if rrna_analysis['mean_rrna_pct'] < 20 else "inverse"
                )

            with col2:
                n_high = len(rrna_analysis['samples_high_rrna'])
                st.metric(
                    "Samples with >50% rRNA",
                    n_high,
                    delta=None if n_high == 0 else "Warning",
                    delta_color="off" if n_high == 0 else "inverse"
                )

            # rRNA percentage by sample
            import plotly.express as px

            rrna_df = pd.DataFrame({
                'Sample': list(rrna_analysis['rrna_percentage'].keys()),
                'rRNA %': list(rrna_analysis['rrna_percentage'].values())
            })

            fig = px.bar(
                rrna_df,
                x='Sample',
                y='rRNA %',
                title='rRNA Contamination by Sample',
                color='rRNA %',
                color_continuous_scale='RdYlGn_r'
            )
            fig.add_hline(y=20, line_dash="dash", line_color="orange",
                         annotation_text="20% threshold")
            fig.update_layout(xaxis_tickangle=-45)
            st.plotly_chart(fig, width="stretch")

            # Warning for high rRNA samples
            if rrna_analysis['samples_high_rrna']:
                st.warning(f"⚠️ High rRNA contamination detected in: {', '.join(rrna_analysis['samples_high_rrna'])}")
                st.markdown("""
                **Recommendations:**
                - Consider filtering rRNA reads before downstream analysis
                - Check library preparation protocol
                - Verify rRNA depletion steps
                """)

            # Option to filter rRNA
            st.divider()
            if st.checkbox("Remove rRNA reads from count matrix"):
                # Filter out rRNA
                rrna_mask = annotations['RNA_category'] != 'rRNA'
                filtered_counts = counts[counts.index.isin(annotations[rrna_mask].index)]

                st.success(f"Filtered: {len(counts)} → {len(filtered_counts)} features")

                if st.button("💾 Save Filtered Counts"):
                    st.session_state.count_matrix = filtered_counts
                    st.success("Filtered count matrix saved!")

        else:
            st.error(rrna_analysis['error'])

        # RNA type distribution
        st.divider()
        st.subheader("RNA Type Distribution")

        if 'RNA_category' in annotations.columns:
            fig = plot_rna_type_distribution(annotations, counts)
            st.plotly_chart(fig, width="stretch")

    else:
        st.info("""
        **To run contamination checks:**
        1. **For FASTQ data:** Go to 'Upload Data' tab, select 'FASTQ files', and run QC analysis first
        2. **For count data:** Upload both count matrix and annotations with RNA_category information
        """)


def analyze_bam_file(bam_file: Path, max_reads: int = 100000) -> Dict:
    """
    Analyze BAM file for post-alignment QC metrics (cached wrapper)

    Args:
        bam_file: Path to BAM file
        max_reads: Maximum reads to sample

    Returns:
        Dictionary with alignment QC metrics
    """
    file_hash = hash_file(bam_file)
    return _analyze_bam_file_cached(file_hash, str(bam_file), max_reads)


@st.cache_data(show_spinner=False, ttl=3600)
def _analyze_bam_file_cached(_file_hash: str, bam_file_path: str, max_reads: int = 100000) -> Dict:
    """Cached implementation of BAM analysis"""
    import subprocess

    bam_file = Path(bam_file_path)
    result = {
        'status': 'success',
        'total_reads': 0,
        'mapped_reads': 0,
        'unmapped_reads': 0,
        'mapping_rate': 0.0,
        'multi_mapped': 0,
        'uniquely_mapped': 0,
        'unique_rate': 0.0,
        'length_distribution': {},
        'mapq_distribution': {},
        'strand_forward': 0,
        'strand_reverse': 0,
        'strand_bias': 0.0,
        'mean_mapq': 0.0,
        'lengths': [],
        'mapqs': []
    }

    try:
        # Get flagstat for overall stats
        flagstat = subprocess.run(
            ['samtools', 'flagstat', str(bam_file)],
            capture_output=True, text=True
        )

        if flagstat.returncode == 0:
            for line in flagstat.stdout.split('\n'):
                if 'in total' in line:
                    result['total_reads'] = int(line.split()[0])
                elif 'mapped (' in line and 'primary' not in line:
                    result['mapped_reads'] = int(line.split()[0])

            if result['total_reads'] > 0:
                result['mapping_rate'] = 100 * result['mapped_reads'] / result['total_reads']
                result['unmapped_reads'] = result['total_reads'] - result['mapped_reads']

        # Sample reads for detailed analysis
        view_cmd = ['samtools', 'view', str(bam_file)]
        view_result = subprocess.run(view_cmd, capture_output=True, text=True)

        if view_result.returncode == 0:
            lengths = []
            mapqs = []
            nh_tags = []  # Number of hits (multi-mapping)
            forward_count = 0
            reverse_count = 0
            read_count = 0

            for line in view_result.stdout.split('\n')[:max_reads]:
                if not line or line.startswith('@'):
                    continue

                fields = line.split('\t')
                if len(fields) < 11:
                    continue

                read_count += 1
                flag = int(fields[1])
                mapq = int(fields[4])
                seq = fields[9]

                # Read length
                lengths.append(len(seq))

                # Mapping quality
                mapqs.append(mapq)

                # Strand
                if flag & 16:  # Reverse strand
                    reverse_count += 1
                else:
                    forward_count += 1

                # Check for NH tag (number of hits)
                for field in fields[11:]:
                    if field.startswith('NH:i:'):
                        nh = int(field.split(':')[2])
                        nh_tags.append(nh)
                        break

            # Calculate statistics
            if lengths:
                result['lengths'] = lengths[:10000]
                result['length_distribution'] = dict(Counter(lengths))
                result['mean_aligned_length'] = float(np.mean(lengths))
                result['median_aligned_length'] = float(np.median(lengths))

            if mapqs:
                result['mapqs'] = mapqs[:10000]
                result['mapq_distribution'] = dict(Counter(mapqs))
                result['mean_mapq'] = float(np.mean(mapqs))

            if nh_tags:
                result['uniquely_mapped'] = sum(1 for nh in nh_tags if nh == 1)
                result['multi_mapped'] = sum(1 for nh in nh_tags if nh > 1)
                if len(nh_tags) > 0:
                    result['unique_rate'] = 100 * result['uniquely_mapped'] / len(nh_tags)

            result['strand_forward'] = forward_count
            result['strand_reverse'] = reverse_count
            if forward_count + reverse_count > 0:
                result['strand_bias'] = abs(0.5 - forward_count / (forward_count + reverse_count))

    except Exception as e:
        result['status'] = 'error'
        result['error'] = str(e)

    return result


def render_post_alignment_qc():
    """Render post-alignment QC section"""
    st.subheader("🎯 Post-Alignment Quality Control")

    # Check for alignment results
    alignment_results = st.session_state.get('alignment_results')
    bam_files = st.session_state.get('alignment_bam_files', [])

    if not alignment_results or not bam_files:
        st.warning("""
        **No alignment data found.**

        Please run the alignment module first:
        1. Go to **🔗 Alignment** in the sidebar
        2. Set up your reference database
        3. Run alignment on your samples
        4. Return here to analyze the results
        """)

        # Option to upload BAM files directly
        st.divider()
        st.markdown("**Or upload BAM files directly:**")
        uploaded_bams = st.file_uploader(
            "Upload BAM files",
            type=['bam'],
            accept_multiple_files=True
        )

        if uploaded_bams:
            st.info(f"Uploaded {len(uploaded_bams)} BAM files")
            # Save to temp and analyze
            temp_dir = Path(tempfile.mkdtemp())
            bam_files = []
            for uf in uploaded_bams:
                bam_path = temp_dir / uf.name
                with open(bam_path, 'wb') as f:
                    f.write(uf.getbuffer())
                bam_files.append(str(bam_path))

            st.session_state.uploaded_bam_files = bam_files

        if not bam_files:
            return

    # Use uploaded BAMs if available
    if not bam_files and st.session_state.get('uploaded_bam_files'):
        bam_files = st.session_state.uploaded_bam_files

    st.success(f"✅ Found {len(bam_files)} BAM files for analysis")

    # Run post-alignment QC
    if st.button("🔍 Run Post-Alignment QC", type="primary"):
        progress = st.progress(0)
        status = st.empty()

        post_qc_results = {}

        for i, bam_file in enumerate(bam_files):
            bam_path = Path(bam_file)
            sample_name = bam_path.stem.replace('_sorted', '').replace('_aligned', '')

            status.text(f"Analyzing {sample_name}...")
            progress.progress((i + 1) / len(bam_files))

            result = analyze_bam_file(bam_path)
            result['sample'] = sample_name
            post_qc_results[sample_name] = result

        st.session_state.post_alignment_qc = post_qc_results
        st.success("✅ Post-alignment QC complete!")
        st.rerun()

    # Display results
    post_qc = st.session_state.get('post_alignment_qc')

    if not post_qc:
        return

    st.divider()

    # Summary table
    st.markdown("### 📊 Alignment Summary")

    summary_data = []
    for sample, result in post_qc.items():
        if result['status'] == 'success':
            summary_data.append({
                'Sample': sample,
                'Total Reads': result.get('total_reads', 0),
                'Mapped': result.get('mapped_reads', 0),
                'Mapping Rate': f"{result.get('mapping_rate', 0):.1f}%",
                'Uniquely Mapped': result.get('uniquely_mapped', 'N/A'),
                'Multi-Mapped': result.get('multi_mapped', 'N/A'),
                'Mean MAPQ': f"{result.get('mean_mapq', 0):.1f}",
                'Strand Bias': f"{result.get('strand_bias', 0):.3f}"
            })

    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        st.dataframe(summary_df, width="stretch")

    # Visualizations
    import plotly.express as px
    import plotly.graph_objects as go

    col1, col2 = st.columns(2)

    with col1:
        # Mapping rate bar chart
        st.markdown("### 📈 Mapping Rates")
        rates_data = [{'Sample': s, 'Mapping Rate': r.get('mapping_rate', 0)}
                      for s, r in post_qc.items() if r['status'] == 'success']

        if rates_data:
            fig = px.bar(
                pd.DataFrame(rates_data),
                x='Sample', y='Mapping Rate',
                color='Mapping Rate',
                color_continuous_scale='RdYlGn',
                title='Mapping Rate by Sample'
            )
            fig.add_hline(y=70, line_dash="dash", line_color="orange",
                         annotation_text="70% threshold")
            fig.update_layout(xaxis_tickangle=-45)
            st.plotly_chart(fig, width="stretch")

    with col2:
        # Unique vs Multi-mapped
        st.markdown("### 🎯 Unique vs Multi-Mapped")
        mapping_data = []
        for sample, result in post_qc.items():
            if result['status'] == 'success':
                unique = result.get('uniquely_mapped', 0)
                multi = result.get('multi_mapped', 0)
                if unique or multi:
                    mapping_data.append({'Sample': sample, 'Type': 'Unique', 'Reads': unique})
                    mapping_data.append({'Sample': sample, 'Type': 'Multi-mapped', 'Reads': multi})

        if mapping_data:
            fig = px.bar(
                pd.DataFrame(mapping_data),
                x='Sample', y='Reads', color='Type',
                barmode='stack',
                title='Unique vs Multi-Mapped Reads'
            )
            fig.update_layout(xaxis_tickangle=-45)
            st.plotly_chart(fig, width="stretch")

    # Aligned read length distribution
    st.markdown("### 📏 Aligned Read Length Distribution")

    # Combine all lengths
    all_lengths = []
    for sample, result in post_qc.items():
        if result['status'] == 'success' and result.get('lengths'):
            for length in result['lengths']:
                all_lengths.append({'Sample': sample, 'Length': length})

    if all_lengths:
        lengths_df = pd.DataFrame(all_lengths)

        fig = px.histogram(
            lengths_df,
            x='Length',
            color='Sample',
            barmode='overlay',
            nbins=50,
            title='Aligned Read Length Distribution',
            opacity=0.7
        )
        fig.update_layout(xaxis_title="Read Length (nt)", yaxis_title="Count")
        st.plotly_chart(fig, width="stretch")

    # Mapping quality distribution
    st.markdown("### 🎚️ Mapping Quality Distribution")

    all_mapqs = []
    for sample, result in post_qc.items():
        if result['status'] == 'success' and result.get('mapqs'):
            for mapq in result['mapqs']:
                all_mapqs.append({'Sample': sample, 'MAPQ': mapq})

    if all_mapqs:
        mapq_df = pd.DataFrame(all_mapqs)

        fig = px.histogram(
            mapq_df,
            x='MAPQ',
            color='Sample',
            barmode='overlay',
            title='Mapping Quality (MAPQ) Distribution',
            opacity=0.7
        )
        fig.update_layout(xaxis_title="Mapping Quality", yaxis_title="Count")
        st.plotly_chart(fig, width="stretch")

    # Strand bias
    st.markdown("### ↔️ Strand Distribution")

    strand_data = []
    for sample, result in post_qc.items():
        if result['status'] == 'success':
            strand_data.append({
                'Sample': sample,
                'Strand': 'Forward (+)',
                'Reads': result.get('strand_forward', 0)
            })
            strand_data.append({
                'Sample': sample,
                'Strand': 'Reverse (-)',
                'Reads': result.get('strand_reverse', 0)
            })

    if strand_data:
        fig = px.bar(
            pd.DataFrame(strand_data),
            x='Sample', y='Reads', color='Strand',
            barmode='group',
            title='Strand Distribution'
        )
        fig.update_layout(xaxis_tickangle=-45)
        st.plotly_chart(fig, width="stretch")

    # Quality warnings
    st.divider()
    st.markdown("### ⚠️ Quality Flags")

    warnings = []
    for sample, result in post_qc.items():
        if result['status'] == 'success':
            # Low mapping rate
            if result.get('mapping_rate', 0) < 50:
                warnings.append(f"⚠️ **{sample}**: Low mapping rate ({result['mapping_rate']:.1f}%)")
            # High strand bias
            if result.get('strand_bias', 0) > 0.2:
                warnings.append(f"⚠️ **{sample}**: High strand bias ({result['strand_bias']:.3f})")
            # Low MAPQ
            if result.get('mean_mapq', 0) < 20:
                warnings.append(f"⚠️ **{sample}**: Low mean mapping quality ({result['mean_mapq']:.1f})")

    if warnings:
        for w in warnings:
            st.markdown(w)
    else:
        st.success("✅ All samples passed quality checks!")

    # Download results
    st.divider()
    csv = summary_df.to_csv(index=False) if summary_data else ""
    st.download_button(
        "📥 Download Post-Alignment QC Summary",
        csv,
        "post_alignment_qc.csv",
        "text/csv"
    )


def render_pre_post_comparison():
    """Render pre vs post alignment comparison"""
    st.subheader("⚖️ Pre-Alignment vs Post-Alignment Comparison")

    # Check for both pre and post alignment data
    pre_qc = st.session_state.get('qc_results')
    post_qc = st.session_state.get('post_alignment_qc')
    alignment_results = st.session_state.get('alignment_results')

    has_pre = pre_qc is not None and len(pre_qc) > 0
    has_post = post_qc is not None and len(post_qc) > 0
    has_alignment = alignment_results is not None and len(alignment_results) > 0

    if not has_pre:
        st.warning("""
        **Pre-alignment QC data not found.**

        Please run FASTQ QC first:
        1. Go to **Pre-Alignment Stats** tab
        2. Upload or select your FASTQ files
        3. Run QC analysis
        """)

    if not has_post and not has_alignment:
        st.warning("""
        **Post-alignment data not found.**

        Please run alignment first:
        1. Go to **🔗 Alignment** module
        2. Run alignment on your samples
        3. Go to **Post-Alignment QC** tab and run analysis
        """)

    if not (has_pre and (has_post or has_alignment)):
        return

    st.success("✅ Both pre-alignment and post-alignment data available!")

    # Build comparison data
    comparison_data = []

    for sample_name, pre_result in pre_qc.items():
        if pre_result.get('status') != 'success':
            continue

        # Find matching post-alignment data
        # Try to match by sample name (may have different suffixes)
        post_result = None
        align_result = None

        # Check post_qc
        if has_post:
            for post_name, post_r in post_qc.items():
                if sample_name in post_name or post_name in sample_name:
                    post_result = post_r
                    break

        # Check alignment_results
        if has_alignment:
            for align_name, align_r in alignment_results.items():
                if sample_name in align_name or align_name in sample_name:
                    align_result = align_r
                    break

        row = {
            'Sample': sample_name,
            # Pre-alignment metrics
            'Input Reads': pre_result.get('total_reads', 0),
            'Pre Mean Length': f"{pre_result.get('mean_length', 0):.1f}",
            'Pre Mean Quality': f"{pre_result.get('mean_quality', 0):.1f}",
            'Pre GC%': f"{pre_result.get('mean_gc', 0):.1f}",
        }

        # Add post-alignment metrics
        if align_result and align_result.get('stats'):
            stats = align_result['stats']
            row['Aligned Reads'] = stats.get('aligned_reads', 0)
            row['Alignment Rate'] = f"{stats.get('alignment_rate', 0):.1f}%"
            row['Not Aligned'] = stats.get('not_aligned', 0)
            row['Suppressed'] = stats.get('suppressed', 0)
        elif post_result:
            row['Aligned Reads'] = post_result.get('mapped_reads', 0)
            row['Alignment Rate'] = f"{post_result.get('mapping_rate', 0):.1f}%"
            row['Not Aligned'] = post_result.get('unmapped_reads', 0)
            row['Suppressed'] = 'N/A'

        if post_result:
            row['Post Mean Length'] = f"{post_result.get('mean_aligned_length', 0):.1f}"
            row['Mean MAPQ'] = f"{post_result.get('mean_mapq', 0):.1f}"
            row['Unique Rate'] = f"{post_result.get('unique_rate', 0):.1f}%"
        else:
            row['Post Mean Length'] = 'N/A'
            row['Mean MAPQ'] = 'N/A'
            row['Unique Rate'] = 'N/A'

        comparison_data.append(row)

    if not comparison_data:
        st.warning("Could not match samples between pre and post alignment data.")
        return

    # Display comparison table
    st.markdown("### 📊 Sample Comparison Table")
    comparison_df = pd.DataFrame(comparison_data)
    st.dataframe(comparison_df, width="stretch")

    # Visualizations
    import plotly.express as px
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    st.divider()

    # Read count flow (Sankey-style bar chart)
    st.markdown("### 📈 Read Count Flow: Input → Aligned")

    flow_data = []
    for row in comparison_data:
        sample = row['Sample']
        input_reads = row.get('Input Reads', 0)
        aligned = row.get('Aligned Reads', 0) if row.get('Aligned Reads') != 'N/A' else 0
        not_aligned = row.get('Not Aligned', 0) if row.get('Not Aligned') != 'N/A' else 0
        suppressed = row.get('Suppressed', 0) if row.get('Suppressed') != 'N/A' else 0

        flow_data.append({'Sample': sample, 'Category': 'Input', 'Reads': input_reads})
        flow_data.append({'Sample': sample, 'Category': 'Aligned', 'Reads': aligned})
        if not_aligned:
            flow_data.append({'Sample': sample, 'Category': 'Not Aligned', 'Reads': not_aligned})
        if suppressed:
            flow_data.append({'Sample': sample, 'Category': 'Suppressed', 'Reads': suppressed})

    if flow_data:
        fig = px.bar(
            pd.DataFrame(flow_data),
            x='Sample', y='Reads', color='Category',
            barmode='group',
            title='Read Count at Each Stage',
            color_discrete_map={
                'Input': '#3498db',
                'Aligned': '#27ae60',
                'Not Aligned': '#e74c3c',
                'Suppressed': '#f39c12'
            }
        )
        fig.update_layout(xaxis_tickangle=-45)
        st.plotly_chart(fig, width="stretch")

    # Length comparison
    st.markdown("### 📏 Read Length: Pre vs Post Alignment")

    col1, col2 = st.columns(2)

    with col1:
        # Pre-alignment length distribution
        st.markdown("**Pre-Alignment (Raw FASTQ)**")
        pre_lengths = []
        for sample, result in pre_qc.items():
            if result.get('status') == 'success' and result.get('lengths'):
                for length in result['lengths'][:5000]:
                    pre_lengths.append({'Sample': sample, 'Length': length})

        if pre_lengths:
            fig = px.histogram(
                pd.DataFrame(pre_lengths),
                x='Length',
                color='Sample',
                barmode='overlay',
                opacity=0.7,
                title='Raw Read Lengths'
            )
            st.plotly_chart(fig, width="stretch")

    with col2:
        # Post-alignment length distribution
        st.markdown("**Post-Alignment (Mapped Reads)**")
        if has_post:
            post_lengths = []
            for sample, result in post_qc.items():
                if result.get('status') == 'success' and result.get('lengths'):
                    for length in result['lengths'][:5000]:
                        post_lengths.append({'Sample': sample, 'Length': length})

            if post_lengths:
                fig = px.histogram(
                    pd.DataFrame(post_lengths),
                    x='Length',
                    color='Sample',
                    barmode='overlay',
                    opacity=0.7,
                    title='Aligned Read Lengths'
                )
                st.plotly_chart(fig, width="stretch")
        else:
            st.info("Run Post-Alignment QC to see aligned read lengths")

    # Summary statistics comparison
    st.divider()
    st.markdown("### 📋 Key Metrics Summary")

    col1, col2, col3, col4 = st.columns(4)

    # Calculate aggregated metrics
    total_input = sum(r.get('Input Reads', 0) for r in comparison_data)
    total_aligned = sum(r.get('Aligned Reads', 0) for r in comparison_data
                        if r.get('Aligned Reads') != 'N/A' and r.get('Aligned Reads'))

    avg_pre_length = np.mean([float(r.get('Pre Mean Length', 0)) for r in comparison_data])
    avg_post_length = np.mean([float(r.get('Post Mean Length', 0)) for r in comparison_data
                               if r.get('Post Mean Length') != 'N/A'])

    with col1:
        st.metric(
            "Total Input Reads",
            f"{total_input:,}",
            help="Total reads before alignment"
        )

    with col2:
        st.metric(
            "Total Aligned Reads",
            f"{total_aligned:,}",
            delta=f"{100*total_aligned/total_input:.1f}%" if total_input > 0 else "0%",
            help="Total reads successfully aligned"
        )

    with col3:
        st.metric(
            "Avg Pre-Alignment Length",
            f"{avg_pre_length:.1f} nt",
            help="Mean read length before alignment"
        )

    with col4:
        if avg_post_length > 0:
            delta = avg_post_length - avg_pre_length
            st.metric(
                "Avg Post-Alignment Length",
                f"{avg_post_length:.1f} nt",
                delta=f"{delta:+.1f} nt",
                help="Mean aligned read length"
            )
        else:
            st.metric("Avg Post-Alignment Length", "N/A")

    # Insights and recommendations
    st.divider()
    st.markdown("### 💡 Insights & Recommendations")

    insights = []

    # Calculate overall alignment rate
    if total_input > 0 and total_aligned > 0:
        overall_rate = 100 * total_aligned / total_input

        if overall_rate >= 70:
            insights.append("✅ **Good alignment rate** - Most reads are mapping to the reference")
        elif overall_rate >= 40:
            insights.append("⚠️ **Moderate alignment rate** - Consider checking adapter trimming or reference database")
        else:
            insights.append("❌ **Low alignment rate** - Check: wrong reference, untrimmed adapters, or organism mismatch")

    # Length changes
    if avg_post_length > 0:
        length_change = avg_post_length - avg_pre_length
        if abs(length_change) < 2:
            insights.append("✅ **Consistent read lengths** - Aligned reads maintain similar lengths")
        elif length_change < -3:
            insights.append("⚠️ **Shorter aligned reads** - Some reads may be partially mapping or soft-clipped")
        elif length_change > 3:
            insights.append("ℹ️ **Longer aligned reads** - Short reads may have been filtered out")

    # Quality recommendations
    low_rate_samples = [r['Sample'] for r in comparison_data
                        if r.get('Alignment Rate', '0%').replace('%', '').replace('N/A', '0')
                        and float(r.get('Alignment Rate', '0%').replace('%', '').replace('N/A', '0')) < 50]

    if low_rate_samples:
        insights.append(f"⚠️ **Low alignment in samples:** {', '.join(low_rate_samples[:3])}")
        insights.append("   → Try different reference database or check library quality")

    if insights:
        for insight in insights:
            st.markdown(insight)
    else:
        st.info("Run both pre-alignment and post-alignment QC to see detailed insights")

    # Download comparison
    st.divider()
    csv = comparison_df.to_csv(index=False)
    st.download_button(
        "📥 Download Pre/Post Comparison",
        csv,
        "pre_post_comparison.csv",
        "text/csv"
    )


def render_multi_sample_qc_overlay():
    """Render multi-sample QC overlays with outlier detection"""
    st.subheader("🔍 Multi-Sample QC Overlay & Outlier Detection")

    st.markdown("""
    Compare QC metrics across all samples simultaneously to identify outliers
    and batch effects. Outliers are detected using the Median Absolute Deviation (MAD)
    method, which is robust to non-normal distributions.
    """)

    # Check for QC data
    qc_results = st.session_state.get('qc_results')
    if not qc_results or 'analysis' not in qc_results:
        st.warning("⚠️ Run QC analysis first to enable multi-sample comparison.")
        return

    analysis = qc_results['analysis']
    successful = {k: v for k, v in analysis.items() if v.get('status') == 'success'}

    if len(successful) < 3:
        st.info("ℹ️ Need at least 3 samples for meaningful outlier detection.")

    # Build metrics DataFrame
    metrics_data = []
    for sample_name, result in successful.items():
        metrics_data.append({
            'Sample': sample_name,
            'Total Reads': result.get('total_reads', 0),
            'Mean Length': result.get('mean_length', 0),
            'Mean Quality': result.get('mean_quality', 0),
            'GC %': result.get('mean_gc', 0),
            'Length Std': result.get('std_length', 0)
        })

    metrics_df = pd.DataFrame(metrics_data)

    if len(metrics_df) == 0:
        st.warning("No QC data available")
        return

    # Outlier detection settings
    st.markdown("#### Outlier Detection Settings")

    col1, col2 = st.columns(2)

    with col1:
        mad_threshold = st.slider(
            "MAD Threshold",
            1.0, 5.0, 3.0, 0.5,
            help="Samples > threshold * MAD from median are flagged as outliers"
        )

    with col2:
        metrics_to_check = st.multiselect(
            "Metrics to check for outliers",
            options=['Total Reads', 'Mean Length', 'Mean Quality', 'GC %'],
            default=['Total Reads', 'Mean Quality']
        )

    # Detect outliers
    outlier_results = detect_qc_outliers(metrics_df, metrics_to_check, mad_threshold)

    # Display outlier summary
    st.divider()
    st.markdown("#### 🚨 Outlier Detection Results")

    n_outliers = len(outlier_results['outlier_samples'])

    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Total Samples", len(metrics_df))
    with col2:
        st.metric("Outlier Samples", n_outliers,
                 delta="Flag" if n_outliers > 0 else "OK",
                 delta_color="inverse" if n_outliers > 0 else "normal")
    with col3:
        st.metric("Pass Rate", f"{100 * (len(metrics_df) - n_outliers) / len(metrics_df):.0f}%")

    # Show outliers
    if outlier_results['outlier_samples']:
        st.warning(f"⚠️ **Outlier samples detected:** {', '.join(outlier_results['outlier_samples'])}")

        # Show which metrics flagged each sample
        st.markdown("**Outlier details:**")
        for sample, flags in outlier_results['outlier_details'].items():
            if flags:
                st.markdown(f"- **{sample}**: {', '.join(flags)}")
    else:
        st.success("✅ No outliers detected with current threshold")

    # Visualizations
    st.divider()
    st.markdown("#### 📊 Multi-Sample QC Overlays")

    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    # Create overlay plots
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=('Total Reads', 'Mean Quality', 'Mean Length', 'GC Content')
    )

    outlier_samples = set(outlier_results['outlier_samples'])

    # Color samples based on outlier status
    colors = ['#e74c3c' if s in outlier_samples else '#3498db' for s in metrics_df['Sample']]

    # Total Reads
    fig.add_trace(
        go.Bar(x=metrics_df['Sample'], y=metrics_df['Total Reads'],
               marker_color=colors, showlegend=False),
        row=1, col=1
    )
    # Add median line
    median_reads = metrics_df['Total Reads'].median()
    fig.add_hline(y=median_reads, line_dash="dash", line_color="gray",
                  annotation_text="Median", row=1, col=1)

    # Mean Quality
    fig.add_trace(
        go.Bar(x=metrics_df['Sample'], y=metrics_df['Mean Quality'],
               marker_color=colors, showlegend=False),
        row=1, col=2
    )
    median_qual = metrics_df['Mean Quality'].median()
    fig.add_hline(y=median_qual, line_dash="dash", line_color="gray", row=1, col=2)

    # Mean Length
    fig.add_trace(
        go.Bar(x=metrics_df['Sample'], y=metrics_df['Mean Length'],
               marker_color=colors, showlegend=False),
        row=2, col=1
    )

    # GC %
    fig.add_trace(
        go.Bar(x=metrics_df['Sample'], y=metrics_df['GC %'],
               marker_color=colors, showlegend=False),
        row=2, col=2
    )

    fig.update_layout(
        title="Multi-Sample QC Metrics (Red = Outliers)",
        height=600,
        showlegend=False
    )
    fig.update_xaxes(tickangle=45)

    st.plotly_chart(fig, use_container_width=True)

    # Distribution plots with overlays
    st.markdown("#### 📈 Sample Distribution Overlays")

    # Read length distribution overlay
    all_lengths = qc_results.get('all_lengths', [])
    if all_lengths:
        st.markdown("##### Read Length Distribution (All Samples)")

        # Create per-sample length distributions
        fig_len = go.Figure()

        for sample_name, result in successful.items():
            if result.get('lengths'):
                sample_lengths = result['lengths'][:5000]  # Limit for performance

                # Determine color based on outlier status
                color = '#e74c3c' if sample_name in outlier_samples else '#3498db'
                opacity = 0.7 if sample_name in outlier_samples else 0.3

                fig_len.add_trace(go.Histogram(
                    x=sample_lengths,
                    name=sample_name,
                    opacity=opacity,
                    marker_color=color,
                    nbinsx=50
                ))

        fig_len.update_layout(
            title="Read Length Distribution Overlay",
            xaxis_title="Read Length (nt)",
            yaxis_title="Count",
            barmode='overlay',
            height=400
        )

        st.plotly_chart(fig_len, use_container_width=True)

    # Quality distribution overlay
    st.markdown("##### Quality Score Distribution (All Samples)")

    fig_qual = go.Figure()

    for sample_name, result in successful.items():
        if result.get('qualities'):
            sample_quals = result['qualities'][:5000]

            color = '#e74c3c' if sample_name in outlier_samples else '#2ecc71'
            opacity = 0.7 if sample_name in outlier_samples else 0.3

            fig_qual.add_trace(go.Histogram(
                x=sample_quals,
                name=sample_name,
                opacity=opacity,
                marker_color=color,
                nbinsx=40
            ))

    fig_qual.update_layout(
        title="Quality Score Distribution Overlay",
        xaxis_title="Phred Quality Score",
        yaxis_title="Count",
        barmode='overlay',
        height=400
    )

    st.plotly_chart(fig_qual, use_container_width=True)

    # PCA plot for batch effect detection
    st.markdown("#### 🎯 Sample Clustering (PCA)")

    if len(metrics_df) >= 3:
        from sklearn.preprocessing import StandardScaler
        from sklearn.decomposition import PCA

        try:
            # Prepare data for PCA
            numeric_cols = ['Total Reads', 'Mean Length', 'Mean Quality', 'GC %']
            X = metrics_df[numeric_cols].values

            # Standardize
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X)

            # PCA
            pca = PCA(n_components=2)
            X_pca = pca.fit_transform(X_scaled)

            pca_df = pd.DataFrame({
                'Sample': metrics_df['Sample'],
                'PC1': X_pca[:, 0],
                'PC2': X_pca[:, 1],
                'Status': ['Outlier' if s in outlier_samples else 'Normal'
                           for s in metrics_df['Sample']]
            })

            import plotly.express as px

            fig_pca = px.scatter(
                pca_df,
                x='PC1', y='PC2',
                color='Status',
                color_discrete_map={'Normal': '#3498db', 'Outlier': '#e74c3c'},
                text='Sample',
                title=f"PCA of QC Metrics (PC1: {pca.explained_variance_ratio_[0]*100:.1f}%, PC2: {pca.explained_variance_ratio_[1]*100:.1f}%)"
            )
            fig_pca.update_traces(textposition='top center')
            fig_pca.update_layout(height=500)

            st.plotly_chart(fig_pca, use_container_width=True)

            st.caption("Samples clustering together have similar QC profiles. Outliers may indicate technical issues or batch effects.")

        except ImportError:
            st.info("Install scikit-learn for PCA visualization: pip install scikit-learn")
        except Exception as e:
            st.warning(f"Could not generate PCA plot: {e}")

    # Download outlier report
    st.divider()
    col1, col2 = st.columns(2)

    with col1:
        # Full metrics table with outlier flags
        metrics_df['Outlier'] = metrics_df['Sample'].apply(
            lambda x: '⚠️ Yes' if x in outlier_samples else '✅ No'
        )
        csv = metrics_df.to_csv(index=False)
        st.download_button(
            "📥 Download QC Metrics with Outlier Flags",
            csv,
            "qc_metrics_outliers.csv",
            "text/csv"
        )

    with col2:
        # Outlier details
        if outlier_results['outlier_details']:
            details = "\n".join([
                f"{sample}: {', '.join(flags)}"
                for sample, flags in outlier_results['outlier_details'].items()
                if flags
            ])
            st.download_button(
                "📥 Download Outlier Details",
                details,
                "outlier_details.txt",
                "text/plain"
            )


def detect_qc_outliers(
    metrics_df: pd.DataFrame,
    metrics_to_check: List[str],
    mad_threshold: float = 3.0
) -> Dict:
    """
    Detect outlier samples using Median Absolute Deviation (MAD)

    Args:
        metrics_df: DataFrame with sample QC metrics
        metrics_to_check: List of metric columns to check
        mad_threshold: Number of MADs from median to flag as outlier

    Returns:
        Dictionary with outlier detection results
    """
    outlier_samples = set()
    outlier_details = defaultdict(list)

    for metric in metrics_to_check:
        if metric not in metrics_df.columns:
            continue

        values = metrics_df[metric].values

        # Calculate median and MAD
        median = np.median(values)
        mad = np.median(np.abs(values - median))

        # Handle case where MAD is 0 (all values same)
        if mad == 0:
            mad = np.std(values)  # Fall back to std dev
            if mad == 0:
                continue

        # Calculate modified Z-scores
        modified_z = 0.6745 * (values - median) / mad

        # Flag outliers
        for i, (sample, z_score) in enumerate(zip(metrics_df['Sample'], modified_z)):
            if abs(z_score) > mad_threshold:
                outlier_samples.add(sample)
                direction = "high" if z_score > 0 else "low"
                outlier_details[sample].append(f"{metric} ({direction})")

    return {
        'outlier_samples': list(outlier_samples),
        'outlier_details': dict(outlier_details),
        'mad_threshold': mad_threshold
    }
