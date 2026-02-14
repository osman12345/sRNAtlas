"""
Post-Alignment QC Module for sRNAtlas
Provides comprehensive quality assessment of aligned small RNA reads
"""
import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
import subprocess
import tempfile
from typing import List, Dict, Optional, Tuple
from collections import Counter, defaultdict
try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))


def check_pysam_installed() -> bool:
    """Check if pysam is available"""
    try:
        import pysam
        return True
    except ImportError:
        return False


def analyze_bam_basic(bam_file: Path, max_reads: int = 100000) -> Dict:
    """
    Perform basic QC analysis on a BAM file

    Args:
        bam_file: Path to sorted, indexed BAM file
        max_reads: Maximum reads to sample for analysis

    Returns:
        Dictionary with QC metrics
    """
    if not PYSAM_AVAILABLE:
        return {
            'status': 'error',
            'error': 'pysam is not installed. Install with: pip install pysam'
        }

    bam_file = Path(bam_file)

    # Validate BAM file exists
    if not bam_file.exists():
        return {
            'status': 'error',
            'error': f'BAM file not found: {bam_file}'
        }

    if bam_file.stat().st_size == 0:
        return {
            'status': 'error',
            'error': f'BAM file is empty: {bam_file}'
        }

    # Check for BAM index
    bai_file = Path(f"{bam_file}.bai")
    csi_file = bam_file.with_suffix('.bam.csi')
    if not bai_file.exists() and not csi_file.exists():
        return {
            'status': 'error',
            'error': f'BAM index not found. Run: samtools index {bam_file}',
            'fix': f'samtools index {bam_file}'
        }

    try:
        import pysam

        stats = {
            'total_reads': 0,
            'mapped_reads': 0,
            'unmapped_reads': 0,
            'unique_reads': 0,
            'multi_mapped': 0,
            'forward_strand': 0,
            'reverse_strand': 0,
            'lengths': [],
            'mapq_scores': [],
            'five_prime_nt': Counter(),
            'three_prime_nt': Counter(),
            'references': Counter(),
        }

        samfile = pysam.AlignmentFile(str(bam_file), "rb")

        read_count = 0
        for read in samfile.fetch(until_eof=True):
            if read_count >= max_reads:
                break

            read_count += 1
            stats['total_reads'] += 1

            if read.is_unmapped:
                stats['unmapped_reads'] += 1
                continue

            stats['mapped_reads'] += 1

            # Strand
            if read.is_reverse:
                stats['reverse_strand'] += 1
            else:
                stats['forward_strand'] += 1

            # Length
            stats['lengths'].append(read.query_length)

            # Mapping quality
            stats['mapq_scores'].append(read.mapping_quality)

            # Multi-mapping (NH tag or MAPQ=0)
            try:
                nh = read.get_tag('NH')
                if nh > 1:
                    stats['multi_mapped'] += 1
                else:
                    stats['unique_reads'] += 1
            except KeyError:
                # No NH tag, use MAPQ as proxy
                if read.mapping_quality == 0:
                    stats['multi_mapped'] += 1
                else:
                    stats['unique_reads'] += 1

            # 5' and 3' nucleotide (important for miRNA - expect U/T enrichment at 5')
            # pysam query_sequence returns the sequence as stored in BAM (already on + strand)
            # For reverse strand reads, BAM stores the reverse complement
            # So the 5' end of the ORIGINAL read is always seq[0] after pysam handles it
            seq = read.query_sequence
            if seq and len(seq) > 0:
                # pysam already handles reverse-complementing for display
                # The first base in query_sequence corresponds to the 5' end of the read
                stats['five_prime_nt'][seq[0]] += 1
                stats['three_prime_nt'][seq[-1]] += 1

            # Reference/chromosome distribution
            if read.reference_name:
                stats['references'][read.reference_name] += 1

        samfile.close()

        # Warn if analysis was truncated
        stats['_truncated'] = read_count >= max_reads
        if stats['_truncated']:
            stats['_truncation_warning'] = (
                f'Analysis based on first {max_reads:,} reads only. '
                f'Increase max_reads for full analysis.'
            )

        # Calculate rates with safe division
        if stats['total_reads'] > 0:
            stats['mapping_rate'] = 100 * stats['mapped_reads'] / stats['total_reads']
        else:
            stats['mapping_rate'] = 0.0

        if stats['mapped_reads'] > 0:
            stats['unique_rate'] = 100 * stats['unique_reads'] / stats['mapped_reads']
            stats['strand_bias'] = abs(stats['forward_strand'] - stats['reverse_strand']) / stats['mapped_reads']
            # Add strand bias QC assessment
            if stats['strand_bias'] < 0.05:
                stats['strand_bias_qc'] = 'PASS'
            elif stats['strand_bias'] < 0.15:
                stats['strand_bias_qc'] = 'WARNING - Moderate strand bias'
            else:
                stats['strand_bias_qc'] = 'FAIL - High strand bias, possible library prep issue'
        else:
            stats['unique_rate'] = 0.0
            stats['strand_bias'] = 0.0
            stats['strand_bias_qc'] = 'N/A - No mapped reads'

        # Length statistics
        if stats['lengths']:
            lengths = np.array(stats['lengths'])
            stats['mean_length'] = float(np.mean(lengths))
            stats['median_length'] = float(np.median(lengths))
            stats['min_length'] = int(np.min(lengths))
            stats['max_length'] = int(np.max(lengths))
            stats['length_distribution'] = dict(Counter(lengths))

        # MAPQ statistics
        if stats['mapq_scores']:
            mapq = np.array(stats['mapq_scores'])
            stats['mean_mapq'] = float(np.mean(mapq))
            stats['mapq_distribution'] = dict(Counter(mapq))

        stats['status'] = 'success'
        return stats

    except Exception as e:
        return {
            'status': 'error',
            'error': str(e)
        }


def get_flagstat(bam_file: Path) -> Dict:
    """
    Get samtools flagstat output for a BAM file

    Args:
        bam_file: Path to BAM file

    Returns:
        Dictionary with flagstat metrics
    """
    try:
        result = subprocess.run(
            ['samtools', 'flagstat', str(bam_file)],
            capture_output=True,
            text=True,
            timeout=300
        )

        if result.returncode != 0:
            return {'status': 'error', 'error': result.stderr}

        stats = {'raw_output': result.stdout}

        # Parse flagstat output
        for line in result.stdout.split('\n'):
            if 'in total' in line:
                stats['total'] = int(line.split()[0])
            elif 'mapped (' in line:
                stats['mapped'] = int(line.split()[0])
                # Extract percentage
                if '(' in line and '%' in line:
                    pct = line.split('(')[1].split('%')[0]
                    stats['mapped_pct'] = float(pct)
            elif 'duplicates' in line:
                stats['duplicates'] = int(line.split()[0])
            elif 'properly paired' in line:
                stats['properly_paired'] = int(line.split()[0])

        stats['status'] = 'success'
        return stats

    except FileNotFoundError:
        return {'status': 'error', 'error': 'samtools not found. Install with: conda install -c bioconda samtools'}
    except subprocess.TimeoutExpired:
        return {'status': 'error', 'error': 'samtools flagstat timed out (>5 minutes)'}
    except Exception as e:
        return {'status': 'error', 'error': str(e)}


def analyze_rna_type_distribution(bam_file: Path, annotation_file: Optional[Path] = None) -> Dict:
    """
    Analyze distribution of reads across RNA types

    Args:
        bam_file: Path to BAM file
        annotation_file: Optional GFF/GTF annotation file

    Returns:
        Dictionary with RNA type distribution
    """
    if not PYSAM_AVAILABLE:
        return {'status': 'error', 'error': 'pysam not installed'}

    # This is a placeholder - full implementation would use annotation file
    # For now, we'll categorize by read length as a proxy for RNA type

    length_categories = {
        'miRNA-like (18-23nt)': 0,
        'piRNA-like (24-32nt)': 0,
        'tRNA-fragment-like (33-40nt)': 0,
        'Other': 0
    }

    try:
        import pysam
        samfile = pysam.AlignmentFile(str(bam_file), "rb")

        for read in samfile.fetch(until_eof=True):
            if read.is_unmapped:
                continue

            length = read.query_length

            # Exclusive categorization (priority-based, no double-counting)
            if 18 <= length <= 23:
                length_categories['miRNA-like (18-23nt)'] += 1
            elif 24 <= length <= 32:
                length_categories['piRNA-like (24-32nt)'] += 1
            elif 33 <= length <= 40:
                length_categories['tRNA-fragment-like (33-40nt)'] += 1
            else:
                length_categories['Other'] += 1

        samfile.close()

        return {
            'status': 'success',
            'categories': length_categories
        }

    except FileNotFoundError:
        return {'status': 'error', 'error': f'BAM file not found: {bam_file}'}
    except ValueError as e:
        return {'status': 'error', 'error': f'Invalid BAM file: {e}'}
    except Exception as e:
        return {'status': 'error', 'error': f'Unexpected error analyzing {bam_file}: {e}'}


def render_post_alignment_qc_page():
    """Render the post-alignment QC page"""
    st.header("🔬 Post-Alignment QC")

    st.markdown("""
    Comprehensive quality assessment of aligned small RNA reads.
    This analysis helps identify potential issues with alignment and
    provides insights into the composition of your small RNA library.
    """)

    # Check for pysam
    if not check_pysam_installed():
        st.error("""
        ❌ **pysam not installed**. Please install it:
        ```bash
        pip install pysam
        # or
        conda install -c bioconda pysam
        ```
        """)
        return

    st.success("✅ pysam installed")

    # Check for alignment results
    bam_files = st.session_state.get('alignment_bam_files', [])

    tab1, tab2, tab3, tab4 = st.tabs([
        "📁 Input",
        "📊 Summary Statistics",
        "📈 Distributions",
        "🧬 RNA Composition"
    ])

    with tab1:
        render_post_qc_input(bam_files)

    with tab2:
        render_summary_stats()

    with tab3:
        render_distributions()

    with tab4:
        render_rna_composition()


def render_post_qc_input(alignment_bam_files: List):
    """Render input section for post-alignment QC"""
    st.subheader("Select BAM Files")

    # Check for files from alignment step
    if alignment_bam_files:
        st.success(f"✅ **{len(alignment_bam_files)} BAM files available from Alignment step**")

        input_source = st.radio(
            "Input Source",
            options=["Use alignment results", "Upload BAM files"],
            index=0,
            horizontal=True
        )

        if input_source == "Use alignment results":
            # Show available files
            file_info = []
            for bam_path in alignment_bam_files:
                bam_path = Path(bam_path)
                if bam_path.exists():
                    size_mb = bam_path.stat().st_size / (1024 * 1024)
                    file_info.append({
                        'File': bam_path.name,
                        'Size (MB)': f"{size_mb:.2f}",
                        'Status': '✅ Ready'
                    })

            if file_info:
                st.dataframe(pd.DataFrame(file_info), width="stretch", hide_index=True)

            st.session_state.post_qc_bam_files = [Path(f) for f in alignment_bam_files if Path(f).exists()]

            if st.button("🔍 Run Post-Alignment QC", type="primary"):
                run_post_alignment_qc()
            return

    # Manual upload
    st.info("Upload sorted, indexed BAM files for QC analysis")

    uploaded_files = st.file_uploader(
        "Upload BAM files",
        type=['bam'],
        accept_multiple_files=True
    )

    if uploaded_files:
        st.success(f"Uploaded {len(uploaded_files)} BAM files")

        # Save to temp directory
        temp_dir = Path(tempfile.mkdtemp())
        saved_files = []

        for uf in uploaded_files:
            file_path = temp_dir / uf.name
            with open(file_path, 'wb') as f:
                f.write(uf.getbuffer())
            saved_files.append(file_path)

            # Try to create index if not present
            try:
                subprocess.run(['samtools', 'index', str(file_path)],
                             capture_output=True, timeout=300)
            except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.CalledProcessError) as e:
                st.warning(f"Could not create BAM index: {e}")

        st.session_state.post_qc_bam_files = saved_files

        if st.button("🔍 Run Post-Alignment QC", type="primary"):
            run_post_alignment_qc()


def run_post_alignment_qc():
    """Run post-alignment QC analysis"""
    bam_files = st.session_state.get('post_qc_bam_files', [])

    if not bam_files:
        st.error("No BAM files available")
        return

    progress = st.progress(0)
    status = st.empty()

    all_results = {}

    for i, bam_file in enumerate(bam_files):
        bam_file = Path(bam_file)
        status.text(f"Analyzing {bam_file.name}...")
        progress.progress((i + 1) / len(bam_files))

        # Basic BAM analysis
        result = analyze_bam_basic(bam_file)
        result['filename'] = bam_file.name

        # Flagstat
        flagstat = get_flagstat(bam_file)
        result['flagstat'] = flagstat

        # RNA type distribution
        rna_dist = analyze_rna_type_distribution(bam_file)
        result['rna_distribution'] = rna_dist

        all_results[bam_file.name] = result

    progress.progress(1.0)
    status.text("Analysis complete!")

    st.session_state.post_qc_results = all_results
    st.success(f"✅ Analyzed {len(all_results)} BAM files")
    st.rerun()


def render_summary_stats():
    """Render summary statistics section"""
    st.subheader("📊 Summary Statistics")

    results = st.session_state.get('post_qc_results')

    if not results:
        st.info("Run Post-Alignment QC first to see results.")
        return

    # Summary table
    summary_data = []
    for name, result in results.items():
        if result.get('status') == 'success':
            summary_data.append({
                'Sample': name,
                'Total Reads': f"{result.get('total_reads', 0):,}",
                'Mapped': f"{result.get('mapped_reads', 0):,}",
                'Mapping %': f"{result.get('mapping_rate', 0):.1f}%",
                'Unique %': f"{result.get('unique_rate', 0):.1f}%",
                'Mean Length': f"{result.get('mean_length', 0):.1f}",
                'Mean MAPQ': f"{result.get('mean_mapq', 0):.1f}",
            })

    if summary_data:
        st.dataframe(pd.DataFrame(summary_data), width="stretch", hide_index=True)

    # Aggregate metrics
    if results:
        col1, col2, col3, col4 = st.columns(4)

        successful = [r for r in results.values() if r.get('status') == 'success']

        if successful:
            total_mapped = sum(r.get('mapped_reads', 0) for r in successful)
            avg_mapping = np.mean([r.get('mapping_rate', 0) for r in successful])
            avg_unique = np.mean([r.get('unique_rate', 0) for r in successful])
            avg_length = np.mean([r.get('mean_length', 0) for r in successful])

            with col1:
                st.metric("Total Mapped Reads", f"{total_mapped:,}")
            with col2:
                st.metric("Avg Mapping Rate", f"{avg_mapping:.1f}%")
            with col3:
                st.metric("Avg Unique Rate", f"{avg_unique:.1f}%")
            with col4:
                st.metric("Avg Read Length", f"{avg_length:.1f} nt")

    # Strand distribution
    st.subheader("Strand Distribution")

    strand_data = []
    for name, result in results.items():
        if result.get('status') == 'success':
            total = result.get('forward_strand', 0) + result.get('reverse_strand', 0)
            if total > 0:
                strand_data.append({
                    'Sample': name,
                    'Forward': result.get('forward_strand', 0),
                    'Reverse': result.get('reverse_strand', 0),
                    'Forward %': f"{100 * result.get('forward_strand', 0) / total:.1f}%",
                    'Reverse %': f"{100 * result.get('reverse_strand', 0) / total:.1f}%",
                })

    if strand_data:
        st.dataframe(pd.DataFrame(strand_data), width="stretch", hide_index=True)


def render_distributions():
    """Render distribution plots"""
    st.subheader("📈 Distributions")

    results = st.session_state.get('post_qc_results')

    if not results:
        st.info("Run Post-Alignment QC first to see results.")
        return

    import plotly.express as px
    import plotly.graph_objects as go

    # Collect all lengths for histogram
    all_lengths = []
    for name, result in results.items():
        if result.get('status') == 'success' and result.get('lengths'):
            all_lengths.extend(result['lengths'])

    if all_lengths:
        # Length distribution
        st.markdown("### Read Length Distribution (Aligned Reads)")

        fig = go.Figure()
        fig.add_trace(go.Histogram(
            x=all_lengths,
            nbinsx=40,
            name="All Samples",
            marker_color='#3498db'
        ))

        # Add RNA type reference regions
        fig.add_vrect(x0=18, x1=25, fillcolor='green', opacity=0.1,
                     annotation_text="miRNA", annotation_position="top left")
        fig.add_vrect(x0=24, x1=32, fillcolor='red', opacity=0.1,
                     annotation_text="piRNA", annotation_position="top right")

        fig.update_layout(
            title="Aligned Read Length Distribution",
            xaxis_title="Read Length (nt)",
            yaxis_title="Count",
            showlegend=False
        )
        st.plotly_chart(fig, width="stretch")

    # 5' Nucleotide bias (important for miRNA - should be enriched for U)
    st.markdown("### 5' Nucleotide Bias")
    st.caption("miRNAs typically show enrichment for Uracil (T in DNA) at the 5' end")

    five_prime_data = []
    for name, result in results.items():
        if result.get('status') == 'success' and result.get('five_prime_nt'):
            total = sum(result['five_prime_nt'].values())
            if total > 0:
                for nt, count in result['five_prime_nt'].items():
                    five_prime_data.append({
                        'Sample': name,
                        'Nucleotide': nt,
                        'Percentage': 100 * count / total
                    })

    if five_prime_data:
        df = pd.DataFrame(five_prime_data)
        fig = px.bar(df, x='Sample', y='Percentage', color='Nucleotide',
                    title="5' Nucleotide Distribution",
                    barmode='group')
        fig.update_layout(xaxis_tickangle=-45)
        st.plotly_chart(fig, width="stretch")

    # MAPQ distribution
    st.markdown("### Mapping Quality Distribution")

    all_mapq = []
    for name, result in results.items():
        if result.get('status') == 'success' and result.get('mapq_scores'):
            all_mapq.extend(result['mapq_scores'])

    if all_mapq:
        fig = px.histogram(x=all_mapq, nbins=50,
                          title="Mapping Quality Score Distribution",
                          labels={'x': 'MAPQ Score', 'y': 'Count'})
        st.plotly_chart(fig, width="stretch")


def render_rna_composition():
    """Render RNA composition analysis"""
    st.subheader("🧬 RNA Type Composition")

    results = st.session_state.get('post_qc_results')

    if not results:
        st.info("Run Post-Alignment QC first to see results.")
        return

    st.markdown("""
    RNA type classification based on read length. This is an approximation -
    for accurate classification, alignment to specific RNA databases is recommended.
    """)

    import plotly.express as px

    # Collect RNA distribution data
    rna_data = []
    for name, result in results.items():
        if result.get('status') == 'success':
            rna_dist = result.get('rna_distribution', {})
            if rna_dist.get('status') == 'success':
                categories = rna_dist.get('categories', {})
                for rna_type, count in categories.items():
                    rna_data.append({
                        'Sample': name,
                        'RNA Type': rna_type,
                        'Count': count
                    })

    if rna_data:
        df = pd.DataFrame(rna_data)

        # Stacked bar chart
        fig = px.bar(df, x='Sample', y='Count', color='RNA Type',
                    title="RNA Type Distribution by Sample",
                    barmode='stack')
        fig.update_layout(xaxis_tickangle=-45)
        st.plotly_chart(fig, width="stretch")

        # Pie chart for combined data
        combined = df.groupby('RNA Type')['Count'].sum().reset_index()
        fig2 = px.pie(combined, values='Count', names='RNA Type',
                     title="Overall RNA Type Composition")
        st.plotly_chart(fig2, width="stretch")

    # Reference distribution
    st.markdown("### Reference/Chromosome Distribution")

    ref_data = []
    for name, result in results.items():
        if result.get('status') == 'success' and result.get('references'):
            refs = result['references']
            total = sum(refs.values())
            # Top 10 references
            for ref, count in sorted(refs.items(), key=lambda x: x[1], reverse=True)[:10]:
                ref_data.append({
                    'Sample': name,
                    'Reference': ref,
                    'Reads': count,
                    'Percentage': 100 * count / total if total > 0 else 0
                })

    if ref_data:
        df = pd.DataFrame(ref_data)
        st.dataframe(df, width="stretch", hide_index=True)
