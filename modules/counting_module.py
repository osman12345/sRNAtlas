"""
Counting Module for sRNAtlas
Extracts read counts from BAM files and annotates with RNAcentral data
"""
import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
from typing import List, Dict, Optional, Tuple
import tempfile
import re

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from config import config
from utils.file_handlers import parse_fasta_annotations, categorize_rna
from utils.plotting import plot_rna_type_distribution
from utils.caching import hash_file
from utils.error_handling import ErrorClassifier, format_error_for_streamlit


def validate_bam_before_counting(bam_file: Path) -> Tuple[bool, str]:
    """
    Validate BAM file before attempting to read counts

    Args:
        bam_file: Path to BAM file

    Returns:
        Tuple of (is_valid, error_message)
    """
    import subprocess

    # Check file exists
    if not bam_file.exists():
        return False, f"BAM file not found: {bam_file}"

    # Check file has content
    if bam_file.stat().st_size == 0:
        return False, "BAM file is empty"

    # Check BAM magic bytes
    try:
        with open(bam_file, 'rb') as f:
            header = f.read(4)
            if header[:2] != b'\x1f\x8b':
                return False, "Invalid BAM file format (not BGZF compressed)"
    except Exception as e:
        return False, f"Cannot read BAM file: {e}"

    # Try samtools quickcheck if available
    try:
        result = subprocess.run(
            ['samtools', 'quickcheck', str(bam_file)],
            capture_output=True,
            text=True,
            timeout=30
        )
        if result.returncode != 0:
            return False, f"BAM file is corrupted: {result.stderr or 'quickcheck failed'}"
    except FileNotFoundError:
        pass  # samtools not available, skip this check
    except subprocess.TimeoutExpired:
        return False, "BAM validation timed out"

    # Try to actually read BAM with pysam to catch BGZF errors
    try:
        import pysam
        with pysam.AlignmentFile(str(bam_file), "rb") as bam:
            # Try to read first 100 reads to verify file is readable
            read_count = 0
            for read in bam:
                read_count += 1
                if read_count >= 100:
                    break
            if read_count == 0:
                return False, "BAM file contains no reads"
    except Exception as e:
        error_str = str(e)
        if 'BGZF' in error_str or 'bgzf' in error_str.lower():
            return False, f"BAM file is corrupted (BGZF error). Delete and re-run alignment."
        return False, f"Cannot read BAM file with pysam: {e}"

    return True, ""


# Counting mode constants
COUNTING_MODE_ALL = "all"  # Count all alignments (original behavior)
COUNTING_MODE_UNIQUE = "unique_only"  # Only unique mappers (MAPQ > threshold)
COUNTING_MODE_FRACTIONAL = "fractional"  # Fractional counts (1/n for n alignments)
COUNTING_MODE_PRIMARY = "primary_only"  # Only primary alignments (SAM flag)

COUNTING_MODES = {
    COUNTING_MODE_ALL: "All alignments (count multi-mappers fully)",
    COUNTING_MODE_UNIQUE: "Unique only (ignore multi-mappers)",
    COUNTING_MODE_FRACTIONAL: "Fractional (1/n weight for n alignments)",
    COUNTING_MODE_PRIMARY: "Primary only (use SAM primary flag)"
}


def count_reads_from_bam(
    bam_file: Path,
    min_mapq: int = 0,
    counting_mode: str = COUNTING_MODE_ALL
) -> Tuple[Dict, Dict]:
    """
    Extract read counts for each RNA from BAM file (cached wrapper)

    Args:
        bam_file: Path to sorted BAM file
        min_mapq: Minimum mapping quality (0 allows multi-mappers)
        counting_mode: How to handle multi-mappers:
            - "all": Count all alignments (inflates multi-mapper counts)
            - "unique_only": Only count reads with MAPQ > threshold
            - "fractional": Each read contributes 1/n to each of n alignments
            - "primary_only": Only count primary alignments (SAM flag)

    Returns:
        Tuple of (counts dict, stats dict)
    """
    file_hash = hash_file(bam_file)
    return _count_reads_from_bam_cached(file_hash, str(bam_file), min_mapq, counting_mode)


@st.cache_data(show_spinner=False, ttl=3600)
def _count_reads_from_bam_cached(
    _file_hash: str,
    bam_file_path: str,
    min_mapq: int = 0,
    counting_mode: str = COUNTING_MODE_ALL
) -> Tuple[Dict, Dict]:
    """Cached implementation of BAM read counting with multi-mapper handling"""
    import pysam

    bam_file = Path(bam_file_path)

    # Validate BAM file first
    is_valid, error_msg = validate_bam_before_counting(bam_file)
    if not is_valid:
        return None, {'error': error_msg}

    counts = defaultdict(float)  # Use float for fractional counting
    total_reads = 0
    mapped_reads = 0
    multi_mapped_reads = 0
    unique_reads = 0
    filtered_reads = 0

    try:
        # Open BAM file - use iteration instead of fetch() to avoid index requirement
        with pysam.AlignmentFile(str(bam_file), "rb") as bam:
            # Check if BAM has any references
            if bam.nreferences == 0:
                return None, {'error': 'BAM file has no reference sequences. Check alignment reference.'}

            # Iterate through all reads (doesn't require index)
            for read in bam:
                total_reads += 1

                if read.is_unmapped:
                    continue

                mapped_reads += 1

                # Get NH tag (number of hits) if available
                try:
                    nh = read.get_tag('NH')
                except KeyError:
                    nh = 1  # Assume unique if NH tag not present

                is_multi_mapper = nh > 1
                if is_multi_mapper:
                    multi_mapped_reads += 1
                else:
                    unique_reads += 1

                # Apply counting mode logic
                if counting_mode == COUNTING_MODE_UNIQUE:
                    # Skip multi-mappers entirely
                    if is_multi_mapper or read.mapping_quality < min_mapq:
                        filtered_reads += 1
                        continue
                    weight = 1.0

                elif counting_mode == COUNTING_MODE_FRACTIONAL:
                    # Each alignment gets 1/NH weight
                    if read.mapping_quality < min_mapq:
                        filtered_reads += 1
                        continue
                    weight = 1.0 / nh

                elif counting_mode == COUNTING_MODE_PRIMARY:
                    # Only count primary alignments
                    if read.is_secondary or read.is_supplementary:
                        filtered_reads += 1
                        continue
                    if read.mapping_quality < min_mapq:
                        filtered_reads += 1
                        continue
                    weight = 1.0

                else:  # COUNTING_MODE_ALL (default)
                    # Count all alignments (original behavior)
                    if read.mapping_quality < min_mapq:
                        filtered_reads += 1
                        continue
                    weight = 1.0

                # Get reference name (RNA ID from FASTA)
                if read.reference_id >= 0:
                    ref_name = bam.get_reference_name(read.reference_id)
                    counts[ref_name] += weight

        # Check if we got any data
        if total_reads == 0:
            return None, {'error': 'BAM file contains no reads'}

        if mapped_reads == 0:
            return None, {'error': f'BAM file has {total_reads} reads but none are mapped'}

        if len(counts) == 0:
            return None, {'error': f'No reads passed filters (mapped={mapped_reads}, min_mapq={min_mapq}, mode={counting_mode})'}

        # Convert to int for non-fractional modes, keep float for fractional
        if counting_mode != COUNTING_MODE_FRACTIONAL:
            counts = {k: int(v) for k, v in counts.items()}

        stats = {
            'total_reads': total_reads,
            'mapped_reads': mapped_reads,
            'unmapped_reads': total_reads - mapped_reads,
            'unique_reads': unique_reads,
            'multi_mapped_reads': multi_mapped_reads,
            'filtered_reads': filtered_reads,
            'unique_rnas': len(counts),
            'alignment_rate': 100 * mapped_reads / total_reads if total_reads > 0 else 0,
            'multi_mapping_rate': 100 * multi_mapped_reads / mapped_reads if mapped_reads > 0 else 0,
            'counting_mode': counting_mode
        }

        return dict(counts), stats

    except Exception as e:
        error_msg = str(e)
        # Provide more helpful error messages
        if 'BGZF' in error_msg or 'bgzf' in error_msg.lower():
            error_msg = f"BAM file is corrupted or truncated. Try re-running alignment. Original error: {e}"
        elif 'truncated' in error_msg.lower():
            error_msg = f"BAM file appears to be truncated. Try re-running alignment. Original error: {e}"
        elif 'no index' in error_msg.lower():
            error_msg = f"BAM index missing. Re-run alignment to generate index. Original error: {e}"
        return None, {'error': error_msg}


def process_bam_files(
    bam_files: List[Path],
    min_mapq: int = 0,
    counting_mode: str = COUNTING_MODE_ALL,
    progress_callback=None
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Process multiple BAM files and create count matrix

    Args:
        bam_files: List of BAM file paths
        min_mapq: Minimum mapping quality threshold
        counting_mode: Multi-mapper handling mode
        progress_callback: Optional callback for progress updates

    Returns:
        Tuple of (count_matrix, sample_stats)
    """
    all_counts = {}
    all_stats = {}

    for i, bam_file in enumerate(bam_files):
        # Extract sample name
        sample_name = bam_file.stem.replace('_sorted', '').replace('_rnacentral', '')

        if progress_callback:
            progress_callback(i / len(bam_files), f"Processing {sample_name}...")

        counts, stats = count_reads_from_bam(bam_file, min_mapq, counting_mode)

        if counts is not None:
            all_counts[sample_name] = counts
            all_stats[sample_name] = stats

    if not all_counts:
        return None, None

    # Create count matrix
    count_df = pd.DataFrame(all_counts).fillna(0).astype(int)

    # Sort rows by total counts (most abundant first)
    count_df['_total'] = count_df.sum(axis=1)
    count_df = count_df.sort_values('_total', ascending=False)
    count_df = count_df.drop('_total', axis=1)

    # Sort columns (samples) alphabetically
    count_df = count_df[sorted(count_df.columns)]

    # Create stats dataframe
    stats_df = pd.DataFrame(all_stats).T
    stats_df.index.name = 'sample'

    return count_df, stats_df


def annotate_count_matrix(
    count_matrix: pd.DataFrame,
    annotations: pd.DataFrame
) -> pd.DataFrame:
    """
    Merge count matrix with RNA annotations

    Args:
        count_matrix: Count matrix with RNA IDs as index
        annotations: Annotation dataframe

    Returns:
        Annotated count matrix
    """
    # Reset index for merging
    counts_reset = count_matrix.reset_index()
    counts_reset.columns = ['full_id'] + list(count_matrix.columns)

    # Merge with annotations
    merged = counts_reset.merge(
        annotations,
        on='full_id',
        how='left'
    )

    # Set index back to full_id
    merged = merged.set_index('full_id')

    return merged


def filter_low_counts(
    count_matrix: pd.DataFrame,
    min_cpm: float = 0.5,
    min_samples: int = 2,
    sample_cols: Optional[List[str]] = None
) -> Tuple[pd.DataFrame, Dict]:
    """
    Filter low-expressed features

    Args:
        count_matrix: Count matrix
        min_cpm: Minimum CPM threshold
        min_samples: Minimum number of samples above threshold
        sample_cols: List of sample column names

    Returns:
        Tuple of (filtered_matrix, filter_stats)
    """
    if sample_cols is None:
        # Auto-detect sample columns (numeric only)
        sample_cols = count_matrix.select_dtypes(include=[np.number]).columns.tolist()

    # Calculate CPM
    lib_sizes = count_matrix[sample_cols].sum()
    cpm = count_matrix[sample_cols].div(lib_sizes) * 1e6

    # Apply filter
    keep = (cpm >= min_cpm).sum(axis=1) >= min_samples

    filtered = count_matrix[keep].copy()

    stats = {
        'before': len(count_matrix),
        'after': len(filtered),
        'removed': len(count_matrix) - len(filtered),
        'min_cpm': min_cpm,
        'min_samples': min_samples
    }

    return filtered, stats


def render_counting_page():
    """Render the counting page"""
    st.header("📈 Read Counting & Annotation")

    tab1, tab2, tab3, tab4 = st.tabs([
        "📁 Input Data",
        "🔢 Count Matrix",
        "🏷️ Annotation",
        "🔍 Filtering"
    ])

    with tab1:
        render_counting_input()

    with tab2:
        render_count_matrix()

    with tab3:
        render_annotation()

    with tab4:
        render_filtering()


def render_counting_input():
    """Render input data section"""
    st.subheader("Input Data")

    input_type = st.radio(
        "Input Type",
        options=["BAM files", "Count matrix (CSV)"]
    )

    if input_type == "BAM files":
        # Check for pysam
        try:
            import pysam
            st.success("✅ pysam available")
        except ImportError:
            st.error("❌ pysam not installed. Install with: pip install pysam")
            return

        # Check if BAM files are available from alignment step
        alignment_bams = st.session_state.get('alignment_bam_files')
        alignment_dir = st.session_state.get('alignment_output_dir')
        bams_from_alignment = alignment_bams and len(alignment_bams) > 0

        bam_files = None
        use_alignment_bams = False

        if bams_from_alignment:
            st.success(f"✅ {len(alignment_bams)} BAM files available from Alignment step")

            input_source = st.radio(
                "Input Source",
                options=["Use BAMs from Alignment", "Upload new BAM files"],
                index=0,
                help="Use BAM files from previous alignment, or upload new ones"
            )

            if input_source == "Use BAMs from Alignment":
                st.info(f"Using BAM files from: {alignment_dir or 'alignment output'}")

                # Show file list
                file_info = [{'Filename': Path(f).name, 'Source': 'Alignment Step'} for f in alignment_bams]
                st.dataframe(pd.DataFrame(file_info), width="stretch")

                use_alignment_bams = True
            else:
                bam_files = st.file_uploader(
                    "Upload BAM files",
                    type=['bam'],
                    accept_multiple_files=True
                )
        else:
            bam_files = st.file_uploader(
                "Upload BAM files",
                type=['bam'],
                accept_multiple_files=True
            )

        if bam_files:
            st.info(f"Uploaded {len(bam_files)} BAM files")

        has_bams = bam_files or use_alignment_bams

        if has_bams:
            st.subheader("⚙️ Counting Settings")

            # Counting mode selector
            counting_mode = st.selectbox(
                "Multi-mapper Handling",
                options=list(COUNTING_MODES.keys()),
                format_func=lambda x: COUNTING_MODES[x],
                index=0,
                help="""
                **All alignments**: Count every alignment (inflates counts for multi-mappers - original behavior)

                **Unique only**: Only count uniquely mapped reads (MAPQ > threshold)

                **Fractional**: Each read contributes 1/n to each of its n alignments (recommended for miRNA families)

                **Primary only**: Only count primary alignments (uses SAM flag)
                """
            )

            # Store in session state for reports
            st.session_state.counting_mode = counting_mode

            min_mapq = st.slider(
                "Minimum Mapping Quality",
                0, 60, config.counting.min_mapq,
                help="0 = include all alignments, higher values filter low-confidence mappings"
            )

            # Show mode-specific info
            if counting_mode == COUNTING_MODE_UNIQUE:
                st.info("💡 **Unique mode**: Multi-mappers will be excluded. Consider increasing MAPQ threshold.")
            elif counting_mode == COUNTING_MODE_FRACTIONAL:
                st.info("💡 **Fractional mode**: A read mapping to 3 locations will contribute 0.33 counts to each.")
            elif counting_mode == COUNTING_MODE_PRIMARY:
                st.info("💡 **Primary mode**: Only the best alignment per read will be counted.")

            if st.button("🔢 Generate Count Matrix", type="primary"):
                # Get BAM file paths
                bam_paths = []

                if use_alignment_bams:
                    # Use BAMs from alignment step
                    bam_paths = [Path(f) for f in alignment_bams]
                else:
                    # Save uploaded files temporarily
                    temp_dir = Path(tempfile.mkdtemp())
                    for bf in bam_files:
                        file_path = temp_dir / bf.name
                        with open(file_path, 'wb') as f:
                            f.write(bf.getbuffer())
                        bam_paths.append(file_path)

                # Process BAM files
                progress_bar = st.progress(0)
                status_text = st.empty()

                def update_progress(progress, message):
                    progress_bar.progress(progress)
                    status_text.text(message)

                with st.spinner("Counting reads..."):
                    count_matrix, sample_stats = process_bam_files(
                        bam_paths,
                        min_mapq,
                        counting_mode,
                        progress_callback=update_progress
                    )

                if count_matrix is not None and not count_matrix.empty:
                    st.session_state.count_matrix = count_matrix
                    st.session_state.sample_stats = sample_stats

                    # Show mode-specific success message
                    mode_msg = COUNTING_MODES.get(counting_mode, counting_mode)
                    st.success(f"Generated count matrix: {count_matrix.shape[0]} RNAs × {count_matrix.shape[1]} samples")
                    st.caption(f"Counting mode: {mode_msg}")
                else:
                    st.error("Failed to generate count matrix.")

                    # Try to show specific errors for each BAM with actionable diagnostics
                    st.markdown("**BAM file status:**")
                    for bam_path in bam_paths:
                        is_valid, error_msg = validate_bam_before_counting(bam_path)
                        if is_valid:
                            st.markdown(f"- ✅ {bam_path.name}: OK")
                        else:
                            st.markdown(f"- ❌ **{bam_path.name}**")
                            # Get actionable error info
                            diag_error = ErrorClassifier.classify_bam_error(error_msg, bam_path)
                            st.markdown(format_error_for_streamlit(diag_error))

    else:  # Count matrix
        count_file = st.file_uploader(
            "Upload count matrix (CSV)",
            type=['csv', 'tsv', 'txt']
        )

        if count_file:
            # Detect delimiter
            content = count_file.getvalue().decode('utf-8')
            sep = '\t' if '\t' in content.split('\n')[0] else ','

            count_file.seek(0)
            try:
                counts = pd.read_csv(count_file, sep=sep, index_col=0)
                st.session_state.count_matrix = counts
                st.success(f"Loaded count matrix: {counts.shape[0]} features × {counts.shape[1]} columns")

                # Show preview
                st.dataframe(counts.head(10), width="stretch")

            except Exception as e:
                st.error(f"Error loading file: {e}")


def render_count_matrix():
    """Render count matrix visualization"""
    st.subheader("Count Matrix")

    if st.session_state.get('count_matrix') is None:
        st.info("Please load data in the 'Input Data' tab first.")
        return

    counts = st.session_state.count_matrix

    # Identify sample columns
    sample_cols = counts.select_dtypes(include=[np.number]).columns.tolist()

    # Summary stats
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric("Total Features", f"{len(counts):,}")

    with col2:
        st.metric("Samples", len(sample_cols))

    with col3:
        total_reads = counts[sample_cols].sum().sum()
        st.metric("Total Reads", f"{total_reads:,.0f}")

    with col4:
        detected = (counts[sample_cols].sum(axis=1) > 0).sum()
        st.metric("Detected Features", f"{detected:,}")

    # Sample statistics
    if st.session_state.get('sample_stats') is not None:
        st.subheader("Sample Statistics")
        st.dataframe(st.session_state.sample_stats, width="stretch")

    # Top features
    st.subheader("Top Abundant Features")

    top_n = st.slider("Number of features to show", 10, 100, 20)
    top_features = counts[sample_cols].sum(axis=1).nlargest(top_n)

    import plotly.express as px

    fig = px.bar(
        x=top_features.index[::-1],
        y=top_features.values[::-1],
        orientation='h',
        title=f"Top {top_n} Most Abundant Features",
        labels={'x': 'RNA ID', 'y': 'Total Reads'}
    )
    fig.update_layout(height=max(400, top_n * 20))
    st.plotly_chart(fig, width="stretch")

    # Download
    st.divider()
    csv = counts.to_csv()
    st.download_button(
        "📥 Download Count Matrix",
        csv,
        "count_matrix.csv",
        "text/csv"
    )


def render_annotation():
    """Render annotation section"""
    st.subheader("RNA Annotation")

    if st.session_state.get('count_matrix') is None:
        st.info("Please load count data first.")
        return

    counts = st.session_state.count_matrix

    annotation_option = st.radio(
        "Annotation Source",
        options=["Upload annotation file", "Parse from FASTA headers"]
    )

    if annotation_option == "Upload annotation file":
        annot_file = st.file_uploader(
            "Upload annotation CSV",
            type=['csv', 'tsv', 'txt']
        )

        if annot_file:
            # Detect delimiter
            content = annot_file.getvalue().decode('utf-8')
            sep = '\t' if '\t' in content.split('\n')[0] else ','

            annot_file.seek(0)
            try:
                annotations = pd.read_csv(annot_file, sep=sep)

                # Check for required columns
                required = ['full_id', 'RNA_category']
                if not all(col in annotations.columns for col in required):
                    # Try to set index column as full_id
                    if annotations.index.name:
                        annotations = annotations.reset_index()
                        annotations.columns = ['full_id'] + list(annotations.columns[1:])

                st.session_state.annotations = annotations
                st.success(f"Loaded {len(annotations)} annotations")

                # Preview
                st.dataframe(annotations.head(10), width="stretch")

            except Exception as e:
                st.error(f"Error loading annotations: {e}")

    else:  # Parse from FASTA headers
        fasta_file = st.file_uploader(
            "Upload Reference FASTA",
            type=['fasta', 'fa', 'fna', 'gz']
        )

        if fasta_file:
            if st.button("🏷️ Parse Annotations"):
                with st.spinner("Parsing FASTA headers..."):
                    # Save temporarily
                    temp_dir = Path(tempfile.mkdtemp())
                    fasta_path = temp_dir / fasta_file.name
                    with open(fasta_path, 'wb') as f:
                        f.write(fasta_file.getbuffer())

                    # Parse
                    annotations = parse_fasta_annotations(fasta_path)

                    st.session_state.annotations = annotations
                    st.success(f"Parsed {len(annotations)} annotations")

    # Merge annotations with counts
    if st.session_state.get('annotations') is not None:
        annotations = st.session_state.annotations

        st.divider()
        st.subheader("Merge Annotations with Counts")

        if st.button("🔗 Merge"):
            with st.spinner("Merging..."):
                annotated = annotate_count_matrix(counts, annotations)
                st.session_state.annotated_counts = annotated
                st.success(f"Created annotated matrix: {annotated.shape}")

        if st.session_state.get('annotated_counts') is not None:
            annotated = st.session_state.annotated_counts

            # Show RNA type distribution
            if 'RNA_category' in annotated.columns:
                st.subheader("RNA Type Distribution")

                sample_cols = annotated.select_dtypes(include=[np.number]).columns.tolist()
                fig = plot_rna_type_distribution(annotated, annotated[sample_cols])
                st.plotly_chart(fig, width="stretch")

                # Type counts table
                type_counts = annotated['RNA_category'].value_counts()
                st.dataframe(type_counts, width="stretch")

            # Download
            csv = annotated.to_csv()
            st.download_button(
                "📥 Download Annotated Counts",
                csv,
                "annotated_counts.csv",
                "text/csv"
            )


def render_filtering():
    """Render filtering section"""
    st.subheader("Low-Count Filtering")

    counts = st.session_state.get('annotated_counts') or st.session_state.get('count_matrix')

    if counts is None:
        st.info("Please load count data first.")
        return

    sample_cols = counts.select_dtypes(include=[np.number]).columns.tolist()

    st.markdown("""
    Filter out lowly-expressed features to improve statistical power in downstream analysis.
    The default settings require at least 0.5 CPM in at least 2 samples.
    """)

    col1, col2 = st.columns(2)

    with col1:
        min_cpm = st.slider(
            "Minimum CPM",
            0.0, 10.0, config.counting.min_cpm,
            step=0.1,
            help="Minimum counts per million threshold"
        )

    with col2:
        # Handle edge case when there's only 1 sample
        if len(sample_cols) <= 1:
            min_samples = 1
            st.info("Only 1 sample - minimum samples set to 1")
        else:
            min_samples = st.slider(
                "Minimum Samples",
                1, len(sample_cols), min(config.counting.min_samples_detected, len(sample_cols)),
                help="Feature must pass threshold in at least this many samples"
            )

    # Preview filter effect
    st.subheader("Filter Preview")

    # Calculate CPM
    lib_sizes = counts[sample_cols].sum()
    cpm = counts[sample_cols].div(lib_sizes) * 1e6

    # Apply filter
    keep = (cpm >= min_cpm).sum(axis=1) >= min_samples

    col1, col2, col3 = st.columns(3)

    with col1:
        st.metric("Before Filtering", len(counts))

    with col2:
        st.metric("After Filtering", keep.sum())

    with col3:
        removed = len(counts) - keep.sum()
        st.metric("Removed", removed, delta=f"-{100*removed/len(counts):.1f}%")

    # Show distribution of features being removed
    import plotly.express as px

    kept_total = counts.loc[keep, sample_cols].sum(axis=1)
    removed_total = counts.loc[~keep, sample_cols].sum(axis=1)

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Kept Features**")
        if len(kept_total) > 0:
            fig = px.histogram(kept_total, nbins=50, title="Total counts distribution (kept)")
            fig.update_layout(showlegend=False)
            st.plotly_chart(fig, width="stretch")

    with col2:
        st.markdown("**Removed Features**")
        if len(removed_total) > 0:
            fig = px.histogram(removed_total, nbins=50, title="Total counts distribution (removed)")
            fig.update_layout(showlegend=False)
            st.plotly_chart(fig, width="stretch")

    # Apply filter
    st.divider()
    if st.button("✅ Apply Filter", type="primary"):
        filtered_counts = counts[keep].copy()

        if st.session_state.get('annotated_counts') is not None:
            st.session_state.annotated_counts = filtered_counts
        else:
            st.session_state.count_matrix = filtered_counts

        st.success(f"Filter applied! {len(filtered_counts)} features remaining.")

        # Download
        csv = filtered_counts.to_csv()
        st.download_button(
            "📥 Download Filtered Counts",
            csv,
            "filtered_counts.csv",
            "text/csv"
        )
