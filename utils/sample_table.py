"""
Sample Table View for sRNAtlas
Unified view of all samples and their processing status
"""
import streamlit as st
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass, field


@dataclass
class SampleStatus:
    """Track processing status for a single sample"""
    name: str
    fastq_file: str = ""
    trimmed_file: str = ""
    bam_file: str = ""
    
    # Metrics
    raw_reads: int = 0
    trimmed_reads: int = 0
    aligned_reads: int = 0
    alignment_rate: float = 0.0
    
    # Status flags
    qc_done: bool = False
    trimming_done: bool = False
    alignment_done: bool = False
    counting_done: bool = False
    
    # Condition info
    condition: str = ""
    batch: str = ""


def init_sample_tracker():
    """Initialize sample tracker in session state"""
    if 'sample_status' not in st.session_state:
        st.session_state.sample_status: Dict[str, SampleStatus] = {}


def add_sample(name: str, fastq_file: str = "", condition: str = "", batch: str = ""):
    """Add a new sample to tracking"""
    init_sample_tracker()
    
    if name not in st.session_state.sample_status:
        st.session_state.sample_status[name] = SampleStatus(
            name=name,
            fastq_file=fastq_file,
            condition=condition,
            batch=batch
        )


def update_sample(name: str, **kwargs):
    """Update sample information"""
    init_sample_tracker()
    
    if name in st.session_state.sample_status:
        sample = st.session_state.sample_status[name]
        for key, value in kwargs.items():
            if hasattr(sample, key):
                setattr(sample, key, value)


def get_sample_dataframe() -> pd.DataFrame:
    """Get all samples as a DataFrame"""
    init_sample_tracker()
    
    if not st.session_state.sample_status:
        return pd.DataFrame()
    
    data = []
    for sample in st.session_state.sample_status.values():
        data.append({
            'Sample': sample.name,
            'Condition': sample.condition or '-',
            'Batch': sample.batch or '-',
            'Raw Reads': f"{sample.raw_reads:,}" if sample.raw_reads else '-',
            'Trimmed': f"{sample.trimmed_reads:,}" if sample.trimmed_reads else '-',
            'Aligned': f"{sample.aligned_reads:,}" if sample.aligned_reads else '-',
            'Align %': f"{sample.alignment_rate:.1f}%" if sample.alignment_rate else '-',
            'QC': '✅' if sample.qc_done else '⬜',
            'Trim': '✅' if sample.trimming_done else '⬜',
            'Align': '✅' if sample.alignment_done else '⬜',
            'Count': '✅' if sample.counting_done else '⬜',
        })
    
    return pd.DataFrame(data)


def render_sample_table():
    """Render the sample status table"""
    init_sample_tracker()
    
    st.subheader("📋 Sample Overview")
    
    if not st.session_state.sample_status:
        st.info("No samples loaded. Upload FASTQ files in the Project module.")
        return
    
    # Summary metrics
    total = len(st.session_state.sample_status)
    qc_done = sum(1 for s in st.session_state.sample_status.values() if s.qc_done)
    trim_done = sum(1 for s in st.session_state.sample_status.values() if s.trimming_done)
    align_done = sum(1 for s in st.session_state.sample_status.values() if s.alignment_done)
    
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Total Samples", total)
    col2.metric("QC Complete", f"{qc_done}/{total}")
    col3.metric("Trimmed", f"{trim_done}/{total}")
    col4.metric("Aligned", f"{align_done}/{total}")
    
    # Sample table
    df = get_sample_dataframe()
    if not df.empty:
        st.dataframe(
            df,
            use_container_width=True,
            hide_index=True,
            column_config={
                'Sample': st.column_config.TextColumn('Sample', width='medium'),
                'Condition': st.column_config.TextColumn('Condition', width='small'),
                'Batch': st.column_config.TextColumn('Batch', width='small'),
                'Raw Reads': st.column_config.TextColumn('Raw Reads', width='small'),
                'Trimmed': st.column_config.TextColumn('Trimmed', width='small'),
                'Aligned': st.column_config.TextColumn('Aligned', width='small'),
                'Align %': st.column_config.TextColumn('Align %', width='small'),
                'QC': st.column_config.TextColumn('QC', width='small'),
                'Trim': st.column_config.TextColumn('Trim', width='small'),
                'Align': st.column_config.TextColumn('Align', width='small'),
                'Count': st.column_config.TextColumn('Count', width='small'),
            }
        )


def render_compact_sample_list():
    """Render compact sample list for sidebar"""
    init_sample_tracker()
    
    if not st.session_state.sample_status:
        st.caption("No samples loaded")
        return
    
    total = len(st.session_state.sample_status)
    conditions = set(s.condition for s in st.session_state.sample_status.values() if s.condition)
    
    st.caption(f"**Samples:** {total}")
    if conditions:
        st.caption(f"**Conditions:** {', '.join(conditions)}")


def export_sample_report() -> str:
    """Export sample status as CSV string"""
    df = get_sample_dataframe()
    return df.to_csv(index=False)


def clear_samples():
    """Clear all sample tracking"""
    if 'sample_status' in st.session_state:
        st.session_state.sample_status = {}
