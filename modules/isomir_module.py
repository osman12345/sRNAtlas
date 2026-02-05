"""
isomiR Analysis Module for sRNAtlas
Detect and analyze miRNA isoforms (isomiRs)
"""
import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import defaultdict
import plotly.express as px
import plotly.graph_objects as go

try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False


def render_isomir_page():
    """Render the isomiR analysis page"""
    st.header("🔬 isomiR Analysis")
    st.markdown("Detect and analyze miRNA isoforms (sequence variants)")
    
    # Info box
    with st.expander("ℹ️ What are isomiRs?", expanded=False):
        st.markdown("""
        **isomiRs** are sequence variants of canonical miRNAs that differ by:
        
        - **5' variants**: Different start position (affects target specificity!)
        - **3' variants**: Different end position (most common)
        - **SNPs**: Internal nucleotide changes
        - **Non-templated additions**: Added nucleotides (often A or U)
        
        isomiRs can have different biological functions and target different genes.
        """)
    
    tab1, tab2, tab3, tab4 = st.tabs([
        "📥 Input",
        "🔬 Analysis",
        "📊 Results",
        "📈 Visualization"
    ])
    
    with tab1:
        render_input_tab()
    
    with tab2:
        render_analysis_tab()
    
    with tab3:
        render_results_tab()
    
    with tab4:
        render_visualization_tab()


def render_input_tab():
    """Input selection for isomiR analysis"""
    st.subheader("Input Data")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("#### Aligned BAM Files")
        
        # Check for BAM files from alignment
        if 'alignment_results' in st.session_state and st.session_state.alignment_results:
            bam_files = []
            for result in st.session_state.alignment_results:
                if 'bam_file' in result and Path(result['bam_file']).exists():
                    bam_files.append(result['bam_file'])
            
            if bam_files:
                st.success(f"Found {len(bam_files)} BAM files from alignment")
                use_bam = st.checkbox("Use BAM files from alignment", value=True)
                if use_bam:
                    st.session_state.isomir_bam_files = bam_files
        else:
            st.info("Run alignment first to generate BAM files")
    
    with col2:
        st.markdown("#### Reference miRNAs")
        
        # Check for reference from database module
        if 'reference_fasta' in st.session_state and st.session_state.reference_fasta:
            st.success("✅ Reference miRNAs available")
            st.session_state.isomir_reference = st.session_state.reference_fasta
        else:
            st.warning("Load reference in Databases module first")


def render_analysis_tab():
    """Run isomiR analysis"""
    st.subheader("Run isomiR Detection")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("#### Detection Settings")
        
        min_reads = st.number_input(
            "Minimum reads per isomiR",
            min_value=1, max_value=100, value=10
        )
        st.session_state.isomir_min_reads = min_reads
        
        max_5p_diff = st.slider(
            "Max 5' position difference",
            0, 5, 2, help="Nucleotides from canonical start"
        )
        st.session_state.isomir_max_5p = max_5p_diff
        
        max_3p_diff = st.slider(
            "Max 3' position difference", 
            0, 10, 4, help="Nucleotides from canonical end"
        )
        st.session_state.isomir_max_3p = max_3p_diff
    
    with col2:
        st.markdown("#### Variant Types")
        
        detect_5p = st.checkbox("Detect 5' variants", value=True)
        detect_3p = st.checkbox("Detect 3' variants", value=True)
        detect_snp = st.checkbox("Detect internal SNPs", value=True)
        detect_nta = st.checkbox("Detect non-templated additions", value=True)
        
        st.session_state.isomir_detect = {
            '5p': detect_5p,
            '3p': detect_3p,
            'snp': detect_snp,
            'nta': detect_nta
        }
    
    st.markdown("---")
    
    # Check prerequisites
    bam_files = st.session_state.get('isomir_bam_files', [])
    reference = st.session_state.get('isomir_reference')
    
    ready = True
    if not bam_files:
        st.error("❌ No BAM files available")
        ready = False
    if not reference:
        st.error("❌ No reference miRNAs loaded")
        ready = False
    if not PYSAM_AVAILABLE:
        st.error("❌ pysam not installed")
        ready = False
    
    if ready:
        if st.button("🚀 Detect isomiRs", type="primary"):
            with st.spinner("Analyzing isomiRs..."):
                results = detect_isomirs(bam_files, reference)
                st.session_state.isomir_results = results
                st.success("✅ Analysis complete!")
                st.rerun()


def detect_isomirs(bam_files: List[str], reference: str) -> Dict:
    """Detect isomiRs from BAM files"""
    results = {
        'isomirs': [],
        'stats': {},
        'by_mirna': defaultdict(list),
        'variant_counts': defaultdict(int)
    }
    
    min_reads = st.session_state.get('isomir_min_reads', 10)
    
    # Track all read sequences per reference
    ref_reads = defaultdict(lambda: defaultdict(int))
    
    for bam_file in bam_files:
        try:
            with pysam.AlignmentFile(bam_file, "rb") as bam:
                for read in bam.fetch():
                    if read.is_unmapped:
                        continue
                    
                    ref_name = read.reference_name
                    seq = read.query_sequence
                    
                    # Count this sequence for this reference
                    ref_reads[ref_name][seq] += 1
        except Exception as e:
            st.warning(f"Error reading {bam_file}: {e}")
    
    # Analyze variants for each miRNA
    for mirna, seq_counts in ref_reads.items():
        # Find canonical (most abundant)
        canonical_seq = max(seq_counts.keys(), key=lambda x: seq_counts[x])
        canonical_count = seq_counts[canonical_seq]
        
        for seq, count in seq_counts.items():
            if count < min_reads:
                continue
            
            # Classify variant type
            variant_type = classify_variant(canonical_seq, seq)
            
            isomir = {
                'mirna': mirna,
                'sequence': seq,
                'length': len(seq),
                'count': count,
                'canonical': seq == canonical_seq,
                'variant_type': variant_type,
                'relative_abundance': count / canonical_count * 100
            }
            
            results['isomirs'].append(isomir)
            results['by_mirna'][mirna].append(isomir)
            results['variant_counts'][variant_type] += 1
    
    # Calculate stats
    total_isomirs = len(results['isomirs'])
    canonical_count = sum(1 for i in results['isomirs'] if i['canonical'])
    
    results['stats'] = {
        'total_isomirs': total_isomirs,
        'canonical': canonical_count,
        'variants': total_isomirs - canonical_count,
        'mirnas_with_variants': len([m for m, v in results['by_mirna'].items() if len(v) > 1])
    }
    
    return results


def classify_variant(canonical: str, variant: str) -> str:
    """Classify the type of isomiR variant"""
    if canonical == variant:
        return "canonical"
    
    # Check for length differences (5'/3' variants)
    len_diff = len(variant) - len(canonical)
    
    if len_diff != 0:
        # Check if it's a 5' or 3' variant
        if variant in canonical or canonical in variant:
            if len_diff > 0:
                return "3p_addition"
            else:
                return "3p_trimming"
    
    # Check for internal mismatches (SNPs)
    if len(canonical) == len(variant):
        mismatches = sum(1 for a, b in zip(canonical, variant) if a != b)
        if mismatches == 1:
            return "snp"
        elif mismatches > 1:
            return "multiple_snp"
    
    # Check for non-templated additions
    if len(variant) > len(canonical):
        if variant.startswith(canonical):
            added = variant[len(canonical):]
            if all(nt in ['A', 'T'] for nt in added):
                return "nta"
    
    return "other"


def render_results_tab():
    """Display isomiR results"""
    st.subheader("isomiR Detection Results")
    
    if 'isomir_results' not in st.session_state:
        st.info("Run analysis to see results")
        return
    
    results = st.session_state.isomir_results
    
    # Summary stats
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Total isomiRs", results['stats']['total_isomirs'])
    col2.metric("Canonical", results['stats']['canonical'])
    col3.metric("Variants", results['stats']['variants'])
    col4.metric("miRNAs with variants", results['stats']['mirnas_with_variants'])
    
    st.markdown("---")
    
    # Variant type breakdown
    st.markdown("#### Variant Types")
    variant_df = pd.DataFrame([
        {'Type': k, 'Count': v}
        for k, v in results['variant_counts'].items()
    ])
    if not variant_df.empty:
        fig = px.pie(variant_df, values='Count', names='Type', title="isomiR Variant Types")
        st.plotly_chart(fig, width="stretch")
    
    st.markdown("---")
    
    # Full results table
    st.markdown("#### All Detected isomiRs")
    
    if results['isomirs']:
        df = pd.DataFrame(results['isomirs'])
        df = df[['mirna', 'sequence', 'length', 'count', 'canonical', 'variant_type', 'relative_abundance']]
        df.columns = ['miRNA', 'Sequence', 'Length', 'Reads', 'Canonical', 'Variant Type', 'Rel. Abundance %']
        df['Rel. Abundance %'] = df['Rel. Abundance %'].round(1)
        df['Canonical'] = df['Canonical'].map({True: '✅', False: ''})
        
        # Filter options
        mirna_filter = st.selectbox(
            "Filter by miRNA",
            ['All'] + sorted(df['miRNA'].unique().tolist())
        )
        
        if mirna_filter != 'All':
            df = df[df['miRNA'] == mirna_filter]
        
        st.dataframe(df, width="stretch", hide_index=True)
        
        # Download
        csv = df.to_csv(index=False)
        st.download_button(
            "📥 Download isomiR Table (CSV)",
            csv,
            file_name="isomir_results.csv",
            mime="text/csv"
        )


def render_visualization_tab():
    """Visualize isomiR results"""
    st.subheader("isomiR Visualization")
    
    if 'isomir_results' not in st.session_state:
        st.info("Run analysis to see visualizations")
        return
    
    results = st.session_state.isomir_results
    
    # Select miRNA to visualize
    mirnas = list(results['by_mirna'].keys())
    if not mirnas:
        st.warning("No isomiRs detected")
        return
    
    selected_mirna = st.selectbox("Select miRNA", mirnas)
    
    if selected_mirna:
        isomirs = results['by_mirna'][selected_mirna]
        
        # Bar chart of isomiR abundances
        df = pd.DataFrame(isomirs)
        df = df.sort_values('count', ascending=False)
        
        fig = px.bar(
            df,
            x='sequence',
            y='count',
            color='variant_type',
            title=f"isomiR Profile: {selected_mirna}",
            labels={'sequence': 'Sequence', 'count': 'Read Count', 'variant_type': 'Variant Type'}
        )
        fig.update_xaxes(tickangle=45)
        st.plotly_chart(fig, width="stretch")
        
        # Sequence alignment view
        st.markdown("#### Sequence Alignment")
        canonical = [i for i in isomirs if i['canonical']]
        if canonical:
            canonical_seq = canonical[0]['sequence']
            st.code(f"Canonical: {canonical_seq}")
            
            for isomir in isomirs[:10]:  # Show top 10
                if not isomir['canonical']:
                    seq = isomir['sequence']
                    # Mark differences
                    marked = ""
                    for i, (c, v) in enumerate(zip(canonical_seq.ljust(len(seq)), seq)):
                        if c == v:
                            marked += v
                        else:
                            marked += f"[{v}]"
                    st.code(f"{isomir['variant_type']:15} {marked} ({isomir['count']} reads)")
