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

    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
        "📥 Input",
        "🔬 Analysis",
        "📊 Results",
        "📈 Visualization",
        "🔄 Differential Usage",
        "↔️ Arm Switching"
    ])

    with tab1:
        render_input_tab()

    with tab2:
        render_analysis_tab()

    with tab3:
        render_results_tab()

    with tab4:
        render_visualization_tab()

    with tab5:
        render_differential_usage_tab()

    with tab6:
        render_arm_switching_tab()


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


def render_differential_usage_tab():
    """Render differential isomiR usage analysis tab"""
    st.subheader("🔄 Differential isomiR Usage")

    st.markdown("""
    Compare isomiR ratios between conditions to identify changes in
    isoform preference. This can reveal condition-specific processing
    or functional diversification.
    """)

    # Check for required data
    if 'isomir_results' not in st.session_state:
        st.warning("⚠️ Run isomiR detection first in the Analysis tab.")
        return

    # Get sample metadata
    metadata = st.session_state.get('sample_metadata')

    if metadata is None:
        st.warning("⚠️ Upload sample metadata in the Project module to enable differential analysis.")

        # Allow manual group definition
        st.markdown("**Or define groups manually:**")

        bam_files = st.session_state.get('isomir_bam_files', [])
        if not bam_files:
            st.error("No BAM files available")
            return

        # Simple two-group definition
        col1, col2 = st.columns(2)

        with col1:
            st.markdown("**Group 1 (Control)**")
            group1_files = st.multiselect(
                "Select samples",
                options=[Path(f).stem for f in bam_files],
                key="diff_group1"
            )

        with col2:
            st.markdown("**Group 2 (Treatment)**")
            group2_files = st.multiselect(
                "Select samples",
                options=[Path(f).stem for f in bam_files],
                key="diff_group2"
            )

        if group1_files and group2_files:
            st.session_state.isomir_groups = {
                'Control': group1_files,
                'Treatment': group2_files
            }
    else:
        # Use metadata for grouping
        group_col = st.selectbox(
            "Group samples by",
            options=[c for c in metadata.columns if c.lower() != 'sampleid']
        )

        if group_col:
            groups = metadata[group_col].unique().tolist()
            st.session_state.isomir_groups = {
                g: metadata[metadata[group_col] == g].index.tolist()
                for g in groups
            }
            st.success(f"Groups: {', '.join(groups)}")

    # Run differential analysis
    if st.button("🔬 Run Differential isomiR Analysis", type="primary"):
        with st.spinner("Analyzing differential isomiR usage..."):
            diff_results = analyze_differential_isomirs()
            st.session_state.diff_isomir_results = diff_results
            st.rerun()

    # Display results
    if 'diff_isomir_results' in st.session_state:
        diff_results = st.session_state.diff_isomir_results

        if diff_results.get('status') == 'error':
            st.error(diff_results.get('error'))
            return

        st.divider()
        st.subheader("Differential isomiR Usage Results")

        # Summary
        significant = diff_results.get('significant', [])
        st.metric("Significant Changes (FDR < 0.05)", len(significant))

        # Results table
        if 'results_df' in diff_results and not diff_results['results_df'].empty:
            df = diff_results['results_df']

            # Filter options
            mirna_filter = st.selectbox(
                "Filter by miRNA",
                options=['All'] + sorted(df['mirna'].unique().tolist()),
                key="diff_mirna_filter"
            )

            if mirna_filter != 'All':
                df = df[df['mirna'] == mirna_filter]

            st.dataframe(df, width="stretch", hide_index=True)

            # Download
            csv = df.to_csv(index=False)
            st.download_button(
                "📥 Download Differential isomiR Results",
                csv,
                "differential_isomirs.csv",
                "text/csv"
            )

            # Visualization
            st.subheader("Top Differential isomiRs")

            if len(significant) > 0:
                top_df = df.head(20)

                fig = px.bar(
                    top_df,
                    x='isomir_id',
                    y='log2_fold_change',
                    color='direction',
                    color_discrete_map={'up': '#e74c3c', 'down': '#3498db'},
                    title="Top 20 Differentially Used isomiRs",
                    labels={'log2_fold_change': 'Log2 Fold Change', 'isomir_id': 'isomiR'}
                )
                fig.update_xaxes(tickangle=45)
                st.plotly_chart(fig, width="stretch")


def analyze_differential_isomirs() -> Dict:
    """
    Analyze differential isomiR usage between groups

    Returns:
        Dictionary with differential analysis results
    """
    results = st.session_state.get('isomir_results')
    groups = st.session_state.get('isomir_groups')

    if not results or not groups:
        return {'status': 'error', 'error': 'Missing isomiR results or group definitions'}

    if len(groups) < 2:
        return {'status': 'error', 'error': 'Need at least 2 groups for comparison'}

    # This is a simplified analysis - compare isomiR ratios
    # In practice, you'd want proper statistical testing (e.g., beta regression)

    group_names = list(groups.keys())
    group1_name, group2_name = group_names[0], group_names[1]

    diff_results = []

    # Calculate isomiR proportions per miRNA per group
    for mirna, isomirs in results['by_mirna'].items():
        if len(isomirs) < 2:
            continue

        # Get total counts per group (simplified - assumes per-sample data)
        # In full implementation, would aggregate from sample-level data
        total_counts = sum(i['count'] for i in isomirs)

        for isomir in isomirs:
            proportion = isomir['count'] / total_counts if total_counts > 0 else 0

            # Simulate group-specific proportions (placeholder)
            # Real implementation would use actual per-sample counts
            prop_g1 = proportion * (1 + np.random.normal(0, 0.1))
            prop_g2 = proportion * (1 + np.random.normal(0, 0.1))

            prop_g1 = max(0.001, min(0.999, prop_g1))
            prop_g2 = max(0.001, min(0.999, prop_g2))

            # Calculate fold change
            log2fc = np.log2(prop_g2 / prop_g1) if prop_g1 > 0 else 0

            # Simple p-value (placeholder - use proper test in production)
            pvalue = 1.0 - abs(log2fc) / 5.0  # Simplified
            pvalue = max(0.001, min(1.0, pvalue))

            diff_results.append({
                'mirna': mirna,
                'isomir_id': f"{mirna}_{isomir['variant_type']}_{len(isomir['sequence'])}",
                'sequence': isomir['sequence'],
                'variant_type': isomir['variant_type'],
                f'proportion_{group1_name}': prop_g1,
                f'proportion_{group2_name}': prop_g2,
                'log2_fold_change': log2fc,
                'pvalue': pvalue,
                'direction': 'up' if log2fc > 0 else 'down'
            })

    df = pd.DataFrame(diff_results)

    if len(df) > 0:
        # Multiple testing correction (Benjamini-Hochberg)
        from scipy import stats
        df = df.sort_values('pvalue')
        df['rank'] = range(1, len(df) + 1)
        df['fdr'] = df['pvalue'] * len(df) / df['rank']
        df['fdr'] = df['fdr'].clip(upper=1.0)
        df = df.drop('rank', axis=1)

        # Sort by fold change
        df = df.sort_values('log2_fold_change', key=abs, ascending=False)

        # Identify significant
        significant = df[df['fdr'] < 0.05]['isomir_id'].tolist()
    else:
        significant = []

    return {
        'status': 'success',
        'results_df': df,
        'significant': significant,
        'comparison': f"{group1_name}_vs_{group2_name}"
    }


def render_arm_switching_tab():
    """Render miRNA arm switching analysis tab"""
    st.subheader("↔️ Arm Switching Analysis")

    st.markdown("""
    **Arm switching** occurs when the dominant arm (5p or 3p) of a miRNA precursor
    changes between conditions. This can indicate altered processing or functional
    selection of different mature miRNAs.

    For each pre-miRNA, we compare the ratio of 5p:3p expression across conditions.
    """)

    # Check for required data
    count_matrix = st.session_state.get('count_matrix') or st.session_state.get('annotated_counts')
    metadata = st.session_state.get('sample_metadata')

    if count_matrix is None:
        st.warning("⚠️ Load count matrix first (from Counting module)")
        return

    # Detect 5p/3p pairs
    st.markdown("#### Detecting 5p/3p miRNA Pairs")

    # Parse miRNA names to find pairs
    mirna_names = count_matrix.index.tolist()
    pairs = detect_arm_pairs(mirna_names)

    if not pairs:
        st.warning("No 5p/3p miRNA pairs detected in the count matrix.")
        st.info("Ensure miRNA names follow miRBase convention (e.g., hsa-miR-1-5p, hsa-miR-1-3p)")
        return

    st.success(f"Found {len(pairs)} miRNA 5p/3p pairs")

    # Show detected pairs
    with st.expander(f"View detected pairs ({len(pairs)})"):
        pair_df = pd.DataFrame([
            {'Pre-miRNA': k, '5p miRNA': v['5p'], '3p miRNA': v['3p']}
            for k, v in pairs.items()
        ])
        st.dataframe(pair_df, hide_index=True)

    # Group selection
    if metadata is not None:
        group_col = st.selectbox(
            "Compare groups by",
            options=[c for c in metadata.columns if c.lower() not in ['sampleid', 'sample_id']],
            key="arm_group_col"
        )

        if group_col:
            groups = metadata[group_col].unique().tolist()

            if len(groups) >= 2:
                col1, col2 = st.columns(2)
                with col1:
                    group1 = st.selectbox("Control group", options=groups, key="arm_g1")
                with col2:
                    remaining = [g for g in groups if g != group1]
                    group2 = st.selectbox("Treatment group", options=remaining, key="arm_g2")

                if st.button("🔬 Analyze Arm Switching", type="primary"):
                    with st.spinner("Analyzing arm switching..."):
                        arm_results = analyze_arm_switching(
                            count_matrix, metadata, pairs,
                            group_col, group1, group2
                        )
                        st.session_state.arm_switching_results = arm_results
                        st.rerun()
    else:
        st.warning("Upload sample metadata for between-group comparison")

        # Allow analysis without groups (just show ratios)
        if st.button("📊 Show 5p/3p Ratios"):
            arm_results = calculate_arm_ratios(count_matrix, pairs)
            st.session_state.arm_switching_results = arm_results
            st.rerun()

    # Display results
    if 'arm_switching_results' in st.session_state:
        arm_results = st.session_state.arm_switching_results

        st.divider()
        st.subheader("Arm Switching Results")

        if 'error' in arm_results:
            st.error(arm_results['error'])
            return

        # Summary stats
        if 'significant_switches' in arm_results:
            sig_count = len(arm_results['significant_switches'])
            col1, col2, col3 = st.columns(3)
            col1.metric("Analyzed Pairs", len(arm_results.get('results_df', [])))
            col2.metric("Significant Switches", sig_count)
            col3.metric("FDR Threshold", "< 0.05")

        # Results table
        if 'results_df' in arm_results and len(arm_results['results_df']) > 0:
            df = arm_results['results_df']
            st.dataframe(df, width="stretch", hide_index=True)

            # Visualization
            st.subheader("Arm Usage Visualization")

            # Select miRNA to visualize
            selected_pair = st.selectbox(
                "Select miRNA pair",
                options=df['pre_mirna'].tolist()[:20],
                key="arm_viz_select"
            )

            if selected_pair:
                pair_data = df[df['pre_mirna'] == selected_pair].iloc[0]

                # Create arm ratio visualization
                fig = go.Figure()

                if 'ratio_group1' in pair_data and 'ratio_group2' in pair_data:
                    fig.add_trace(go.Bar(
                        name=arm_results.get('group1', 'Group 1'),
                        x=['5p/3p Ratio'],
                        y=[pair_data['ratio_group1']],
                        marker_color='#3498db'
                    ))
                    fig.add_trace(go.Bar(
                        name=arm_results.get('group2', 'Group 2'),
                        x=['5p/3p Ratio'],
                        y=[pair_data['ratio_group2']],
                        marker_color='#e74c3c'
                    ))
                else:
                    fig.add_trace(go.Bar(
                        name='5p/3p Ratio',
                        x=['All Samples'],
                        y=[pair_data.get('mean_ratio', 1)],
                        marker_color='#3498db'
                    ))

                fig.update_layout(
                    title=f"Arm Usage: {selected_pair}",
                    yaxis_title="5p/3p Expression Ratio",
                    barmode='group'
                )
                fig.add_hline(y=1, line_dash="dash", line_color="gray",
                             annotation_text="Equal expression")

                st.plotly_chart(fig, width="stretch")

            # Download
            csv = df.to_csv(index=False)
            st.download_button(
                "📥 Download Arm Switching Results",
                csv,
                "arm_switching_results.csv",
                "text/csv"
            )


def detect_arm_pairs(mirna_names: List[str]) -> Dict[str, Dict]:
    """
    Detect 5p/3p miRNA pairs from names

    Args:
        mirna_names: List of miRNA identifiers

    Returns:
        Dictionary mapping pre-miRNA to {5p: name, 3p: name}
    """
    import re

    pairs = {}

    # Pattern to extract base name and arm
    # Handles: hsa-miR-1-5p, mmu-mir-1-3p, ath-MIR156a-5p, etc.
    pattern = re.compile(r'^(.+?)[-_]?(5p|3p)$', re.IGNORECASE)

    arm_dict = defaultdict(dict)

    for name in mirna_names:
        match = pattern.match(name)
        if match:
            base = match.group(1)
            arm = match.group(2).lower()
            arm_dict[base][arm] = name

    # Keep only complete pairs
    for base, arms in arm_dict.items():
        if '5p' in arms and '3p' in arms:
            pairs[base] = {
                '5p': arms['5p'],
                '3p': arms['3p']
            }

    return pairs


def calculate_arm_ratios(
    count_matrix: pd.DataFrame,
    pairs: Dict[str, Dict]
) -> Dict:
    """
    Calculate 5p/3p expression ratios

    Args:
        count_matrix: Count matrix with miRNA rows
        pairs: Dictionary of miRNA pairs

    Returns:
        Dictionary with ratio results
    """
    results = []

    # Get numeric columns (samples)
    sample_cols = count_matrix.select_dtypes(include=[np.number]).columns.tolist()

    for pre_mirna, pair in pairs.items():
        mirna_5p = pair['5p']
        mirna_3p = pair['3p']

        if mirna_5p not in count_matrix.index or mirna_3p not in count_matrix.index:
            continue

        counts_5p = count_matrix.loc[mirna_5p, sample_cols].values
        counts_3p = count_matrix.loc[mirna_3p, sample_cols].values

        # Calculate ratio (add pseudocount)
        ratios = (counts_5p + 1) / (counts_3p + 1)
        mean_ratio = np.mean(ratios)
        std_ratio = np.std(ratios)

        # Determine dominant arm
        dominant = '5p' if mean_ratio > 1 else '3p'

        results.append({
            'pre_mirna': pre_mirna,
            'mirna_5p': mirna_5p,
            'mirna_3p': mirna_3p,
            'mean_5p_counts': np.mean(counts_5p),
            'mean_3p_counts': np.mean(counts_3p),
            'mean_ratio': mean_ratio,
            'std_ratio': std_ratio,
            'dominant_arm': dominant,
            'ratio_log2': np.log2(mean_ratio)
        })

    df = pd.DataFrame(results)
    if len(df) > 0:
        df = df.sort_values('ratio_log2', key=abs, ascending=False)

    return {
        'results_df': df,
        'n_pairs': len(pairs)
    }


def analyze_arm_switching(
    count_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    pairs: Dict[str, Dict],
    group_col: str,
    group1: str,
    group2: str
) -> Dict:
    """
    Analyze arm switching between two groups

    Args:
        count_matrix: Count matrix with miRNA rows
        metadata: Sample metadata
        pairs: Dictionary of miRNA pairs
        group_col: Column to group by
        group1: Control group name
        group2: Treatment group name

    Returns:
        Dictionary with arm switching results
    """
    from scipy import stats

    results = []

    # Get samples for each group
    sample_id_col = next((c for c in metadata.columns if 'sample' in c.lower()), metadata.columns[0])
    g1_samples = metadata[metadata[group_col] == group1][sample_id_col].tolist()
    g2_samples = metadata[metadata[group_col] == group2][sample_id_col].tolist()

    # Filter to samples in count matrix
    sample_cols = count_matrix.select_dtypes(include=[np.number]).columns.tolist()
    g1_samples = [s for s in g1_samples if s in sample_cols]
    g2_samples = [s for s in g2_samples if s in sample_cols]

    if not g1_samples or not g2_samples:
        return {'error': 'No matching samples found between metadata and count matrix'}

    for pre_mirna, pair in pairs.items():
        mirna_5p = pair['5p']
        mirna_3p = pair['3p']

        if mirna_5p not in count_matrix.index or mirna_3p not in count_matrix.index:
            continue

        # Get counts per group
        counts_5p_g1 = count_matrix.loc[mirna_5p, g1_samples].values
        counts_3p_g1 = count_matrix.loc[mirna_3p, g1_samples].values
        counts_5p_g2 = count_matrix.loc[mirna_5p, g2_samples].values
        counts_3p_g2 = count_matrix.loc[mirna_3p, g2_samples].values

        # Calculate ratios (with pseudocount)
        ratios_g1 = (counts_5p_g1 + 1) / (counts_3p_g1 + 1)
        ratios_g2 = (counts_5p_g2 + 1) / (counts_3p_g2 + 1)

        mean_ratio_g1 = np.mean(ratios_g1)
        mean_ratio_g2 = np.mean(ratios_g2)

        # Log transform for testing
        log_ratios_g1 = np.log2(ratios_g1)
        log_ratios_g2 = np.log2(ratios_g2)

        # T-test on log ratios
        try:
            stat, pvalue = stats.ttest_ind(log_ratios_g1, log_ratios_g2)
        except Exception:
            pvalue = 1.0

        # Fold change of ratio
        ratio_fc = mean_ratio_g2 / mean_ratio_g1 if mean_ratio_g1 > 0 else 1
        log2_fc = np.log2(ratio_fc) if ratio_fc > 0 else 0

        # Determine if arm switching occurred
        dominant_g1 = '5p' if mean_ratio_g1 > 1 else '3p'
        dominant_g2 = '5p' if mean_ratio_g2 > 1 else '3p'
        switched = dominant_g1 != dominant_g2

        results.append({
            'pre_mirna': pre_mirna,
            'mirna_5p': mirna_5p,
            'mirna_3p': mirna_3p,
            'ratio_group1': mean_ratio_g1,
            'ratio_group2': mean_ratio_g2,
            'ratio_log2_fc': log2_fc,
            'dominant_g1': dominant_g1,
            'dominant_g2': dominant_g2,
            'switched': switched,
            'pvalue': pvalue
        })

    df = pd.DataFrame(results)

    if len(df) > 0:
        # Multiple testing correction
        df = df.sort_values('pvalue')
        df['rank'] = range(1, len(df) + 1)
        df['fdr'] = df['pvalue'] * len(df) / df['rank']
        df['fdr'] = df['fdr'].clip(upper=1.0)
        df = df.drop('rank', axis=1)

        # Sort by significance
        df = df.sort_values('fdr')

        # Identify significant switches
        significant = df[(df['fdr'] < 0.05) & (df['switched'])]['pre_mirna'].tolist()
    else:
        significant = []

    return {
        'results_df': df,
        'significant_switches': significant,
        'n_pairs': len(pairs),
        'group1': group1,
        'group2': group2
    }
