"""
Differential Expression Module for sRNAtlas
Uses pyDESeq2 for analysis (Python-only, Streamlit Cloud compatible)
"""
import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import time
import io

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from config import config
from utils.plotting import (
    plot_pca,
    plot_volcano,
    plot_ma,
    plot_heatmap,
    plot_sample_correlation
)


def load_demo_dataset():
    """Load the demo dataset for testing"""
    demo_dir = Path(__file__).parent.parent / "data" / "demo"

    counts_file = demo_dir / "demo_counts.csv"
    metadata_file = demo_dir / "demo_metadata.csv"

    if counts_file.exists() and metadata_file.exists():
        counts = pd.read_csv(counts_file, index_col=0)
        metadata = pd.read_csv(metadata_file)
        return counts, metadata
    else:
        return None, None


def render_figure_download_buttons(fig, name: str):
    """Render download buttons for a figure"""
    st.markdown("**📥 Download Figure:**")

    col1, col2, col3, col4 = st.columns(4)

    # HTML (always available)
    html_str = fig.to_html(include_plotlyjs='cdn')
    with col1:
        st.download_button(
            "📄 HTML",
            html_str,
            f"{name}.html",
            "text/html",
            key=f"dl_html_{name}"
        )

    # Try PNG/SVG/PDF (requires kaleido)
    try:
        png_bytes = fig.to_image(format="png", width=1200, height=900, scale=2)
        with col2:
            st.download_button(
                "🖼️ PNG",
                png_bytes,
                f"{name}.png",
                "image/png",
                key=f"dl_png_{name}"
            )
    except Exception:
        with col2:
            st.caption("PNG: Install kaleido")

    try:
        svg_bytes = fig.to_image(format="svg", width=1200, height=900)
        with col3:
            st.download_button(
                "📐 SVG",
                svg_bytes,
                f"{name}.svg",
                "image/svg+xml",
                key=f"dl_svg_{name}"
            )
    except Exception:
        with col3:
            st.caption("SVG: Install kaleido")

    try:
        pdf_bytes = fig.to_image(format="pdf", width=1200, height=900)
        with col4:
            st.download_button(
                "📑 PDF",
                pdf_bytes,
                f"{name}.pdf",
                "application/pdf",
                key=f"dl_pdf_{name}"
            )
    except Exception:
        with col4:
            st.caption("PDF: Install kaleido")


def run_pydeseq2_analysis(
    count_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    design_factor: str,
    comparisons: List[Dict],
    alpha: float = 0.05,
    progress_callback=None
) -> Dict:
    """
    Run differential expression analysis using pyDESeq2

    Args:
        count_matrix: Count matrix (genes x samples)
        metadata: Sample metadata
        design_factor: Main factor for design
        comparisons: List of comparison dictionaries
        alpha: FDR threshold
        progress_callback: Optional callback function for progress updates

    Returns:
        Dictionary with analysis results
    """
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        return {
            'status': 'error',
            'error': 'pyDESeq2 not installed. Install with: pip install pydeseq2'
        }

    def update_progress(pct, msg):
        if progress_callback:
            progress_callback(pct, msg)

    try:
        update_progress(5, "Preparing data...")

        # Get sample columns (numeric only)
        sample_cols = count_matrix.select_dtypes(include=[np.number]).columns.tolist()

        # Prepare counts: pyDESeq2 expects samples x genes
        counts = count_matrix[sample_cols].T
        counts = counts.round().astype(int)

        # Ensure all counts are non-negative
        counts = counts.clip(lower=0)

        update_progress(10, "Matching samples...")

        # Prepare metadata - ensure it matches count columns
        if 'SampleID' in metadata.columns:
            metadata_indexed = metadata.set_index('SampleID')
        else:
            metadata_indexed = metadata.copy()

        # Filter to only samples present in counts
        common_samples = counts.index.intersection(metadata_indexed.index)
        if len(common_samples) == 0:
            return {
                'status': 'error',
                'error': 'No matching samples between count matrix and metadata. Check that sample names match.'
            }

        counts = counts.loc[common_samples]
        metadata_indexed = metadata_indexed.loc[common_samples]

        # Ensure design factor exists and has correct type
        if design_factor not in metadata_indexed.columns:
            return {
                'status': 'error',
                'error': f'Design factor "{design_factor}" not found in metadata'
            }

        metadata_indexed[design_factor] = metadata_indexed[design_factor].astype(str)

        results = {}
        dds = None

        total_comparisons = len(comparisons)

        for i, comp in enumerate(comparisons):
            factor = comp.get('factor', design_factor)
            ref = str(comp.get('ref'))
            test = str(comp.get('test'))
            comp_name = comp['name']

            # Calculate progress percentage for this comparison
            base_pct = 15 + int((i / total_comparisons) * 60)
            update_progress(base_pct, f"Running: {comp_name}...")

            try:
                # Create DeseqDataSet
                dds = DeseqDataSet(
                    counts=counts,
                    metadata=metadata_indexed,
                    design_factors=factor,
                    refit_cooks=True
                )

                update_progress(base_pct + 10, f"Fitting model for {comp_name}...")

                # Run DESeq2
                dds.deseq2()

                update_progress(base_pct + 20, f"Computing statistics for {comp_name}...")

                # Get statistics for this contrast
                stat_res = DeseqStats(
                    dds,
                    contrast=[factor, test, ref],
                    alpha=alpha
                )
                stat_res.summary()

                # Get results DataFrame
                res_df = stat_res.results_df.copy()

                # Rename columns to match expected format
                if 'log2FoldChange' not in res_df.columns and 'lfc' in res_df.columns:
                    res_df = res_df.rename(columns={'lfc': 'log2FoldChange'})
                if 'pvalue' not in res_df.columns and 'pval' in res_df.columns:
                    res_df = res_df.rename(columns={'pval': 'pvalue'})

                results[comp_name] = res_df

            except Exception as e:
                st.warning(f"Comparison {comp_name} failed: {str(e)}")
                continue

        if not results:
            return {
                'status': 'error',
                'error': 'All comparisons failed'
            }

        update_progress(80, "Extracting normalized counts...")

        # Get normalized counts from the last successful run
        norm_counts_df = None
        if dds is not None:
            try:
                # Get normalized counts
                norm_counts = dds.layers['normed_counts']
                norm_counts_df = pd.DataFrame(
                    norm_counts.T,
                    index=counts.columns,  # genes
                    columns=counts.index   # samples
                )
            except Exception:
                # Fallback: use simple size factor normalization
                size_factors = counts.sum(axis=1) / counts.sum(axis=1).median()
                norm_counts_df = (counts.T / size_factors).T
                norm_counts_df = pd.DataFrame(
                    norm_counts_df.T.values,
                    index=counts.columns,
                    columns=counts.index
                )

        update_progress(90, "Finalizing results...")

        return {
            'status': 'success',
            'results': results,
            'normalized_counts': norm_counts_df,
        }

    except Exception as e:
        import traceback
        return {
            'status': 'error',
            'error': f'{str(e)}\n{traceback.format_exc()}'
        }


def apply_global_fdr(results_dict: Dict[str, pd.DataFrame], progress_callback=None) -> Dict[str, pd.DataFrame]:
    """
    Apply global FDR correction across all comparisons
    """
    from scipy.stats import false_discovery_control

    if progress_callback:
        progress_callback(85, "Applying global FDR correction...")

    # Collect all p-values
    all_pvals = []
    all_info = []  # (comparison, gene)

    for comp_name, df in results_dict.items():
        pval_col = 'pvalue' if 'pvalue' in df.columns else 'pval'
        if pval_col not in df.columns:
            continue

        for gene in df.index:
            pval = df.loc[gene, pval_col]
            if not pd.isna(pval) and pval > 0:
                all_pvals.append(pval)
                all_info.append((comp_name, gene))

    # Apply BH correction
    if all_pvals:
        global_fdr = false_discovery_control(all_pvals, method='bh')

        # Map back to results
        for i, (comp_name, gene) in enumerate(all_info):
            if 'global_FDR' not in results_dict[comp_name].columns:
                results_dict[comp_name]['global_FDR'] = np.nan
            results_dict[comp_name].loc[gene, 'global_FDR'] = global_fdr[i]

    return results_dict


def render_de_page():
    """Render the DE analysis page"""
    st.header("🔬 Differential Expression Analysis")

    st.info("This tool uses **pyDESeq2** for differential expression analysis - a Python implementation of the DESeq2 algorithm.")

    tab1, tab2, tab3, tab4 = st.tabs([
        "📁 Setup",
        "▶️ Run Analysis",
        "📊 Results",
        "📈 Visualizations"
    ])

    with tab1:
        render_de_setup()

    with tab2:
        render_run_de()

    with tab3:
        render_de_results()

    with tab4:
        render_de_visualizations()


def render_de_setup():
    """Render DE setup section"""
    st.subheader("Setup")

    # Demo dataset button
    st.markdown("### 🎯 Quick Start with Demo Data")

    col1, col2 = st.columns([1, 2])
    with col1:
        if st.button("📦 Load Demo Dataset", type="secondary", width="stretch"):
            with st.spinner("Loading demo dataset..."):
                counts, metadata = load_demo_dataset()
                if counts is not None and metadata is not None:
                    st.session_state.count_matrix = counts
                    st.session_state.sample_metadata = metadata
                    st.success("✅ Demo dataset loaded!")
                    st.rerun()
                else:
                    st.error("Demo dataset files not found.")

    with col2:
        st.caption("Load a sample miRNA dataset to try the tool immediately. Contains 50 miRNAs across 6 samples (Control vs Treatment).")

    st.divider()

    # Count matrix
    st.markdown("### 📊 Count Matrix")
    counts = st.session_state.get('annotated_counts') or st.session_state.get('count_matrix')

    if counts is not None:
        st.success(f"✅ Count matrix loaded: {counts.shape[0]} features × {counts.shape[1]} columns")

        with st.expander("Preview count matrix"):
            st.dataframe(counts.head(10), width="stretch")

        if st.button("🔄 Replace Count Matrix"):
            st.session_state.count_matrix = None
            st.session_state.annotated_counts = None
            st.rerun()
    else:
        st.warning("No count matrix loaded. Please upload one below or load the demo dataset.")

        count_file = st.file_uploader("Upload count matrix (CSV/TSV)", type=['csv', 'tsv', 'txt'])
        if count_file:
            try:
                # Detect delimiter
                content = count_file.getvalue().decode('utf-8')
                sep = '\t' if '\t' in content.split('\n')[0] else ','
                count_file.seek(0)

                counts = pd.read_csv(count_file, sep=sep, index_col=0)
                st.session_state.count_matrix = counts
                st.success(f"✅ Loaded: {counts.shape[0]} features × {counts.shape[1]} columns")
                st.rerun()
            except Exception as e:
                st.error(f"Error loading file: {e}")

    # Metadata
    st.divider()
    st.markdown("### 📋 Sample Metadata")

    metadata = st.session_state.get('sample_metadata')

    if metadata is not None:
        st.success(f"✅ Metadata loaded: {len(metadata)} samples")
        st.dataframe(metadata, width="stretch", height=200)

        if st.button("🔄 Replace Metadata"):
            st.session_state.sample_metadata = None
            st.rerun()
    else:
        st.warning("No metadata loaded. Please upload one below.")

        st.markdown("""
        **Expected format:** CSV with columns:
        - `SampleID`: Sample names (must match column names in count matrix)
        - Treatment/Condition columns (e.g., `Condition`, `Treatment`, `Group`)
        """)

        metadata_file = st.file_uploader("Upload metadata CSV", type=['csv'])

        if metadata_file:
            try:
                metadata = pd.read_csv(metadata_file)
                st.session_state.sample_metadata = metadata
                st.success(f"✅ Loaded metadata: {len(metadata)} samples")
                st.dataframe(metadata, width="stretch")
                st.rerun()
            except Exception as e:
                st.error(f"Error loading metadata: {e}")


def render_run_de():
    """Render run DE analysis section"""
    st.subheader("Run Differential Expression Analysis")

    # Check data
    counts = st.session_state.get('annotated_counts') or st.session_state.get('count_matrix')
    metadata = st.session_state.get('sample_metadata')

    if counts is None:
        st.warning("⚠️ Please load count matrix in Setup tab.")
        return

    if metadata is None:
        st.warning("⚠️ Please load metadata in Setup tab.")
        return

    # Show data summary
    sample_cols = counts.select_dtypes(include=[np.number]).columns.tolist()
    st.success(f"Ready: {len(counts)} features, {len(sample_cols)} samples")

    # Experimental Design
    st.subheader("Experimental Design")

    # Get available factors from metadata
    factor_cols = [c for c in metadata.columns if c.lower() not in ['sampleid', 'sample', 'name']]

    if not factor_cols:
        st.error("No factor columns found in metadata. Please ensure metadata has treatment/condition columns.")
        return

    # Select main factor
    main_factor = st.selectbox(
        "Main Factor (grouping variable)",
        options=factor_cols,
        help="Select the main variable for comparison (e.g., Treatment, Condition)"
    )

    # Show factor levels
    levels = metadata[main_factor].unique().tolist()
    st.info(f"**{main_factor} levels:** {', '.join(map(str, levels))}")

    if len(levels) < 2:
        st.error("Main factor must have at least 2 levels for comparison.")
        return

    # Define comparisons
    st.subheader("Define Comparisons")

    comparison_mode = st.radio(
        "Comparison Mode",
        options=["vs Reference (recommended)", "All pairwise", "Custom"],
        horizontal=True
    )

    comparisons = []

    if comparison_mode == "vs Reference (recommended)":
        ref_level = st.selectbox(
            "Reference Level (control group)",
            options=levels,
            help="Select the control/baseline group"
        )
        for level in levels:
            if str(level) != str(ref_level):
                comparisons.append({
                    'name': f"{level}_vs_{ref_level}",
                    'factor': main_factor,
                    'test': str(level),
                    'ref': str(ref_level)
                })

    elif comparison_mode == "All pairwise":
        from itertools import combinations
        for l1, l2 in combinations(levels, 2):
            comparisons.append({
                'name': f"{l2}_vs_{l1}",
                'factor': main_factor,
                'test': str(l2),
                'ref': str(l1)
            })

    else:  # Custom
        n_comparisons = st.number_input("Number of comparisons", 1, 10, 1)

        for i in range(int(n_comparisons)):
            col1, col2 = st.columns(2)
            with col1:
                test = st.selectbox(f"Test group {i+1}", options=levels, key=f"test_{i}")
            with col2:
                ref = st.selectbox(f"Reference group {i+1}", options=levels, key=f"ref_{i}")

            if str(test) != str(ref):
                comparisons.append({
                    'name': f"{test}_vs_{ref}",
                    'factor': main_factor,
                    'test': str(test),
                    'ref': str(ref)
                })

    # Show planned comparisons
    if comparisons:
        st.markdown(f"**Comparisons to perform ({len(comparisons)}):**")
        for comp in comparisons:
            st.text(f"  • {comp['name']}")
    else:
        st.warning("No valid comparisons defined.")
        return

    # Analysis settings
    st.divider()
    st.subheader("Analysis Settings")

    col1, col2 = st.columns(2)

    with col1:
        alpha = st.slider("FDR Threshold", 0.01, 0.20, 0.05, 0.01)

    with col2:
        apply_global = st.checkbox("Apply global FDR correction", value=True,
                                   help="Correct for multiple testing across all comparisons")

    # Run analysis
    st.divider()

    if st.button("🚀 Run Differential Expression Analysis", type="primary", width="stretch"):
        # Progress tracking
        progress_bar = st.progress(0)
        status_text = st.empty()

        def progress_callback(pct, msg):
            progress_bar.progress(min(pct, 100))
            status_text.text(msg)

        progress_callback(5, "Initializing analysis...")

        result = run_pydeseq2_analysis(
            counts, metadata, main_factor, comparisons, alpha,
            progress_callback=progress_callback
        )

        if result['status'] == 'success':
            # Apply global FDR if requested
            if apply_global and result['results']:
                progress_callback(85, "Applying global FDR correction...")
                result['results'] = apply_global_fdr(result['results'], progress_callback)

            progress_callback(100, "Complete!")
            time.sleep(0.5)

            st.session_state.de_results = result

            st.success("✅ Analysis complete!")

            # Summary
            st.subheader("Summary")
            cols = st.columns(len(result['results']))
            for i, (comp_name, res_df) in enumerate(result['results'].items()):
                padj_col = 'padj' if 'padj' in res_df.columns else 'global_FDR'
                if padj_col in res_df.columns:
                    n_sig = (res_df[padj_col] < alpha).sum()
                else:
                    n_sig = 0
                with cols[i]:
                    st.metric(comp_name, f"{n_sig} significant")

            st.balloons()

        else:
            progress_bar.progress(100)
            st.error(f"❌ Analysis failed: {result['error']}")


def render_de_results():
    """Render DE results section"""
    st.subheader("Differential Expression Results")

    if st.session_state.get('de_results') is None:
        st.info("No results yet. Run analysis in the 'Run Analysis' tab.")
        return

    results = st.session_state.de_results['results']

    # Select comparison
    comparison = st.selectbox(
        "Select Comparison",
        options=list(results.keys())
    )

    res_df = results[comparison]

    # Determine column names
    padj_col = 'padj' if 'padj' in res_df.columns else ('global_FDR' if 'global_FDR' in res_df.columns else None)
    lfc_col = 'log2FoldChange' if 'log2FoldChange' in res_df.columns else 'lfc'
    pval_col = 'pvalue' if 'pvalue' in res_df.columns else 'pval'

    # Summary metrics
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric("Total Genes", len(res_df))

    with col2:
        if padj_col:
            n_sig = (res_df[padj_col] < 0.05).sum()
        else:
            n_sig = 0
        st.metric("Significant (FDR<0.05)", n_sig)

    with col3:
        if padj_col and lfc_col in res_df.columns:
            n_up = ((res_df[padj_col] < 0.05) & (res_df[lfc_col] > 0)).sum()
        else:
            n_up = 0
        st.metric("Up-regulated", n_up)

    with col4:
        if padj_col and lfc_col in res_df.columns:
            n_down = ((res_df[padj_col] < 0.05) & (res_df[lfc_col] < 0)).sum()
        else:
            n_down = 0
        st.metric("Down-regulated", n_down)

    # Results table
    st.subheader("Results Table")

    col1, col2 = st.columns(2)
    with col1:
        available_sort_cols = [c for c in [padj_col, lfc_col, 'baseMean'] if c and c in res_df.columns]
        if available_sort_cols:
            sort_by = st.selectbox("Sort by", options=available_sort_cols)
        else:
            sort_by = res_df.columns[0]
    with col2:
        ascending = st.checkbox("Ascending", value=True)

    show_significant = st.checkbox("Show only significant genes (FDR < 0.05)", value=False)

    display_df = res_df.copy()
    if show_significant and padj_col:
        display_df = display_df[display_df[padj_col] < 0.05]

    if sort_by in display_df.columns:
        display_df = display_df.sort_values(sort_by, ascending=ascending)

    st.dataframe(display_df, width="stretch", height=400)

    # Download
    csv = res_df.to_csv()
    st.download_button(
        f"📥 Download {comparison} Results",
        csv,
        f"{comparison}_results.csv",
        "text/csv"
    )

    # Download all results
    st.divider()
    if st.button("📥 Download All Results (ZIP)"):
        import zipfile

        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, 'w', zipfile.ZIP_DEFLATED) as zf:
            for comp_name, df in results.items():
                zf.writestr(f"{comp_name}_results.csv", df.to_csv())

        st.download_button(
            "📥 Download ZIP",
            buffer.getvalue(),
            "de_results.zip",
            "application/zip"
        )


def render_de_visualizations():
    """Render DE visualizations section"""
    st.subheader("Visualizations")

    if st.session_state.get('de_results') is None:
        st.info("No results yet. Run analysis first.")
        return

    de_results = st.session_state.de_results
    results = de_results['results']

    # Tip about publication-ready exports
    st.info("💡 **Tip:** Each visualization includes download buttons for publication-ready figures (PNG, SVG, PDF). For high-resolution exports, install `kaleido`: `pip install kaleido`")

    viz_type = st.selectbox(
        "Visualization Type",
        options=["Volcano Plot", "MA Plot", "PCA", "Sample Correlation", "Heatmap"]
    )

    if viz_type == "Volcano Plot":
        comparison = st.selectbox("Select Comparison", options=list(results.keys()), key="volcano_comp")
        res_df = results[comparison]

        lfc_col = 'log2FoldChange' if 'log2FoldChange' in res_df.columns else 'lfc'
        pval_col = 'pvalue' if 'pvalue' in res_df.columns else 'pval'
        padj_col = 'padj' if 'padj' in res_df.columns else 'global_FDR'

        col1, col2 = st.columns(2)
        with col1:
            lfc_thresh = st.slider("LFC Threshold", 0.0, 3.0, 0.585, 0.1, key="volcano_lfc")
        with col2:
            fdr_thresh = st.slider("FDR Threshold", 0.01, 0.20, 0.05, 0.01, key="volcano_fdr")

        fig = plot_volcano(
            res_df,
            lfc_col=lfc_col,
            pval_col=pval_col,
            padj_col=padj_col,
            lfc_threshold=lfc_thresh,
            padj_threshold=fdr_thresh,
            title=f"Volcano Plot: {comparison}"
        )
        st.plotly_chart(fig, width="stretch")

        # Download buttons
        render_figure_download_buttons(fig, f"volcano_{comparison}")

    elif viz_type == "MA Plot":
        comparison = st.selectbox("Select Comparison", options=list(results.keys()), key="ma_comp")
        res_df = results[comparison]

        basemean_col = 'baseMean' if 'baseMean' in res_df.columns else (res_df.columns[0] if len(res_df.columns) > 0 else None)
        lfc_col = 'log2FoldChange' if 'log2FoldChange' in res_df.columns else 'lfc'
        padj_col = 'padj' if 'padj' in res_df.columns else 'global_FDR'

        if basemean_col and lfc_col in res_df.columns and padj_col in res_df.columns:
            fig = plot_ma(
                res_df,
                basemean_col=basemean_col,
                lfc_col=lfc_col,
                padj_col=padj_col,
                title=f"MA Plot: {comparison}"
            )
            st.plotly_chart(fig, width="stretch")

            render_figure_download_buttons(fig, f"ma_{comparison}")
        else:
            st.warning("Required columns not available for MA plot")

    elif viz_type == "PCA":
        norm_counts = de_results.get('normalized_counts')
        metadata = st.session_state.get('sample_metadata')

        if norm_counts is not None and metadata is not None:
            color_options = [c for c in metadata.columns if c.lower() not in ['sampleid', 'sample']]
            if color_options:
                color_by = st.selectbox("Color by", options=color_options, key="pca_color")
                fig = plot_pca(norm_counts, metadata, color_by=color_by)
                st.plotly_chart(fig, width="stretch")

                render_figure_download_buttons(fig, f"pca_{color_by}")
            else:
                st.warning("No grouping columns available in metadata")
        else:
            st.warning("Normalized counts or metadata not available for PCA")

    elif viz_type == "Sample Correlation":
        norm_counts = de_results.get('normalized_counts')

        if norm_counts is not None:
            method = st.selectbox("Correlation Method", options=['spearman', 'pearson'], key="corr_method")
            fig = plot_sample_correlation(norm_counts, method=method)
            st.plotly_chart(fig, width="stretch")

            render_figure_download_buttons(fig, f"correlation_{method}")
        else:
            st.warning("Normalized counts not available")

    elif viz_type == "Heatmap":
        comparison = st.selectbox("Select Comparison", options=list(results.keys()), key="heatmap_comp")
        res_df = results[comparison]
        norm_counts = de_results.get('normalized_counts')

        if norm_counts is not None:
            padj_col = 'padj' if 'padj' in res_df.columns else 'global_FDR'

            if padj_col in res_df.columns:
                sig_genes = res_df[res_df[padj_col] < 0.05].index.tolist()

                if sig_genes:
                    max_genes = min(100, len(sig_genes))
                    n_genes = st.slider("Number of genes", 10, max_genes, min(50, max_genes), key="heatmap_n")

                    top_genes = res_df.loc[sig_genes].nsmallest(n_genes, padj_col).index
                    heatmap_data = norm_counts.loc[norm_counts.index.isin(top_genes)]

                    if len(heatmap_data) > 0:
                        # Z-score normalize
                        heatmap_z = (heatmap_data.T - heatmap_data.mean(axis=1)) / (heatmap_data.std(axis=1) + 1e-10)
                        heatmap_z = heatmap_z.T

                        fig = plot_heatmap(
                            heatmap_z,
                            title=f"Top {n_genes} DE Genes: {comparison}"
                        )
                        st.plotly_chart(fig, width="stretch")

                        render_figure_download_buttons(fig, f"heatmap_{comparison}")
                    else:
                        st.warning("No matching genes found in normalized counts")
                else:
                    st.warning("No significant genes found")
            else:
                st.warning("FDR column not available")
        else:
            st.warning("Normalized counts not available")
