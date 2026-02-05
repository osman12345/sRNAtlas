"""
GO/Pathway Enrichment Module for sRNAtlas
Functional enrichment analysis using g:Profiler (Python, no R required)
"""
import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional, Tuple

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from config import config
from utils.plotting import plot_enrichment_dotplot


# Organism mapping for KEGG
KEGG_ORGANISMS = {
    'Arabidopsis thaliana': 'ath',
    'Medicago truncatula': 'mtr',
    'Oryza sativa': 'osa',
    'Zea mays': 'zma',
    'Glycine max': 'gmx',
    'Homo sapiens': 'hsa',
    'Mus musculus': 'mmu',
    'Drosophila melanogaster': 'dme',
    'Caenorhabditis elegans': 'cel',
}


def run_gprofiler_enrichment(
    genes: List[str],
    organism: str = 'athaliana',
    sources: List[str] = ['GO:BP', 'GO:MF', 'GO:CC', 'KEGG'],
    background: Optional[List[str]] = None
) -> Dict:
    """
    Run enrichment analysis using g:Profiler API (no R required)

    Args:
        genes: List of gene IDs
        organism: Organism code (e.g., 'athaliana', 'hsapiens')
        sources: Data sources to query
        background: Background gene list

    Returns:
        Dictionary with enrichment results
    """
    try:
        from gprofiler import GProfiler
    except ImportError:
        return {
            'status': 'error',
            'error': 'gprofiler-official not installed. Install with: pip install gprofiler-official'
        }

    try:
        gp = GProfiler(return_dataframe=True)

        results = gp.profile(
            organism=organism,
            query=genes,
            background=background,
            sources=sources,
            user_threshold=0.05,
            all_results=False,
            ordered=False,
            no_iea=False,
            combined=False
        )

        if results.empty:
            return {
                'status': 'success',
                'results': {},
                'message': 'No significant enrichment found'
            }

        # Organize results by source
        organized = {}
        for source in results['source'].unique():
            source_df = results[results['source'] == source].copy()
            organized[source] = source_df

        return {
            'status': 'success',
            'results': organized
        }

    except Exception as e:
        return {
            'status': 'error',
            'error': str(e)
        }


def extract_target_genes(
    de_results: pd.DataFrame,
    mirna_targets: Optional[pd.DataFrame] = None,
    padj_threshold: float = 0.05,
    lfc_threshold: float = 0.585,
    direction: str = 'both'
) -> List[str]:
    """
    Extract target genes for enrichment analysis

    For miRNAs, this would involve finding their target genes.
    For other ncRNAs, this extracts the significant ones directly.

    Args:
        de_results: DE results DataFrame
        mirna_targets: Optional miRNA-target mapping
        padj_threshold: FDR threshold
        lfc_threshold: Log2FC threshold
        direction: 'up', 'down', or 'both'

    Returns:
        List of gene IDs for enrichment
    """
    padj_col = 'padj' if 'padj' in de_results.columns else 'global_FDR'
    lfc_col = 'log2FoldChange' if 'log2FoldChange' in de_results.columns else 'lfc'

    # Filter by significance
    sig_mask = de_results[padj_col] < padj_threshold

    if direction == 'up':
        sig_mask = sig_mask & (de_results[lfc_col] > lfc_threshold)
    elif direction == 'down':
        sig_mask = sig_mask & (de_results[lfc_col] < -lfc_threshold)
    else:  # both
        sig_mask = sig_mask & (de_results[lfc_col].abs() > lfc_threshold)

    sig_genes = de_results[sig_mask].index.tolist()

    # If miRNA targets provided, map miRNAs to their targets
    if mirna_targets is not None and len(sig_genes) > 0:
        target_genes = []
        for mirna in sig_genes:
            if mirna in mirna_targets.index:
                targets = mirna_targets.loc[mirna, 'targets']
                if isinstance(targets, str):
                    target_genes.extend(targets.split(','))
                elif isinstance(targets, list):
                    target_genes.extend(targets)
        return list(set(target_genes))

    return sig_genes


def render_enrichment_page():
    """Render the enrichment analysis page"""
    st.header("🧬 GO & Pathway Enrichment Analysis")

    # Check requirements
    with st.expander("📦 Check Requirements"):
        try:
            import gprofiler
            st.success("✅ g:Profiler (Python) - Ready for enrichment analysis")
        except ImportError:
            st.error("❌ g:Profiler not installed. Install with: pip install gprofiler-official")

    tab1, tab2, tab3 = st.tabs([
        "📁 Input",
        "▶️ Run Analysis",
        "📊 Results"
    ])

    with tab1:
        render_enrichment_input()

    with tab2:
        render_run_enrichment()

    with tab3:
        render_enrichment_results()


def render_enrichment_input():
    """Render enrichment input section"""
    st.subheader("Input Gene List")

    input_source = st.radio(
        "Gene List Source",
        options=["From DE Results", "Upload Gene List", "Manual Entry"]
    )

    if input_source == "From DE Results":
        if st.session_state.get('de_results') is None:
            st.warning("No DE results available. Run DE analysis first.")
            return

        results = st.session_state.de_results['results']

        comparison = st.selectbox(
            "Select Comparison",
            options=list(results.keys())
        )

        res_df = results[comparison]

        col1, col2 = st.columns(2)

        with col1:
            padj_threshold = st.slider("FDR Threshold", 0.01, 0.20, 0.05, 0.01)
            lfc_threshold = st.slider("LFC Threshold", 0.0, 2.0, 0.585, 0.1)

        with col2:
            direction = st.selectbox(
                "Direction",
                options=["Both", "Up-regulated", "Down-regulated"]
            )

        # Extract genes
        direction_map = {"Both": "both", "Up-regulated": "up", "Down-regulated": "down"}
        genes = extract_target_genes(
            res_df,
            padj_threshold=padj_threshold,
            lfc_threshold=lfc_threshold,
            direction=direction_map[direction]
        )

        st.info(f"Selected {len(genes)} genes for enrichment analysis")

        if genes:
            with st.expander("View Gene List"):
                st.text('\n'.join(genes[:100]))
                if len(genes) > 100:
                    st.text(f"... and {len(genes) - 100} more")

            st.session_state.enrichment_genes = genes

    elif input_source == "Upload Gene List":
        gene_file = st.file_uploader(
            "Upload gene list (one gene per line)",
            type=['txt', 'csv']
        )

        if gene_file:
            content = gene_file.getvalue().decode('utf-8')
            genes = [g.strip() for g in content.split('\n') if g.strip()]
            st.info(f"Loaded {len(genes)} genes")
            st.session_state.enrichment_genes = genes

    else:  # Manual Entry
        genes_text = st.text_area(
            "Enter genes (one per line)",
            height=200
        )

        if genes_text:
            genes = [g.strip() for g in genes_text.split('\n') if g.strip()]
            st.info(f"Entered {len(genes)} genes")
            st.session_state.enrichment_genes = genes

    # Background genes (optional)
    st.divider()
    st.subheader("Background Genes (Optional)")

    use_background = st.checkbox("Use custom background")

    if use_background:
        bg_source = st.radio(
            "Background Source",
            options=["All detected genes", "Upload list"]
        )

        if bg_source == "All detected genes":
            counts = st.session_state.get('count_matrix')
            if counts is not None:
                background = counts.index.tolist()
                st.info(f"Using {len(background)} genes as background")
                st.session_state.enrichment_background = background
            else:
                st.warning("No count matrix available")
        else:
            bg_file = st.file_uploader(
                "Upload background gene list",
                type=['txt', 'csv'],
                key='bg_upload'
            )

            if bg_file:
                content = bg_file.getvalue().decode('utf-8')
                background = [g.strip() for g in content.split('\n') if g.strip()]
                st.info(f"Loaded {len(background)} background genes")
                st.session_state.enrichment_background = background


def render_run_enrichment():
    """Render run enrichment section"""
    st.subheader("Run Enrichment Analysis")

    genes = st.session_state.get('enrichment_genes')

    if not genes:
        st.warning("Please select genes in the Input tab first.")
        return

    st.success(f"Ready to analyze {len(genes)} genes")

    # Organism selection
    organism = st.selectbox(
        "Organism",
        options=list(KEGG_ORGANISMS.keys()),
        index=0
    )

    # Analysis options
    st.subheader("Analysis Options")

    col1, col2 = st.columns(2)

    with col1:
        run_go_bp = st.checkbox("GO Biological Process", value=True)
        run_go_mf = st.checkbox("GO Molecular Function", value=True)
        run_go_cc = st.checkbox("GO Cellular Component", value=True)

    with col2:
        run_kegg = st.checkbox("KEGG Pathways", value=True)
        pvalue_cutoff = st.slider("P-value cutoff", 0.01, 0.10, 0.05, 0.01)

    # Run analysis
    st.divider()
    if st.button("🚀 Run Enrichment Analysis", type="primary", width="stretch"):
        background = st.session_state.get('enrichment_background')

        with st.spinner("Running enrichment analysis with g:Profiler..."):
            # Map organism name to g:Profiler code
            org_map = {
                'Arabidopsis thaliana': 'athaliana',
                'Medicago truncatula': 'mtruncatula',
                'Homo sapiens': 'hsapiens',
                'Mus musculus': 'mmusculus',
                'Oryza sativa': 'osativa',
                'Zea mays': 'zmays',
                'Glycine max': 'gmax',
                'Drosophila melanogaster': 'dmelanogaster',
                'Caenorhabditis elegans': 'celegans',
            }
            org_code = org_map.get(organism, 'athaliana')

            # Build sources list
            sources = []
            if run_go_bp:
                sources.append('GO:BP')
            if run_go_mf:
                sources.append('GO:MF')
            if run_go_cc:
                sources.append('GO:CC')
            if run_kegg:
                sources.append('KEGG')

            result = run_gprofiler_enrichment(
                genes=genes,
                organism=org_code,
                sources=sources,
                background=background
            )

            if result['status'] == 'success':
                st.session_state.enrichment_results = result
                st.success("Enrichment analysis complete!")

                # Summary
                if result.get('results'):
                    for source, df in result['results'].items():
                        if isinstance(df, pd.DataFrame) and len(df) > 0:
                            st.metric(source, f"{len(df)} terms")
                else:
                    st.info("No significant enrichment found")
            else:
                st.error(f"Analysis failed: {result.get('error', 'Unknown error')}")


def render_enrichment_results():
    """Render enrichment results section"""
    st.subheader("Enrichment Results")

    if st.session_state.get('enrichment_results') is None:
        st.info("No results yet. Run analysis first.")
        return

    results = st.session_state.enrichment_results.get('results', {})

    if not results:
        st.warning("No significant enrichment found")
        return

    # Select source
    source = st.selectbox(
        "Select Category",
        options=list(results.keys())
    )

    df = results[source]

    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)

    if df.empty:
        st.warning(f"No results for {source}")
        return

    # Summary metrics
    col1, col2, col3 = st.columns(3)

    with col1:
        st.metric("Total Terms", len(df))

    with col2:
        if 'p_value' in df.columns:
            sig = (df['p_value'] < 0.05).sum()
        elif 'pvalue' in df.columns:
            sig = (df['pvalue'] < 0.05).sum()
        else:
            sig = len(df)
        st.metric("Significant (p<0.05)", sig)

    with col3:
        if 'intersection_size' in df.columns:
            total_genes = df['intersection_size'].sum()
        elif 'Count' in df.columns:
            total_genes = df['Count'].sum()
        else:
            total_genes = "N/A"
        st.metric("Total Gene Hits", total_genes)

    # Visualization
    st.subheader("Visualization")

    viz_type = st.selectbox(
        "Plot Type",
        options=["Dot Plot", "Bar Plot", "Table"]
    )

    # Determine column names (varies by method)
    term_col = next((c for c in ['name', 'Description', 'term_name'] if c in df.columns), df.columns[0])
    pval_col = next((c for c in ['p_value', 'pvalue', 'p.adjust'] if c in df.columns), None)
    count_col = next((c for c in ['intersection_size', 'Count', 'gene_count'] if c in df.columns), None)
    ratio_col = next((c for c in ['precision', 'GeneRatio', 'gene_ratio'] if c in df.columns), None)

    top_n = st.slider("Number of terms to show", 5, 50, 20)

    if viz_type == "Dot Plot" and pval_col and count_col:
        # Sort by p-value
        plot_df = df.nsmallest(top_n, pval_col).copy()

        if ratio_col:
            fig = plot_enrichment_dotplot(
                plot_df,
                term_col=term_col,
                pval_col=pval_col,
                count_col=count_col,
                generatio_col=ratio_col,
                top_n=top_n,
                title=f"{source} Enrichment"
            )
        else:
            # Create simple dot plot with plotly
            import plotly.express as px

            plot_df['neg_log10_p'] = -np.log10(plot_df[pval_col].clip(lower=1e-50))

            fig = px.scatter(
                plot_df,
                x='neg_log10_p',
                y=term_col,
                size=count_col,
                color='neg_log10_p',
                color_continuous_scale='Reds',
                title=f"{source} Enrichment"
            )
            fig.update_layout(
                yaxis=dict(categoryorder='total ascending'),
                height=max(400, top_n * 25)
            )

        st.plotly_chart(fig, width="stretch")

    elif viz_type == "Bar Plot" and pval_col:
        import plotly.express as px

        plot_df = df.nsmallest(top_n, pval_col).copy()
        plot_df['neg_log10_p'] = -np.log10(plot_df[pval_col].clip(lower=1e-50))

        fig = px.bar(
            plot_df,
            x='neg_log10_p',
            y=term_col,
            orientation='h',
            color='neg_log10_p',
            color_continuous_scale='Reds',
            title=f"{source} Enrichment"
        )
        fig.update_layout(
            yaxis=dict(categoryorder='total ascending'),
            height=max(400, top_n * 25)
        )
        st.plotly_chart(fig, width="stretch")

    else:  # Table
        st.dataframe(df.head(top_n), width="stretch")

    # Download
    st.divider()
    csv = df.to_csv(index=False)
    st.download_button(
        f"📥 Download {source} Results",
        csv,
        f"{source}_enrichment.csv",
        "text/csv"
    )

    # Download all results
    if len(results) > 1:
        if st.button("📥 Download All Results"):
            import io
            import zipfile

            buffer = io.BytesIO()
            with zipfile.ZipFile(buffer, 'w', zipfile.ZIP_DEFLATED) as zf:
                for src, src_df in results.items():
                    if isinstance(src_df, pd.DataFrame):
                        zf.writestr(f"{src}_enrichment.csv", src_df.to_csv(index=False))

            st.download_button(
                "📥 Download ZIP",
                buffer.getvalue(),
                "enrichment_results.zip",
                "application/zip"
            )
