"""
Multi-Sample Comparison for sRNAtlas
Compare expression across multiple conditions (>2 groups)
"""
import streamlit as st
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple
import plotly.express as px
import plotly.graph_objects as go
from itertools import combinations

try:
    from scipy import stats
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False


def run_anova_analysis(
    count_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    condition_col: str = 'condition',
    min_count: int = 10
) -> pd.DataFrame:
    """
    Run one-way ANOVA for each feature across multiple groups
    
    Args:
        count_matrix: Expression matrix (features x samples)
        metadata: Sample metadata with condition column
        condition_col: Column name for conditions
        min_count: Minimum total count to include feature
        
    Returns:
        DataFrame with ANOVA results
    """
    if not SCIPY_AVAILABLE:
        raise ImportError("scipy required for ANOVA analysis")
    
    # Ensure sample order matches
    samples = [s for s in count_matrix.columns if s in metadata['sample'].values]
    count_matrix = count_matrix[samples]
    metadata = metadata[metadata['sample'].isin(samples)]
    
    # Get unique conditions
    conditions = metadata[condition_col].unique()
    
    if len(conditions) < 2:
        raise ValueError("Need at least 2 conditions for comparison")
    
    results = []
    
    for feature in count_matrix.index:
        row = count_matrix.loc[feature]
        
        # Skip low-count features
        if row.sum() < min_count:
            continue
        
        # Get values for each condition
        groups = []
        for cond in conditions:
            cond_samples = metadata[metadata[condition_col] == cond]['sample']
            values = row[cond_samples].values
            groups.append(values)
        
        # Run ANOVA
        try:
            f_stat, p_value = stats.f_oneway(*groups)
            
            # Calculate group means
            means = {cond: row[metadata[metadata[condition_col] == cond]['sample']].mean() 
                    for cond in conditions}
            
            result = {
                'feature': feature,
                'f_statistic': f_stat,
                'p_value': p_value,
                'mean_expression': row.mean(),
                **{f'mean_{cond}': means[cond] for cond in conditions}
            }
            results.append(result)
            
        except Exception:
            continue
    
    results_df = pd.DataFrame(results)
    
    # Multiple testing correction (Benjamini-Hochberg)
    if not results_df.empty:
        results_df = results_df.sort_values('p_value')
        n = len(results_df)
        results_df['rank'] = range(1, n + 1)
        results_df['padj'] = results_df['p_value'] * n / results_df['rank']
        results_df['padj'] = results_df['padj'].clip(upper=1.0)
        results_df = results_df.drop('rank', axis=1)
    
    return results_df


def run_pairwise_comparisons(
    count_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    condition_col: str = 'condition',
    method: str = 'ttest'
) -> Dict[str, pd.DataFrame]:
    """
    Run pairwise comparisons between all condition pairs
    
    Returns:
        Dictionary mapping comparison name to results DataFrame
    """
    conditions = metadata[condition_col].unique()
    pairs = list(combinations(conditions, 2))
    
    results = {}
    
    for cond1, cond2 in pairs:
        comparison_name = f"{cond1}_vs_{cond2}"
        
        # Get samples for each condition
        samples1 = metadata[metadata[condition_col] == cond1]['sample']
        samples2 = metadata[metadata[condition_col] == cond2]['sample']
        
        pair_results = []
        
        for feature in count_matrix.index:
            vals1 = count_matrix.loc[feature, samples1].values
            vals2 = count_matrix.loc[feature, samples2].values
            
            mean1 = np.mean(vals1)
            mean2 = np.mean(vals2)
            
            # Log2 fold change
            log2fc = np.log2((mean2 + 1) / (mean1 + 1))
            
            # Statistical test
            if method == 'ttest':
                try:
                    _, p_value = stats.ttest_ind(vals1, vals2)
                except:
                    p_value = 1.0
            elif method == 'mannwhitney':
                try:
                    _, p_value = stats.mannwhitneyu(vals1, vals2)
                except:
                    p_value = 1.0
            else:
                p_value = 1.0
            
            pair_results.append({
                'feature': feature,
                f'mean_{cond1}': mean1,
                f'mean_{cond2}': mean2,
                'log2FoldChange': log2fc,
                'pvalue': p_value
            })
        
        df = pd.DataFrame(pair_results)
        
        # FDR correction
        if not df.empty:
            df = df.sort_values('pvalue')
            n = len(df)
            df['rank'] = range(1, n + 1)
            df['padj'] = df['pvalue'] * n / df['rank']
            df['padj'] = df['padj'].clip(upper=1.0)
            df = df.drop('rank', axis=1)
        
        results[comparison_name] = df
    
    return results


def plot_multi_group_boxplot(
    count_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    feature: str,
    condition_col: str = 'condition'
) -> go.Figure:
    """Create boxplot for a feature across multiple groups"""
    
    # Prepare data
    plot_data = []
    for _, row in metadata.iterrows():
        sample = row['sample']
        condition = row[condition_col]
        if sample in count_matrix.columns:
            plot_data.append({
                'Sample': sample,
                'Condition': condition,
                'Expression': count_matrix.loc[feature, sample]
            })
    
    df = pd.DataFrame(plot_data)
    
    fig = px.box(
        df,
        x='Condition',
        y='Expression',
        points='all',
        title=f"Expression of {feature}",
        color='Condition'
    )
    
    return fig


def plot_multi_group_heatmap(
    anova_results: pd.DataFrame,
    count_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    condition_col: str = 'condition',
    top_n: int = 50,
    padj_threshold: float = 0.05
) -> go.Figure:
    """Create heatmap of top significant features"""
    
    # Get significant features
    sig_results = anova_results[anova_results['padj'] < padj_threshold]
    sig_results = sig_results.nlargest(top_n, 'f_statistic')
    
    if sig_results.empty:
        return None
    
    # Get expression data
    features = sig_results['feature'].tolist()
    data = count_matrix.loc[features]
    
    # Log transform and z-score
    data_log = np.log2(data + 1)
    data_z = (data_log.T - data_log.mean(axis=1)) / data_log.std(axis=1)
    data_z = data_z.T
    
    # Order samples by condition
    sample_order = metadata.sort_values(condition_col)['sample'].tolist()
    sample_order = [s for s in sample_order if s in data_z.columns]
    data_z = data_z[sample_order]
    
    # Create heatmap
    fig = px.imshow(
        data_z,
        title=f"Top {len(features)} Significant Features (ANOVA p<{padj_threshold})",
        color_continuous_scale='RdBu_r',
        aspect='auto',
        labels=dict(color="Z-score")
    )
    
    fig.update_layout(height=600)
    
    return fig


def render_multi_comparison_ui(count_matrix: pd.DataFrame, metadata: pd.DataFrame):
    """Render multi-sample comparison UI"""
    st.subheader("📊 Multi-Group Comparison")
    
    if count_matrix is None or metadata is None:
        st.warning("Load count matrix and metadata first")
        return
    
    if not SCIPY_AVAILABLE:
        st.error("scipy required for statistical analysis")
        return
    
    # Settings
    col1, col2 = st.columns(2)
    
    with col1:
        condition_col = st.selectbox(
            "Condition column",
            [c for c in metadata.columns if c != 'sample']
        )
        
        conditions = metadata[condition_col].unique()
        st.info(f"Found {len(conditions)} groups: {', '.join(conditions)}")
    
    with col2:
        min_count = st.number_input("Minimum count filter", 1, 100, 10)
        padj_threshold = st.slider("FDR threshold", 0.01, 0.1, 0.05)
    
    # Run analysis
    if st.button("🔍 Run Multi-Group Analysis", type="primary"):
        with st.spinner("Running ANOVA analysis..."):
            # ANOVA
            anova_results = run_anova_analysis(
                count_matrix, metadata, condition_col, min_count
            )
            st.session_state.anova_results = anova_results
            
            # Pairwise
            pairwise_results = run_pairwise_comparisons(
                count_matrix, metadata, condition_col
            )
            st.session_state.pairwise_results = pairwise_results
            
            st.success("✅ Analysis complete!")
    
    # Display results
    if 'anova_results' in st.session_state:
        anova_results = st.session_state.anova_results
        
        # Summary
        sig_count = (anova_results['padj'] < padj_threshold).sum()
        st.metric("Significant Features", f"{sig_count} / {len(anova_results)}")
        
        # Results table
        st.markdown("#### ANOVA Results")
        display_df = anova_results[anova_results['padj'] < padj_threshold].head(50)
        st.dataframe(display_df, use_container_width=True)
        
        # Download
        csv = anova_results.to_csv(index=False)
        st.download_button(
            "📥 Download Full Results",
            csv,
            file_name="anova_results.csv",
            mime="text/csv"
        )
        
        # Heatmap
        st.markdown("#### Heatmap of Top Significant Features")
        heatmap = plot_multi_group_heatmap(
            anova_results, count_matrix, metadata, condition_col,
            top_n=30, padj_threshold=padj_threshold
        )
        if heatmap:
            st.plotly_chart(heatmap, use_container_width=True)
        
        # Pairwise comparisons
        if 'pairwise_results' in st.session_state:
            st.markdown("#### Pairwise Comparisons")
            
            for comp_name, comp_df in st.session_state.pairwise_results.items():
                with st.expander(f"📊 {comp_name}"):
                    sig = (comp_df['padj'] < padj_threshold).sum()
                    st.caption(f"Significant: {sig} features")
                    st.dataframe(comp_df[comp_df['padj'] < padj_threshold].head(20))
