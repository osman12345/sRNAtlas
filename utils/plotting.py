"""
Plotting utilities for sRNAtlas
Interactive visualizations using Plotly
"""
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from typing import List, Dict, Optional, Tuple, Union
from scipy import stats

# =============================================================================
# COLOR PALETTES FOR CONSISTENT VISUALIZATION
# =============================================================================

# Professional color palette for categorical data (colorblind-friendly)
CATEGORICAL_COLORS = [
    '#E64B35',  # Red
    '#4DBBD5',  # Cyan
    '#00A087',  # Teal
    '#3C5488',  # Dark blue
    '#F39B7F',  # Salmon
    '#8491B4',  # Lavender
    '#91D1C2',  # Mint
    '#DC0000',  # Bright red
    '#7E6148',  # Brown
    '#B09C85',  # Tan
]

# Alternative vibrant palette
VIBRANT_COLORS = [
    '#E41A1C',  # Red
    '#377EB8',  # Blue
    '#4DAF4A',  # Green
    '#984EA3',  # Purple
    '#FF7F00',  # Orange
    '#FFFF33',  # Yellow
    '#A65628',  # Brown
    '#F781BF',  # Pink
]

# Paired colors for two-group comparisons
PAIRED_COLORS = {
    'Control': '#3C5488',
    'Treatment': '#E64B35',
    'WT': '#3C5488',
    'Mutant': '#E64B35',
    'Normal': '#4DBBD5',
    'Tumor': '#E64B35',
    'Untreated': '#00A087',
    'Treated': '#F39B7F',
}

# Significance colors
SIGNIFICANCE_COLORS = {
    'Up': '#E64B35',       # Red for upregulated
    'Down': '#3C5488',     # Blue for downregulated
    'Not Significant': '#CCCCCC',  # Gray for non-significant
}


def plot_read_length_distribution(
    lengths: List[int],
    title: str = "Read Length Distribution",
    expected_ranges: Optional[Dict[str, Tuple[int, int]]] = None
) -> go.Figure:
    """
    Create interactive read length distribution plot

    Args:
        lengths: List of read lengths
        title: Plot title
        expected_ranges: Dict of RNA type -> (min, max) length ranges

    Returns:
        Plotly figure
    """
    fig = go.Figure()

    # Histogram
    fig.add_trace(go.Histogram(
        x=lengths,
        nbinsx=50,
        name='Read Lengths',
        marker_color='#1E88E5',
        opacity=0.7
    ))

    # Add expected ranges if provided
    if expected_ranges:
        colors = px.colors.qualitative.Set2
        for i, (rna_type, (min_len, max_len)) in enumerate(expected_ranges.items()):
            fig.add_vrect(
                x0=min_len, x1=max_len,
                fillcolor=colors[i % len(colors)],
                opacity=0.2,
                line_width=0,
                annotation_text=rna_type,
                annotation_position="top left"
            )

    fig.update_layout(
        title=title,
        xaxis_title="Read Length (nt)",
        yaxis_title="Count",
        template="plotly_white",
        showlegend=True
    )

    return fig


def plot_quality_distribution(
    qualities: List[float],
    title: str = "Quality Score Distribution"
) -> go.Figure:
    """
    Create quality score distribution plot

    Args:
        qualities: List of quality scores
        title: Plot title

    Returns:
        Plotly figure
    """
    fig = go.Figure()

    fig.add_trace(go.Histogram(
        x=qualities,
        nbinsx=40,
        name='Quality Scores',
        marker_color='#43A047',
        opacity=0.7
    ))

    # Add quality threshold lines
    fig.add_vline(x=20, line_dash="dash", line_color="orange",
                  annotation_text="Q20")
    fig.add_vline(x=30, line_dash="dash", line_color="green",
                  annotation_text="Q30")

    fig.update_layout(
        title=title,
        xaxis_title="Phred Quality Score",
        yaxis_title="Count",
        template="plotly_white"
    )

    return fig


def plot_rna_type_distribution(
    annotations: pd.DataFrame,
    count_matrix: Optional[pd.DataFrame] = None,
    column: str = 'RNA_category'
) -> go.Figure:
    """
    Create RNA type distribution plot (by count or by reads)

    Args:
        annotations: DataFrame with RNA annotations
        count_matrix: Optional count matrix for read-based distribution
        column: Column to use for categorization

    Returns:
        Plotly figure
    """
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=["By Number of RNAs", "By Total Reads"],
        specs=[[{"type": "pie"}, {"type": "pie"}]]
    )

    # Distribution by number of RNAs
    type_counts = annotations[column].value_counts()

    fig.add_trace(
        go.Pie(
            labels=type_counts.index,
            values=type_counts.values,
            hole=0.4,
            name="RNA Count"
        ),
        row=1, col=1
    )

    # Distribution by reads (if count matrix provided)
    if count_matrix is not None:
        # Match annotations to counts
        merged = annotations.set_index('full_id' if 'full_id' in annotations.columns else annotations.index)
        sample_cols = [c for c in count_matrix.columns if c not in annotations.columns]

        read_counts = {}
        for rna_type in type_counts.index:
            mask = merged[column] == rna_type
            rnas_of_type = merged[mask].index.intersection(count_matrix.index)
            total_reads = count_matrix.loc[rnas_of_type, sample_cols].sum().sum()
            read_counts[rna_type] = total_reads

        fig.add_trace(
            go.Pie(
                labels=list(read_counts.keys()),
                values=list(read_counts.values()),
                hole=0.4,
                name="Reads"
            ),
            row=1, col=2
        )

    fig.update_layout(
        title="RNA Type Distribution",
        template="plotly_white"
    )

    return fig


def plot_pca(
    data: pd.DataFrame,
    metadata: pd.DataFrame,
    color_by: str = 'Treatment',
    shape_by: Optional[str] = None,
    n_components: int = 2,
    show_labels: bool = True,
    marker_size: int = 12
) -> go.Figure:
    """
    Create interactive PCA plot with improved visualization

    Args:
        data: Normalized count matrix (samples x genes or genes x samples)
        metadata: Sample metadata
        color_by: Column to use for coloring
        shape_by: Column to use for marker shapes
        n_components: Number of PCA components
        show_labels: Whether to show sample labels
        marker_size: Size of markers

    Returns:
        Plotly figure
    """
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    # Prepare data (transpose if needed - want samples as rows)
    data_for_pca = data.copy()
    if data_for_pca.shape[0] > data_for_pca.shape[1]:
        # Rows are genes, transpose to get samples as rows
        data_for_pca = data_for_pca.T

    # Handle missing values
    data_for_pca = data_for_pca.fillna(0)

    # Standardize
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data_for_pca)

    # PCA
    pca = PCA(n_components=min(n_components, min(data_scaled.shape)))
    pca_result = pca.fit_transform(data_scaled)

    # Create DataFrame for plotting
    plot_df = pd.DataFrame(
        pca_result,
        columns=[f'PC{i+1}' for i in range(pca_result.shape[1])],
        index=data_for_pca.index
    )
    plot_df['Sample'] = plot_df.index

    # Match metadata to samples
    if metadata is not None and len(metadata) > 0:
        # Determine the sample ID column in metadata
        if 'SampleID' in metadata.columns:
            meta_indexed = metadata.set_index('SampleID')
        elif 'sample' in metadata.columns:
            meta_indexed = metadata.set_index('sample')
        else:
            meta_indexed = metadata

        # Add metadata columns
        for col in metadata.columns:
            if col.lower() in ['sampleid', 'sample']:
                continue
            try:
                # Match samples between plot_df and metadata
                matched_values = []
                for sample in plot_df.index:
                    if sample in meta_indexed.index:
                        matched_values.append(meta_indexed.loc[sample, col])
                    else:
                        matched_values.append('Unknown')
                plot_df[col] = matched_values
            except Exception:
                continue

    # Determine color column
    color_col = None
    if color_by and color_by in plot_df.columns:
        color_col = color_by

    # Get unique groups for coloring
    if color_col:
        unique_groups = plot_df[color_col].unique()
        n_groups = len(unique_groups)

        # Assign colors - use predefined colors for common names, otherwise use palette
        color_map = {}
        palette_idx = 0
        for group in unique_groups:
            group_str = str(group)
            if group_str in PAIRED_COLORS:
                color_map[group] = PAIRED_COLORS[group_str]
            else:
                color_map[group] = CATEGORICAL_COLORS[palette_idx % len(CATEGORICAL_COLORS)]
                palette_idx += 1

    # Create figure
    fig = go.Figure()

    if color_col:
        # Plot each group separately for better control
        for group in unique_groups:
            group_data = plot_df[plot_df[color_col] == group]

            fig.add_trace(go.Scatter(
                x=group_data['PC1'],
                y=group_data['PC2'],
                mode='markers+text' if show_labels else 'markers',
                name=str(group),
                text=group_data['Sample'] if show_labels else None,
                textposition='top center',
                textfont=dict(size=10),
                marker=dict(
                    size=marker_size,
                    color=color_map[group],
                    line=dict(width=1, color='white'),
                    opacity=0.85
                ),
                hovertemplate=(
                    '<b>%{text}</b><br>' +
                    f'{color_col}: {group}<br>' +
                    'PC1: %{x:.2f}<br>' +
                    'PC2: %{y:.2f}<extra></extra>'
                )
            ))
    else:
        # Single group - no color distinction
        fig.add_trace(go.Scatter(
            x=plot_df['PC1'],
            y=plot_df['PC2'],
            mode='markers+text' if show_labels else 'markers',
            name='Samples',
            text=plot_df['Sample'] if show_labels else None,
            textposition='top center',
            textfont=dict(size=10),
            marker=dict(
                size=marker_size,
                color=CATEGORICAL_COLORS[0],
                line=dict(width=1, color='white'),
                opacity=0.85
            ),
            hovertemplate=(
                '<b>%{text}</b><br>' +
                'PC1: %{x:.2f}<br>' +
                'PC2: %{y:.2f}<extra></extra>'
            )
        ))

    # Calculate variance explained
    var_pc1 = pca.explained_variance_ratio_[0] * 100
    var_pc2 = pca.explained_variance_ratio_[1] * 100 if len(pca.explained_variance_ratio_) > 1 else 0

    # Update layout
    fig.update_layout(
        title=dict(
            text=f"<b>PCA Plot</b><br><sup>PC1: {var_pc1:.1f}% | PC2: {var_pc2:.1f}% variance explained</sup>",
            x=0.5,
            xanchor='center'
        ),
        xaxis_title=f"PC1 ({var_pc1:.1f}% variance)",
        yaxis_title=f"PC2 ({var_pc2:.1f}% variance)",
        template="plotly_white",
        legend=dict(
            title=color_col if color_col else "",
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99,
            bgcolor='rgba(255,255,255,0.8)',
            bordercolor='rgba(0,0,0,0.1)',
            borderwidth=1
        ),
        width=800,
        height=600,
        hovermode='closest'
    )

    # Make axes equal scale for proper PCA representation
    fig.update_xaxes(
        showgrid=True,
        gridwidth=1,
        gridcolor='rgba(0,0,0,0.1)',
        zeroline=True,
        zerolinewidth=1,
        zerolinecolor='rgba(0,0,0,0.2)'
    )
    fig.update_yaxes(
        showgrid=True,
        gridwidth=1,
        gridcolor='rgba(0,0,0,0.1)',
        zeroline=True,
        zerolinewidth=1,
        zerolinecolor='rgba(0,0,0,0.2)',
        scaleanchor="x",
        scaleratio=1
    )

    return fig


def plot_volcano(
    de_results: pd.DataFrame,
    lfc_col: str = 'log2FoldChange',
    pval_col: str = 'pvalue',
    padj_col: str = 'padj',
    lfc_threshold: float = 0.585,
    padj_threshold: float = 0.05,
    label_col: Optional[str] = None,
    n_labels: int = 10,
    title: str = "Volcano Plot"
) -> go.Figure:
    """
    Create interactive volcano plot with improved visualization

    Args:
        de_results: DataFrame with DE results
        lfc_col: Column name for log2 fold change
        pval_col: Column name for p-value
        padj_col: Column name for adjusted p-value
        lfc_threshold: LFC threshold for significance
        padj_threshold: FDR threshold for significance
        label_col: Column to use for labels
        n_labels: Number of top genes to label
        title: Plot title

    Returns:
        Plotly figure
    """
    # Prepare data
    df = de_results.copy()
    df['neg_log10_pval'] = -np.log10(df[pval_col].clip(lower=1e-300))

    # Categorize significance
    df['significance'] = 'Not Significant'
    df.loc[(df[padj_col] < padj_threshold) & (df[lfc_col] > lfc_threshold), 'significance'] = 'Up'
    df.loc[(df[padj_col] < padj_threshold) & (df[lfc_col] < -lfc_threshold), 'significance'] = 'Down'

    # Count significant
    n_up = (df['significance'] == 'Up').sum()
    n_down = (df['significance'] == 'Down').sum()
    n_total = len(df)

    # Create figure
    fig = go.Figure()

    # Plot order: non-significant first (background), then significant on top
    plot_order = ['Not Significant', 'Down', 'Up']

    for sig_type in plot_order:
        subset = df[df['significance'] == sig_type]
        if len(subset) == 0:
            continue

        # Marker settings based on significance
        if sig_type == 'Not Significant':
            marker_settings = dict(
                color='rgba(180, 180, 180, 0.4)',
                size=5,
                line=dict(width=0)
            )
        elif sig_type == 'Up':
            marker_settings = dict(
                color=SIGNIFICANCE_COLORS['Up'],
                size=8,
                line=dict(width=1, color='white'),
                opacity=0.8
            )
        else:  # Down
            marker_settings = dict(
                color=SIGNIFICANCE_COLORS['Down'],
                size=8,
                line=dict(width=1, color='white'),
                opacity=0.8
            )

        fig.add_trace(go.Scatter(
            x=subset[lfc_col],
            y=subset['neg_log10_pval'],
            mode='markers',
            name=f"{sig_type} ({len(subset)})",
            marker=marker_settings,
            text=subset.index if label_col is None else subset[label_col],
            hovertemplate=(
                '<b>%{text}</b><br>' +
                'Log2FC: %{x:.3f}<br>' +
                '-log10(p): %{y:.2f}<br>' +
                f'Status: {sig_type}<extra></extra>'
            )
        ))

    # Add threshold lines
    fig.add_hline(
        y=-np.log10(padj_threshold),
        line_dash="dash",
        line_color="rgba(100, 100, 100, 0.5)",
        annotation_text=f"FDR={padj_threshold}",
        annotation_position="right"
    )
    fig.add_vline(x=lfc_threshold, line_dash="dash", line_color="rgba(100, 100, 100, 0.5)")
    fig.add_vline(x=-lfc_threshold, line_dash="dash", line_color="rgba(100, 100, 100, 0.5)")

    # Add labels for top genes
    if n_labels > 0:
        sig_genes = df[df['significance'] != 'Not Significant'].nlargest(n_labels, 'neg_log10_pval')
        for idx, row in sig_genes.iterrows():
            label_text = str(idx)[:25] if label_col is None else str(row[label_col])[:25]
            fig.add_annotation(
                x=row[lfc_col],
                y=row['neg_log10_pval'],
                text=label_text,
                showarrow=True,
                arrowhead=2,
                arrowsize=1,
                arrowwidth=1,
                arrowcolor='rgba(0,0,0,0.5)',
                font=dict(size=9),
                bgcolor='rgba(255,255,255,0.8)',
                borderpad=2
            )

    fig.update_layout(
        title=dict(
            text=f"<b>{title}</b><br><sup>↑ Up: {n_up} | ↓ Down: {n_down} | Total: {n_total}</sup>",
            x=0.5,
            xanchor='center'
        ),
        xaxis_title="Log2 Fold Change",
        yaxis_title="-Log10 P-value",
        template="plotly_white",
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99,
            bgcolor='rgba(255,255,255,0.9)',
            bordercolor='rgba(0,0,0,0.1)',
            borderwidth=1
        ),
        width=800,
        height=600,
        hovermode='closest'
    )

    # Update axes
    fig.update_xaxes(
        showgrid=True,
        gridwidth=1,
        gridcolor='rgba(0,0,0,0.05)',
        zeroline=True,
        zerolinewidth=1,
        zerolinecolor='rgba(0,0,0,0.2)'
    )
    fig.update_yaxes(
        showgrid=True,
        gridwidth=1,
        gridcolor='rgba(0,0,0,0.05)'
    )

    return fig


def plot_ma(
    de_results: pd.DataFrame,
    basemean_col: str = 'baseMean',
    lfc_col: str = 'log2FoldChange',
    padj_col: str = 'padj',
    padj_threshold: float = 0.05,
    lfc_threshold: float = 0.585,
    title: str = "MA Plot"
) -> go.Figure:
    """
    Create interactive MA plot with improved visualization

    Args:
        de_results: DataFrame with DE results
        basemean_col: Column name for mean expression
        lfc_col: Column name for log2 fold change
        padj_col: Column name for adjusted p-value
        padj_threshold: FDR threshold for significance
        lfc_threshold: LFC threshold for coloring
        title: Plot title

    Returns:
        Plotly figure
    """
    df = de_results.copy()
    df['log10_mean'] = np.log10(df[basemean_col].clip(lower=1))

    # Categorize by direction
    df['status'] = 'Not Significant'
    df.loc[(df[padj_col] < padj_threshold) & (df[lfc_col] > lfc_threshold), 'status'] = 'Up'
    df.loc[(df[padj_col] < padj_threshold) & (df[lfc_col] < -lfc_threshold), 'status'] = 'Down'

    n_up = (df['status'] == 'Up').sum()
    n_down = (df['status'] == 'Down').sum()

    fig = go.Figure()

    # Plot order for layering
    for status in ['Not Significant', 'Down', 'Up']:
        subset = df[df['status'] == status]
        if len(subset) == 0:
            continue

        if status == 'Not Significant':
            marker_settings = dict(
                color='rgba(180, 180, 180, 0.4)',
                size=5,
                line=dict(width=0)
            )
        elif status == 'Up':
            marker_settings = dict(
                color=SIGNIFICANCE_COLORS['Up'],
                size=7,
                line=dict(width=0.5, color='white'),
                opacity=0.8
            )
        else:  # Down
            marker_settings = dict(
                color=SIGNIFICANCE_COLORS['Down'],
                size=7,
                line=dict(width=0.5, color='white'),
                opacity=0.8
            )

        fig.add_trace(go.Scatter(
            x=subset['log10_mean'],
            y=subset[lfc_col],
            mode='markers',
            name=f'{status} ({len(subset)})',
            marker=marker_settings,
            text=subset.index,
            hovertemplate=(
                '<b>%{text}</b><br>' +
                'Mean (log10): %{x:.2f}<br>' +
                'Log2FC: %{y:.3f}<extra></extra>'
            )
        ))

    # Add horizontal line at y=0
    fig.add_hline(y=0, line_color="rgba(0,0,0,0.5)", line_width=1.5)

    fig.update_layout(
        title=dict(
            text=f"<b>{title}</b><br><sup>↑ Up: {n_up} | ↓ Down: {n_down}</sup>",
            x=0.5,
            xanchor='center'
        ),
        xaxis_title="Log10 Mean Expression",
        yaxis_title="Log2 Fold Change",
        template="plotly_white",
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99,
            bgcolor='rgba(255,255,255,0.9)'
        ),
        width=800,
        height=600
    )

    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='rgba(0,0,0,0.05)')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='rgba(0,0,0,0.05)')

    return fig


def plot_heatmap(
    data: pd.DataFrame,
    row_labels: Optional[List[str]] = None,
    col_labels: Optional[List[str]] = None,
    title: str = "Heatmap",
    cluster_rows: bool = True,
    cluster_cols: bool = True,
    color_scale: str = "RdBu_r",
    center_zero: bool = True
) -> go.Figure:
    """
    Create interactive heatmap with improved visualization

    Args:
        data: Matrix data
        row_labels: Labels for rows
        col_labels: Labels for columns
        title: Plot title
        cluster_rows: Whether to cluster rows
        cluster_cols: Whether to cluster columns
        color_scale: Plotly color scale
        center_zero: Whether to center the color scale at zero

    Returns:
        Plotly figure
    """
    from scipy.cluster.hierarchy import linkage, leaves_list

    plot_data = data.copy()

    # Handle missing values
    plot_data = plot_data.fillna(0)

    # Cluster if requested
    if cluster_rows and len(data) > 2:
        try:
            row_linkage = linkage(plot_data, method='ward')
            row_order = leaves_list(row_linkage)
            plot_data = plot_data.iloc[row_order]
        except Exception:
            pass  # Skip clustering if it fails

    if cluster_cols and len(data.columns) > 2:
        try:
            col_linkage = linkage(plot_data.T, method='ward')
            col_order = leaves_list(col_linkage)
            plot_data = plot_data.iloc[:, col_order]
        except Exception:
            pass

    # Determine color scale centering
    z_min, z_max = plot_data.values.min(), plot_data.values.max()
    if center_zero and z_min < 0 and z_max > 0:
        z_abs_max = max(abs(z_min), abs(z_max))
        z_min, z_max = -z_abs_max, z_abs_max

    fig = go.Figure(data=go.Heatmap(
        z=plot_data.values,
        x=col_labels if col_labels else plot_data.columns.tolist(),
        y=row_labels if row_labels else plot_data.index.tolist(),
        colorscale=color_scale,
        zmin=z_min,
        zmax=z_max,
        hovertemplate='<b>%{y}</b><br>Sample: %{x}<br>Value: %{z:.3f}<extra></extra>',
        colorbar=dict(
            title='Z-score' if center_zero else 'Value',
            titleside='right'
        )
    ))

    fig.update_layout(
        title=dict(text=f"<b>{title}</b>", x=0.5, xanchor='center'),
        template="plotly_white",
        height=max(500, len(plot_data) * 18 + 150),
        width=max(600, len(plot_data.columns) * 60 + 200),
        xaxis=dict(tickangle=-45),
        yaxis=dict(tickfont=dict(size=10))
    )

    return fig


def plot_enrichment_dotplot(
    enrichment_results: pd.DataFrame,
    term_col: str = 'Description',
    pval_col: str = 'p.adjust',
    count_col: str = 'Count',
    generatio_col: str = 'GeneRatio',
    top_n: int = 20,
    title: str = "GO Enrichment"
) -> go.Figure:
    """
    Create enrichment dot plot with improved visualization (like clusterProfiler)

    Args:
        enrichment_results: DataFrame with enrichment results
        term_col: Column for term names
        pval_col: Column for adjusted p-values
        count_col: Column for gene counts
        generatio_col: Column for gene ratio
        top_n: Number of top terms to show
        title: Plot title

    Returns:
        Plotly figure
    """
    # Get top terms
    df = enrichment_results.nsmallest(top_n, pval_col).copy()

    if len(df) == 0:
        fig = go.Figure()
        fig.add_annotation(text="No significant enrichment results",
                          xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False)
        return fig

    # Parse gene ratio if needed
    if generatio_col in df.columns:
        if df[generatio_col].dtype == object:
            df['ratio'] = df[generatio_col].apply(
                lambda x: eval(x) if isinstance(x, str) and '/' in x else float(x) if x else 0
            )
        else:
            df['ratio'] = df[generatio_col]
    else:
        df['ratio'] = 0.1  # Default if no ratio column

    df['neg_log10_padj'] = -np.log10(df[pval_col].clip(lower=1e-50))

    # Truncate long term names
    df['term_display'] = df[term_col].apply(lambda x: x[:50] + '...' if len(str(x)) > 50 else x)

    # Custom colorscale (Viridis-like but more saturated)
    colorscale = [
        [0.0, '#440154'],
        [0.25, '#3B528B'],
        [0.5, '#21918C'],
        [0.75, '#5DC863'],
        [1.0, '#FDE725']
    ]

    fig = go.Figure()

    # Calculate size scaling
    max_count = df[count_col].max() if count_col in df.columns else 10
    min_count = df[count_col].min() if count_col in df.columns else 1

    fig.add_trace(go.Scatter(
        x=df['ratio'],
        y=df['term_display'],
        mode='markers',
        marker=dict(
            size=df[count_col] if count_col in df.columns else 10,
            sizemode='area',
            sizeref=2. * max_count / (35.**2),
            sizemin=8,
            color=df['neg_log10_padj'],
            colorscale=colorscale,
            showscale=True,
            colorbar=dict(
                title=dict(text='-log10(FDR)', side='right'),
                thickness=15
            ),
            line=dict(width=1, color='white')
        ),
        text=[f"Count: {c}<br>FDR: {p:.2e}" for c, p in zip(df[count_col], df[pval_col])],
        hovertemplate=(
            '<b>%{y}</b><br>' +
            'Gene Ratio: %{x:.3f}<br>' +
            '%{text}<extra></extra>'
        )
    ))

    fig.update_layout(
        title=dict(
            text=f"<b>{title}</b><br><sup>Top {len(df)} enriched terms</sup>",
            x=0.5,
            xanchor='center'
        ),
        xaxis_title="Gene Ratio",
        yaxis_title="",
        template="plotly_white",
        height=max(450, len(df) * 28 + 100),
        width=800,
        yaxis=dict(
            categoryorder='total ascending',
            tickfont=dict(size=10)
        ),
        margin=dict(l=250)  # More space for term names
    )

    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='rgba(0,0,0,0.05)')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='rgba(0,0,0,0.05)')

    return fig


def plot_sample_correlation(
    count_matrix: pd.DataFrame,
    metadata: Optional[pd.DataFrame] = None,
    method: str = 'spearman',
    title: str = "Sample Correlation"
) -> go.Figure:
    """
    Create sample correlation heatmap with improved visualization

    Args:
        count_matrix: Count matrix (genes x samples)
        metadata: Optional metadata for annotation
        method: Correlation method ('pearson' or 'spearman')
        title: Plot title

    Returns:
        Plotly figure
    """
    # Calculate correlation
    corr = count_matrix.corr(method=method)

    # Custom colorscale that emphasizes high correlations
    colorscale = [
        [0.0, '#3C5488'],   # Low correlation - blue
        [0.5, '#FFFFFF'],   # Mid - white
        [0.7, '#FFF5F0'],   # Light pink
        [0.85, '#FC9272'],  # Medium salmon
        [1.0, '#E64B35']    # High correlation - red
    ]

    # Add annotation text
    annotations = []
    for i, row in enumerate(corr.index):
        for j, col in enumerate(corr.columns):
            annotations.append(dict(
                x=col,
                y=row,
                text=f'{corr.iloc[i, j]:.2f}',
                font=dict(size=9, color='black' if 0.3 < corr.iloc[i, j] < 0.9 else 'white'),
                showarrow=False
            ))

    fig = go.Figure(data=go.Heatmap(
        z=corr.values,
        x=corr.columns,
        y=corr.index,
        colorscale=colorscale,
        zmin=0,
        zmax=1,
        hovertemplate='<b>%{y}</b> vs <b>%{x}</b><br>Correlation: %{z:.3f}<extra></extra>',
        colorbar=dict(
            title=f'{method.capitalize()}<br>Correlation',
            titleside='right'
        )
    ))

    # Add correlation values as text (only if not too many samples)
    if len(corr) <= 15:
        fig.update_layout(annotations=annotations)

    n_samples = len(corr)
    fig.update_layout(
        title=dict(
            text=f"<b>{title}</b><br><sup>{method.capitalize()} correlation | {n_samples} samples</sup>",
            x=0.5,
            xanchor='center'
        ),
        template="plotly_white",
        width=max(600, n_samples * 50 + 150),
        height=max(600, n_samples * 50 + 150),
        xaxis=dict(tickangle=-45, side='bottom'),
        yaxis=dict(autorange='reversed')  # To match traditional heatmap orientation
    )

    return fig
