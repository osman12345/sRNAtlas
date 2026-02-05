"""
Cluster Analysis for sRNAtlas
Hierarchical clustering and heatmap visualization
"""
import streamlit as st
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple
import plotly.express as px
import plotly.graph_objects as go
import plotly.figure_factory as ff

try:
    from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
    from scipy.spatial.distance import pdist
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

try:
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False


def run_hierarchical_clustering(
    data: pd.DataFrame,
    method: str = 'average',
    metric: str = 'euclidean',
    n_clusters: Optional[int] = None
) -> Dict:
    """
    Perform hierarchical clustering on expression data
    
    Args:
        data: Expression matrix (features x samples)
        method: Linkage method ('single', 'complete', 'average', 'ward')
        metric: Distance metric ('euclidean', 'correlation', 'cosine')
        n_clusters: Number of clusters to cut tree into
    
    Returns:
        Dictionary with clustering results
    """
    if not SCIPY_AVAILABLE:
        return {'error': 'scipy not installed'}
    
    results = {}
    
    # Standardize data
    if SKLEARN_AVAILABLE:
        scaler = StandardScaler()
        data_scaled = pd.DataFrame(
            scaler.fit_transform(data),
            index=data.index,
            columns=data.columns
        )
    else:
        # Simple z-score normalization
        data_scaled = (data - data.mean()) / data.std()
    
    # Cluster samples (columns)
    sample_dist = pdist(data_scaled.T, metric=metric)
    sample_linkage = linkage(sample_dist, method=method)
    results['sample_linkage'] = sample_linkage
    
    # Cluster features (rows)
    if len(data) > 1:
        feature_dist = pdist(data_scaled, metric=metric)
        feature_linkage = linkage(feature_dist, method=method)
        results['feature_linkage'] = feature_linkage
    
    # Cut tree if n_clusters specified
    if n_clusters:
        results['sample_clusters'] = fcluster(sample_linkage, n_clusters, criterion='maxclust')
        if 'feature_linkage' in results:
            results['feature_clusters'] = fcluster(feature_linkage, n_clusters, criterion='maxclust')
    
    results['data_scaled'] = data_scaled
    
    return results


def plot_clustered_heatmap(
    data: pd.DataFrame,
    cluster_results: Dict,
    title: str = "Clustered Heatmap",
    color_scale: str = "RdBu_r",
    show_dendrograms: bool = True
) -> go.Figure:
    """
    Create clustered heatmap with dendrograms
    """
    if not SCIPY_AVAILABLE:
        # Fallback to simple heatmap
        fig = px.imshow(
            data,
            title=title,
            color_continuous_scale=color_scale,
            aspect='auto'
        )
        return fig
    
    data_scaled = cluster_results.get('data_scaled', data)
    
    # Create dendrogram figure
    if show_dendrograms and 'sample_linkage' in cluster_results:
        fig = ff.create_dendrogram(
            data_scaled.T,
            orientation='bottom',
            labels=data_scaled.columns.tolist()
        )
        
        # Get reordered indices
        dendro_leaves = fig['layout']['xaxis']['ticktext']
        
        # Reorder data
        data_reordered = data_scaled[dendro_leaves]
        
        # Create heatmap
        heatmap = go.Heatmap(
            z=data_reordered.values,
            x=data_reordered.columns.tolist(),
            y=data_reordered.index.tolist(),
            colorscale=color_scale,
            colorbar=dict(title="Z-score")
        )
        
        fig.add_trace(heatmap)
        fig.update_layout(
            title=title,
            xaxis=dict(tickangle=45),
            height=600 + len(data) * 5
        )
    else:
        # Simple heatmap
        fig = px.imshow(
            data_scaled,
            title=title,
            color_continuous_scale=color_scale,
            aspect='auto',
            labels=dict(color="Z-score")
        )
        fig.update_layout(height=600)
    
    return fig


def plot_sample_dendrogram(cluster_results: Dict, labels: List[str] = None) -> go.Figure:
    """Plot sample dendrogram"""
    if not SCIPY_AVAILABLE or 'sample_linkage' not in cluster_results:
        return None
    
    # Create dendrogram
    fig = ff.create_dendrogram(
        cluster_results['data_scaled'].T,
        labels=labels or cluster_results['data_scaled'].columns.tolist(),
        orientation='bottom'
    )
    
    fig.update_layout(
        title="Sample Clustering Dendrogram",
        xaxis_title="Samples",
        yaxis_title="Distance",
        height=400
    )
    
    return fig


def run_pca_analysis(
    data: pd.DataFrame,
    n_components: int = 2
) -> Dict:
    """
    Run PCA on expression data
    """
    if not SKLEARN_AVAILABLE:
        return {'error': 'scikit-learn not installed'}
    
    # Standardize
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data.T)  # Transpose: samples as rows
    
    # PCA
    pca = PCA(n_components=min(n_components, data_scaled.shape[1]))
    pca_result = pca.fit_transform(data_scaled)
    
    # Create result DataFrame
    pca_df = pd.DataFrame(
        pca_result,
        columns=[f'PC{i+1}' for i in range(pca_result.shape[1])],
        index=data.columns
    )
    
    return {
        'pca_df': pca_df,
        'explained_variance': pca.explained_variance_ratio_,
        'loadings': pd.DataFrame(
            pca.components_.T,
            columns=[f'PC{i+1}' for i in range(pca.n_components_)],
            index=data.index
        )
    }


def plot_pca(
    pca_results: Dict,
    metadata: pd.DataFrame = None,
    color_by: str = None,
    title: str = "PCA Plot"
) -> go.Figure:
    """Create PCA scatter plot"""
    pca_df = pca_results['pca_df'].copy()
    var_explained = pca_results['explained_variance']
    
    # Add metadata if provided
    if metadata is not None and color_by:
        pca_df = pca_df.join(metadata[[color_by]])
    
    # Create plot
    if color_by and color_by in pca_df.columns:
        fig = px.scatter(
            pca_df,
            x='PC1',
            y='PC2',
            color=color_by,
            text=pca_df.index,
            title=title
        )
    else:
        fig = px.scatter(
            pca_df,
            x='PC1',
            y='PC2',
            text=pca_df.index,
            title=title
        )
    
    fig.update_traces(textposition='top center')
    fig.update_layout(
        xaxis_title=f"PC1 ({var_explained[0]*100:.1f}%)",
        yaxis_title=f"PC2 ({var_explained[1]*100:.1f}%)",
        height=500
    )
    
    return fig


def render_cluster_analysis_ui(count_matrix: pd.DataFrame, metadata: pd.DataFrame = None):
    """Render cluster analysis UI in Streamlit"""
    st.subheader("🔬 Cluster Analysis")
    
    if count_matrix is None or count_matrix.empty:
        st.warning("No count matrix available. Run counting first.")
        return
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("#### Clustering Settings")
        
        method = st.selectbox(
            "Linkage method",
            ['average', 'complete', 'single', 'ward'],
            help="Ward minimizes variance, average/complete for general use"
        )
        
        metric = st.selectbox(
            "Distance metric",
            ['euclidean', 'correlation', 'cosine'],
            help="Correlation often best for expression data"
        )
        
        n_clusters = st.slider(
            "Number of clusters (optional)",
            2, 10, 3,
            help="Cut dendrogram into this many clusters"
        )
    
    with col2:
        st.markdown("#### Visualization Options")
        
        top_n = st.slider(
            "Top variable features",
            10, 500, 50,
            help="Use most variable features for clustering"
        )
        
        color_scale = st.selectbox(
            "Color scale",
            ['RdBu_r', 'Viridis', 'Plasma', 'Blues', 'Reds']
        )
        
        show_dendro = st.checkbox("Show dendrograms", value=True)
    
    if st.button("🔍 Run Clustering", type="primary"):
        with st.spinner("Running cluster analysis..."):
            # Select top variable features
            variances = count_matrix.var(axis=1)
            top_features = variances.nlargest(top_n).index
            data_subset = count_matrix.loc[top_features]
            
            # Log transform
            data_log = np.log2(data_subset + 1)
            
            # Run clustering
            cluster_results = run_hierarchical_clustering(
                data_log,
                method=method,
                metric=metric,
                n_clusters=n_clusters
            )
            
            if 'error' in cluster_results:
                st.error(cluster_results['error'])
                return
            
            # Store results
            st.session_state.cluster_results = cluster_results
            st.session_state.cluster_data = data_log
            
            st.success("✅ Clustering complete!")
    
    # Display results
    if 'cluster_results' in st.session_state:
        cluster_results = st.session_state.cluster_results
        data_log = st.session_state.cluster_data
        
        # Heatmap
        st.markdown("#### Clustered Heatmap")
        fig = plot_clustered_heatmap(
            data_log,
            cluster_results,
            color_scale=color_scale,
            show_dendrograms=show_dendro
        )
        st.plotly_chart(fig, use_container_width=True)
        
        # Sample dendrogram
        if show_dendro:
            st.markdown("#### Sample Dendrogram")
            dendro_fig = plot_sample_dendrogram(cluster_results)
            if dendro_fig:
                st.plotly_chart(dendro_fig, use_container_width=True)
        
        # PCA
        st.markdown("#### PCA Analysis")
        pca_results = run_pca_analysis(data_log)
        
        if 'error' not in pca_results:
            color_by = None
            if metadata is not None and not metadata.empty:
                color_cols = [c for c in metadata.columns if c != 'sample']
                if color_cols:
                    color_by = st.selectbox("Color by", ['None'] + color_cols)
                    if color_by == 'None':
                        color_by = None
            
            pca_fig = plot_pca(pca_results, metadata, color_by)
            st.plotly_chart(pca_fig, use_container_width=True)
            
            # Variance explained
            var_exp = pca_results['explained_variance']
            st.caption(f"Variance explained: PC1={var_exp[0]*100:.1f}%, PC2={var_exp[1]*100:.1f}%")
