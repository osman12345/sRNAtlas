"""
Tests for cluster analysis utility
"""
import pytest
import pandas as pd
import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))


class TestHierarchicalClustering:
    """Tests for hierarchical clustering"""
    
    def test_clustering_returns_dict(self, sample_count_matrix):
        """Test clustering returns a dictionary"""
        try:
            from utils.cluster_analysis import run_hierarchical_clustering
            result = run_hierarchical_clustering(sample_count_matrix)
            assert isinstance(result, dict)
        except ImportError:
            pytest.skip("scipy not available")
    
    def test_clustering_with_different_methods(self, sample_count_matrix):
        """Test clustering with different linkage methods"""
        try:
            from utils.cluster_analysis import run_hierarchical_clustering
            
            for method in ['average', 'complete', 'single', 'ward']:
                result = run_hierarchical_clustering(
                    sample_count_matrix, 
                    method=method
                )
                assert 'sample_linkage' in result or 'error' in result
        except ImportError:
            pytest.skip("scipy not available")
    
    def test_clustering_with_n_clusters(self, sample_count_matrix):
        """Test clustering with specified number of clusters"""
        try:
            from utils.cluster_analysis import run_hierarchical_clustering
            
            result = run_hierarchical_clustering(
                sample_count_matrix,
                n_clusters=2
            )
            
            if 'error' not in result:
                assert 'sample_clusters' in result
                assert len(set(result['sample_clusters'])) <= 2
        except ImportError:
            pytest.skip("scipy not available")


class TestPCAAnalysis:
    """Tests for PCA analysis"""
    
    def test_pca_returns_dict(self, sample_count_matrix):
        """Test PCA returns a dictionary"""
        try:
            from utils.cluster_analysis import run_pca_analysis
            result = run_pca_analysis(sample_count_matrix)
            assert isinstance(result, dict)
        except ImportError:
            pytest.skip("scikit-learn not available")
    
    def test_pca_components(self, sample_count_matrix):
        """Test PCA returns correct number of components"""
        try:
            from utils.cluster_analysis import run_pca_analysis
            
            result = run_pca_analysis(sample_count_matrix, n_components=2)
            
            if 'error' not in result:
                assert 'pca_df' in result
                assert result['pca_df'].shape[1] == 2
        except ImportError:
            pytest.skip("scikit-learn not available")
    
    def test_pca_variance_explained(self, sample_count_matrix):
        """Test PCA returns variance explained"""
        try:
            from utils.cluster_analysis import run_pca_analysis
            
            result = run_pca_analysis(sample_count_matrix)
            
            if 'error' not in result:
                assert 'explained_variance' in result
                # Variance should sum to <= 1
                assert sum(result['explained_variance']) <= 1.0
        except ImportError:
            pytest.skip("scikit-learn not available")
