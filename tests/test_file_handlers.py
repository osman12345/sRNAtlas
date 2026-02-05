"""
Tests for file handler utilities
"""
import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))


class TestCountMatrix:
    """Tests for count matrix handling"""
    
    def test_count_matrix_shape(self, sample_count_matrix):
        """Test count matrix has expected shape"""
        assert sample_count_matrix.shape == (100, 4)
    
    def test_count_matrix_columns(self, sample_count_matrix):
        """Test count matrix has expected columns"""
        expected = ['Sample1', 'Sample2', 'Sample3', 'Sample4']
        assert list(sample_count_matrix.columns) == expected
    
    def test_count_matrix_values_non_negative(self, sample_count_matrix):
        """Test all count values are non-negative"""
        assert (sample_count_matrix >= 0).all().all()
    
    def test_count_matrix_index(self, sample_count_matrix):
        """Test count matrix has proper index"""
        assert all(idx.startswith('miRNA_') for idx in sample_count_matrix.index)


class TestMetadata:
    """Tests for metadata handling"""
    
    def test_metadata_columns(self, sample_metadata):
        """Test metadata has required columns"""
        required = ['sample', 'condition']
        assert all(col in sample_metadata.columns for col in required)
    
    def test_metadata_sample_count(self, sample_metadata):
        """Test metadata has expected samples"""
        assert len(sample_metadata) == 4
    
    def test_metadata_conditions(self, sample_metadata):
        """Test metadata has expected conditions"""
        conditions = set(sample_metadata['condition'])
        assert conditions == {'control', 'treatment'}


class TestFastqParsing:
    """Tests for FASTQ file parsing"""
    
    def test_fastq_file_exists(self, temp_fastq_file):
        """Test temporary FASTQ file was created"""
        assert temp_fastq_file.exists()
    
    def test_fastq_content_format(self, temp_fastq_file):
        """Test FASTQ file has correct format"""
        with open(temp_fastq_file) as f:
            lines = f.readlines()
        
        # FASTQ has 4 lines per read
        assert len(lines) % 4 == 0
        
        # First line should start with @
        assert lines[0].startswith('@')
        
        # Third line should be +
        assert lines[2].strip() == '+'
    
    def test_fastq_sequence_length(self, temp_fastq_file):
        """Test sequences are in expected length range"""
        with open(temp_fastq_file) as f:
            lines = f.readlines()
        
        for i in range(1, len(lines), 4):
            seq = lines[i].strip()
            assert 18 <= len(seq) <= 30, f"Sequence length {len(seq)} out of range"


class TestFastaParsing:
    """Tests for FASTA file parsing"""
    
    def test_fasta_file_exists(self, temp_fasta_file):
        """Test temporary FASTA file was created"""
        assert temp_fasta_file.exists()
    
    def test_fasta_content_format(self, temp_fasta_file):
        """Test FASTA file has correct format"""
        with open(temp_fasta_file) as f:
            content = f.read()
        
        # Should have header lines starting with >
        assert '>' in content
        
        # Count sequences
        headers = [l for l in content.split('\n') if l.startswith('>')]
        assert len(headers) == 3
    
    def test_fasta_mirna_names(self, temp_fasta_file):
        """Test FASTA contains expected miRNA names"""
        with open(temp_fasta_file) as f:
            content = f.read()
        
        assert 'hsa-miR-21-5p' in content
        assert 'hsa-let-7a-5p' in content
