"""
Pytest configuration and fixtures for sRNAtlas tests
"""
import pytest
import pandas as pd
import numpy as np
import tempfile
from pathlib import Path


@pytest.fixture
def sample_count_matrix():
    """Create a sample count matrix for testing"""
    np.random.seed(42)
    data = {
        'Sample1': np.random.randint(0, 1000, 100),
        'Sample2': np.random.randint(0, 1000, 100),
        'Sample3': np.random.randint(0, 1000, 100),
        'Sample4': np.random.randint(0, 1000, 100),
    }
    index = [f'miRNA_{i}' for i in range(100)]
    return pd.DataFrame(data, index=index)


@pytest.fixture
def sample_metadata():
    """Create sample metadata for testing"""
    return pd.DataFrame({
        'sample': ['Sample1', 'Sample2', 'Sample3', 'Sample4'],
        'condition': ['control', 'control', 'treatment', 'treatment'],
        'batch': ['1', '1', '2', '2']
    })


@pytest.fixture
def sample_fastq_content():
    """Create sample FASTQ content for testing"""
    return """@SEQ_1
TGAGGTAGTAGGTTGTATAGTT
+
IIIIIIIIIIIIIIIIIIIII
@SEQ_2
TAGCTTATCAGACTGATGTTGA
+
IIIIIIIIIIIIIIIIIIIII
@SEQ_3
TGAGGTAGTAGGTTGTATAGTT
+
IIIIIIIIIIIIIIIIIIIII
"""


@pytest.fixture
def temp_fastq_file(sample_fastq_content):
    """Create a temporary FASTQ file"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
        f.write(sample_fastq_content)
        return Path(f.name)


@pytest.fixture
def sample_fasta_content():
    """Create sample FASTA content for testing"""
    return """>hsa-miR-21-5p
TAGCTTATCAGACTGATGTTGA
>hsa-let-7a-5p
TGAGGTAGTAGGTTGTATAGTT
>hsa-miR-155-5p
TTAATGCTAATCGTGATAGGGGT
"""


@pytest.fixture
def temp_fasta_file(sample_fasta_content):
    """Create a temporary FASTA file"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(sample_fasta_content)
        return Path(f.name)
