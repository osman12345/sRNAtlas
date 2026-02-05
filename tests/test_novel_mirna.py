"""
Unit tests for Novel miRNA Discovery Module
"""
import pytest
import tempfile
import gzip
from pathlib import Path


class TestNovelMiRNADiscovery:
    """Test novel miRNA discovery functionality"""

    def test_sequence_counting_from_fastq(self, tmp_path):
        """Test reading and counting sequences from FASTQ"""
        # Create test FASTQ file
        fastq_content = """@read1
TGAGGTAGTAGGTTGTATAGTT
+
IIIIIIIIIIIIIIIIIIIIII
@read2
TGAGGTAGTAGGTTGTATAGTT
+
IIIIIIIIIIIIIIIIIIIIII
@read3
AAAGGTAGTAGGTTGTATAGTT
+
IIIIIIIIIIIIIIIIIIIIII
@read4
TGAGGTAGTAGGTTGTATAGTT
+
IIIIIIIIIIIIIIIIIIIIII
"""
        fastq_file = tmp_path / "test.fastq"
        fastq_file.write_text(fastq_content)

        # Count sequences (mimicking run_novel_discovery logic)
        seq_counts = {}
        min_len, max_len = 20, 24

        with open(fastq_file, 'r') as f:
            line_num = 0
            for line in f:
                line_num += 1
                if line_num % 4 == 2:  # Sequence line
                    seq = line.strip()
                    if min_len <= len(seq) <= max_len:
                        seq_counts[seq] = seq_counts.get(seq, 0) + 1

        assert len(seq_counts) == 2  # Two unique sequences
        assert seq_counts["TGAGGTAGTAGGTTGTATAGTT"] == 3  # Appears 3 times
        assert seq_counts["AAAGGTAGTAGGTTGTATAGTT"] == 1  # Appears 1 time

    def test_sequence_counting_from_gzip(self, tmp_path):
        """Test reading sequences from gzipped FASTQ"""
        fastq_content = """@read1
TGAGGTAGTAGGTTGTATAGTT
+
IIIIIIIIIIIIIIIIIIIIII
@read2
TGAGGTAGTAGGTTGTATAGTT
+
IIIIIIIIIIIIIIIIIIIIII
"""
        fastq_file = tmp_path / "test.fastq.gz"
        with gzip.open(fastq_file, 'wt') as f:
            f.write(fastq_content)

        seq_counts = {}
        min_len, max_len = 20, 24

        with gzip.open(fastq_file, 'rt') as f:
            line_num = 0
            for line in f:
                line_num += 1
                if line_num % 4 == 2:
                    seq = line.strip()
                    if min_len <= len(seq) <= max_len:
                        seq_counts[seq] = seq_counts.get(seq, 0) + 1

        assert seq_counts["TGAGGTAGTAGGTTGTATAGTT"] == 2

    def test_length_filtering(self):
        """Test that sequences are filtered by length"""
        sequences = [
            "ATCG",  # Too short (4 nt)
            "TGAGGTAGTAGGTTGTATAG",  # Valid (20 nt)
            "TGAGGTAGTAGGTTGTATAGTT",  # Valid (22 nt)
            "TGAGGTAGTAGGTTGTATAGTTAA",  # Valid (24 nt)
            "TGAGGTAGTAGGTTGTATAGTTAAGG",  # Too long (26 nt)
        ]

        min_len, max_len = 20, 24
        filtered = [s for s in sequences if min_len <= len(s) <= max_len]

        assert len(filtered) == 3
        assert "ATCG" not in filtered
        assert "TGAGGTAGTAGGTTGTATAGTTAAGG" not in filtered

    def test_u_bias_detection(self):
        """Test 5' U/T bias detection"""
        test_cases = [
            ("TGAGGTAGTAGGTTGTATAGTT", True),  # Starts with T
            ("UGAGGTAGTAGGTTGTATAGTT", True),  # Starts with U
            ("AGAGGTAGTAGGTTGTATAGTT", False),  # Starts with A
            ("CGAGGTAGTAGGTTGTATAGTT", False),  # Starts with C
            ("GGAGGTAGTAGGTTGTATAGTT", False),  # Starts with G
        ]

        for seq, expected_u_bias in test_cases:
            has_u_bias = seq[0] in ['T', 'U']
            assert has_u_bias == expected_u_bias, f"Failed for {seq}"

    def test_gc_content_calculation(self):
        """Test GC content calculation"""
        test_cases = [
            ("AAAAAAAAAA", 0.0),    # All A
            ("GGGGGGGGGG", 100.0),  # All G
            ("CCCCCCCCCC", 100.0),  # All C
            ("ATCGATCGAT", 40.0),   # 4 G/C out of 10
            ("GCGCGCGCGC", 100.0),  # All G/C
        ]

        for seq, expected_gc in test_cases:
            gc_content = (seq.count('G') + seq.count('C')) / len(seq) * 100
            assert gc_content == expected_gc, f"Failed for {seq}"

    def test_candidate_scoring(self):
        """Test candidate scoring logic"""
        # Simulate candidate creation
        def create_candidate(seq, count, prefer_u=True):
            candidate = {
                'sequence': seq,
                'length': len(seq),
                'count': count,
                'first_nt': seq[0],
                'u_bias': seq[0] in ['T', 'U'],
                'gc_content': (seq.count('G') + seq.count('C')) / len(seq) * 100
            }

            # Apply U preference scoring
            if prefer_u and not candidate['u_bias']:
                candidate['score'] = 0.5
            else:
                candidate['score'] = 1.0

            # Read count contribution
            candidate['score'] *= min(count / 1000, 1.0)

            return candidate

        # Test U bias scoring
        u_biased = create_candidate("TGAGGTAGTAGGTTGTATAGTT", 1000)
        non_u = create_candidate("AGAGGTAGTAGGTTGTATAGTT", 1000)

        assert u_biased['score'] > non_u['score']
        assert u_biased['score'] == 1.0
        assert non_u['score'] == 0.5

        # Test count contribution
        high_count = create_candidate("TGAGGTAGTAGGTTGTATAGTT", 2000)
        low_count = create_candidate("TGAGGTAGTAGGTTGTATAGTT", 100)

        assert high_count['score'] == 1.0  # Capped at 1.0
        assert low_count['score'] == 0.1   # 100/1000


class TestMiRNACharacteristics:
    """Test miRNA-like characteristic detection"""

    def test_typical_mirna_length(self):
        """Test typical miRNA length range (20-24 nt)"""
        typical_lengths = range(20, 25)

        for length in typical_lengths:
            seq = "T" + "A" * (length - 1)
            assert 20 <= len(seq) <= 24

    def test_mirna_vs_sirna_vs_pirna_length(self):
        """Test different small RNA length characteristics"""
        size_ranges = {
            'miRNA': (20, 24),
            'siRNA': (21, 24),
            'piRNA': (24, 32),
        }

        # 22-nt could be miRNA or siRNA
        assert size_ranges['miRNA'][0] <= 22 <= size_ranges['miRNA'][1]
        assert size_ranges['siRNA'][0] <= 22 <= size_ranges['siRNA'][1]

        # 28-nt is piRNA
        assert size_ranges['piRNA'][0] <= 28 <= size_ranges['piRNA'][1]
        assert not (size_ranges['miRNA'][0] <= 28 <= size_ranges['miRNA'][1])
