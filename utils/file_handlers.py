"""
File handling utilities for sRNAtlas
"""
import os
import gzip
import tempfile
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Union
import pandas as pd
from Bio import SeqIO
from loguru import logger


def save_uploaded_file(uploaded_file, destination_dir: Path) -> Path:
    """
    Save an uploaded file to the destination directory

    Args:
        uploaded_file: Streamlit uploaded file object
        destination_dir: Directory to save the file

    Returns:
        Path to the saved file
    """
    destination_dir.mkdir(parents=True, exist_ok=True)
    file_path = destination_dir / uploaded_file.name

    with open(file_path, 'wb') as f:
        f.write(uploaded_file.getbuffer())

    logger.info(f"Saved uploaded file: {file_path}")
    return file_path


def read_fastq_stats(fastq_file: Union[str, Path]) -> Dict:
    """
    Read basic statistics from a FASTQ file

    Args:
        fastq_file: Path to FASTQ file (can be gzipped)

    Returns:
        Dictionary with read statistics
    """
    fastq_file = Path(fastq_file)

    # Determine if file is gzipped
    if str(fastq_file).endswith('.gz'):
        handle = gzip.open(fastq_file, 'rt')
    else:
        handle = open(fastq_file, 'r')

    read_count = 0
    total_length = 0
    lengths = []
    quality_scores = []

    try:
        for record in SeqIO.parse(handle, 'fastq'):
            read_count += 1
            length = len(record.seq)
            total_length += length
            lengths.append(length)

            # Calculate mean quality for this read
            if record.letter_annotations.get('phred_quality'):
                mean_qual = sum(record.letter_annotations['phred_quality']) / length
                quality_scores.append(mean_qual)

            # Only process first 100k reads for speed
            if read_count >= 100000:
                break

    finally:
        handle.close()

    stats = {
        'read_count': read_count,
        'total_bases': total_length,
        'mean_length': total_length / read_count if read_count > 0 else 0,
        'min_length': min(lengths) if lengths else 0,
        'max_length': max(lengths) if lengths else 0,
        'mean_quality': sum(quality_scores) / len(quality_scores) if quality_scores else 0,
        'length_distribution': lengths,
        'sampled': read_count >= 100000
    }

    return stats


def read_count_matrix(file_path: Union[str, Path], **kwargs) -> pd.DataFrame:
    """
    Read a count matrix from CSV/TSV file

    Args:
        file_path: Path to count matrix file
        **kwargs: Additional arguments for pandas read_csv

    Returns:
        DataFrame with count matrix
    """
    file_path = Path(file_path)

    # Detect delimiter
    with open(file_path, 'r') as f:
        first_line = f.readline()
        if '\t' in first_line:
            sep = '\t'
        else:
            sep = ','

    df = pd.read_csv(file_path, sep=sep, index_col=0, **kwargs)
    logger.info(f"Loaded count matrix: {df.shape[0]} features x {df.shape[1]} samples")

    return df


def read_metadata(file_path: Union[str, Path]) -> pd.DataFrame:
    """
    Read sample metadata from CSV file

    Args:
        file_path: Path to metadata CSV

    Returns:
        DataFrame with sample metadata
    """
    file_path = Path(file_path)

    df = pd.read_csv(file_path)

    # Standardize column names
    df.columns = df.columns.str.strip()

    # Check for required columns
    required_cols = ['SampleID']
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        # Try to infer SampleID
        if 'sample' in df.columns.str.lower().tolist():
            sample_col = df.columns[df.columns.str.lower() == 'sample'][0]
            df = df.rename(columns={sample_col: 'SampleID'})
        else:
            raise ValueError(f"Missing required columns: {missing}")

    logger.info(f"Loaded metadata: {len(df)} samples")
    return df


def parse_fasta_annotations(fasta_file: Union[str, Path]) -> pd.DataFrame:
    """
    Parse annotations from FASTA headers (supports miRBase and RNAcentral formats)

    Supported formats:
    - miRBase: >mtr-miR162_MIMAT0001640_Medicago_truncatula_miR162
    - RNAcentral: >mtr-miR169l-5p_URS00000036DB_3880

    Args:
        fasta_file: Path to FASTA file

    Returns:
        DataFrame with annotations
    """
    fasta_file = Path(fasta_file)

    annotations = []

    # Handle gzipped files
    if str(fasta_file).endswith('.gz'):
        handle = gzip.open(fasta_file, 'rt')
    else:
        handle = open(fasta_file, 'r')

    try:
        for record in SeqIO.parse(handle, 'fasta'):
            header = record.id
            full_description = record.description  # Full header line

            # Try to detect format and parse accordingly
            annotation = _parse_header(header, full_description, len(record.seq))
            annotations.append(annotation)

    finally:
        handle.close()

    df = pd.DataFrame(annotations)
    logger.info(f"Parsed {len(df)} annotations from FASTA")

    return df


def _parse_header(header: str, full_description: str, seq_length: int) -> Dict:
    """
    Parse a single FASTA header into annotation fields

    Supports multiple formats:
    - miRBase mature: mtr-miR162_MIMAT0001640_Medicago_truncatula_miR162
    - miRBase stem-loop: mtr-MIR162_MI0001738_Medicago_truncatula_miR162_stem-loop
    - RNAcentral: URS00000045E3_3880 Medicago truncatula tRNA-Lys (TTT)
    - Generic: any_id description text
    """
    import re

    # Check for miRBase mature miRNA format: has MIMAT accession
    mimat_match = re.search(r'MIMAT\d+', header)
    if mimat_match:
        mimat_id = mimat_match.group()
        mimat_idx = header.find(mimat_id)
        mirna_name = header[:mimat_idx].rstrip('_')

        after_mimat = header[mimat_idx + len(mimat_id):].lstrip('_')
        after_parts = after_mimat.split('_')

        if len(after_parts) >= 3:
            organism = f"{after_parts[0]}_{after_parts[1]}"
            mirna_id = '_'.join(after_parts[2:])
        elif len(after_parts) == 2:
            organism = after_parts[0]
            mirna_id = after_parts[1]
        else:
            organism = after_parts[0] if after_parts else "unknown"
            mirna_id = mirna_name

        return {
            'full_id': header,
            'miRNA_name': mirna_name,
            'miRBase_accession': mimat_id,
            'organism': organism,
            'miRNA_id': mirna_id,
            'description': mirna_name,
            'RNA_category': 'miRNA',
            'RNA_subtype': 'mature_miRNA',
            'sequence_length': seq_length,
            'source': 'miRBase'
        }

    # Check for miRBase stem-loop/precursor format: has MI accession (not MIMAT)
    mi_match = re.search(r'MI\d{7}', header)  # MI followed by exactly 7 digits
    if mi_match:
        mi_id = mi_match.group()
        mi_idx = header.find(mi_id)
        mir_name = header[:mi_idx].rstrip('_')

        after_mi = header[mi_idx + len(mi_id):].lstrip('_')
        after_parts = after_mi.split('_')

        if len(after_parts) >= 2:
            organism = f"{after_parts[0]}_{after_parts[1]}"
            extra_info = '_'.join(after_parts[2:]) if len(after_parts) > 2 else ''
        else:
            organism = after_parts[0] if after_parts else "unknown"
            extra_info = ''

        return {
            'full_id': header,
            'miRNA_name': mir_name,
            'miRBase_accession': mi_id,
            'organism': organism,
            'description': mir_name,
            'RNA_category': 'miRNA',
            'RNA_subtype': 'stem-loop' if 'stem-loop' in header else 'precursor',
            'sequence_length': seq_length,
            'source': 'miRBase'
        }

    # Check for RNAcentral format: URS*_taxid description
    urs_match = re.search(r'URS[0-9A-F]+', header)
    if urs_match:
        urs_id = urs_match.group()

        # Split header to get taxid
        parts = header.split('_')
        if len(parts) >= 2:
            taxid = parts[1].split()[0]  # Get taxid before any space
        else:
            taxid = "unknown"

        # Use full_description for categorization (contains RNA type info)
        rna_category, rna_subtype = categorize_rna(full_description)

        return {
            'full_id': header,
            'URS_ID': urs_id,
            'description': full_description,
            'RNA_category': rna_category,
            'RNA_subtype': rna_subtype,
            'sequence_length': seq_length,
            'taxid': taxid,
            'source': 'RNAcentral'
        }

    # Generic format - try to extract what we can
    rna_category, rna_subtype = categorize_rna(full_description)

    return {
        'full_id': header,
        'description': full_description if full_description != header else header,
        'RNA_category': rna_category,
        'RNA_subtype': rna_subtype,
        'sequence_length': seq_length,
        'source': 'custom'
    }


def categorize_rna(description: str) -> Tuple[str, str]:
    """
    Categorize RNA type based on description

    Args:
        description: RNA description string

    Returns:
        Tuple of (category, subtype)
    """
    if pd.isna(description) or description == '':
        return 'unannotated', 'unannotated'

    desc_lower = str(description).lower()

    # miRNA
    if any(p in desc_lower for p in ['mir-', '-mir', 'mirna', 'microrna']):
        return 'miRNA', description

    # lncRNA
    if any(p in desc_lower for p in ['long_non-coding', 'lncrna', 'lincrna']):
        return 'lncRNA', 'lncRNA'

    # tRNA
    if 'trna' in desc_lower or 'transfer' in desc_lower:
        return 'tRNA', 'tRNA'

    # rRNA
    if 'rrna' in desc_lower or 'ribosomal' in desc_lower:
        if '18s' in desc_lower:
            return 'rRNA', '18S rRNA'
        elif '28s' in desc_lower or '25s' in desc_lower:
            return 'rRNA', '28S rRNA'
        elif '5.8s' in desc_lower:
            return 'rRNA', '5.8S rRNA'
        elif '5s' in desc_lower:
            return 'rRNA', '5S rRNA'
        return 'rRNA', 'rRNA'

    # snoRNA
    if 'snorna' in desc_lower or 'small_nucleolar' in desc_lower:
        return 'snoRNA', 'snoRNA'

    # snRNA
    if 'snrna' in desc_lower or 'small_nuclear' in desc_lower:
        return 'snRNA', 'snRNA'

    # siRNA
    if 'sirna' in desc_lower:
        return 'siRNA', 'siRNA'

    # piRNA
    if 'pirna' in desc_lower or 'piwi' in desc_lower:
        return 'piRNA', 'piRNA'

    return 'other_ncRNA', description


def validate_bam_file(bam_file: Union[str, Path]) -> Dict:
    """
    Validate a BAM file and return basic statistics

    Args:
        bam_file: Path to BAM file

    Returns:
        Dictionary with validation results
    """
    import pysam

    bam_file = Path(bam_file)
    result = {
        'valid': False,
        'indexed': False,
        'stats': None,
        'error': None
    }

    try:
        # Check if file exists
        if not bam_file.exists():
            result['error'] = f"File not found: {bam_file}"
            return result

        # Check if index exists
        index_file = Path(str(bam_file) + '.bai')
        result['indexed'] = index_file.exists()

        # Try to open and get stats
        with pysam.AlignmentFile(str(bam_file), 'rb') as bam:
            stats = {
                'mapped_reads': bam.mapped,
                'unmapped_reads': bam.unmapped,
                'references': len(bam.references),
            }
            result['stats'] = stats
            result['valid'] = True

    except Exception as e:
        result['error'] = str(e)

    return result


def create_project_structure(project_dir: Union[str, Path], project_name: str) -> Dict[str, Path]:
    """
    Create the directory structure for a new project

    Args:
        project_dir: Base directory for projects
        project_name: Name of the project

    Returns:
        Dictionary mapping directory names to paths
    """
    project_dir = Path(project_dir)
    project_path = project_dir / project_name

    directories = {
        'root': project_path,
        'data': project_path / 'data',
        'fastq': project_path / 'data' / 'fastq',
        'bam': project_path / 'data' / 'bam',
        'reference': project_path / 'data' / 'reference',
        'results': project_path / 'results',
        'qc': project_path / 'results' / 'qc',
        'alignment': project_path / 'results' / 'alignment',
        'counts': project_path / 'results' / 'counts',
        'de_analysis': project_path / 'results' / 'de_analysis',
        'enrichment': project_path / 'results' / 'enrichment',
        'reports': project_path / 'results' / 'reports',
        'logs': project_path / 'logs',
    }

    for name, path in directories.items():
        path.mkdir(parents=True, exist_ok=True)
        logger.debug(f"Created directory: {path}")

    logger.info(f"Created project structure for: {project_name}")
    return directories
