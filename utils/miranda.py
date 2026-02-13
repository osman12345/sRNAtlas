"""
miRanda integration for animal miRNA target prediction
Provides wrapper functions for running miRanda and parsing results
"""
import subprocess
import tempfile
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd
import re
from dataclasses import dataclass


@dataclass
class MirandaHit:
    """A single miRanda prediction hit"""
    mirna_id: str
    target_id: str
    score: float
    energy: float
    query_start: int
    query_end: int
    ref_start: int
    ref_end: int
    alignment_length: int
    identity: float
    similarity: float
    query_seq: str
    match_seq: str
    ref_seq: str


def check_miranda_installed() -> Tuple[bool, str]:
    """
    Check if miRanda is installed and accessible

    Returns:
        Tuple of (is_installed, version_or_error)
    """
    try:
        result = subprocess.run(
            ['miranda', '-h'],
            capture_output=True,
            text=True,
            timeout=10
        )
        # miRanda outputs to stderr on -h
        output = result.stderr or result.stdout
        if 'miranda' in output.lower() or result.returncode == 0:
            # Try to extract version
            version_match = re.search(r'v?(\d+\.?\d*)', output)
            version = version_match.group(0) if version_match else "unknown"
            return True, f"miRanda {version}"
        return False, "miRanda not found"
    except FileNotFoundError:
        return False, "miRanda executable not found in PATH"
    except subprocess.TimeoutExpired:
        return False, "miRanda check timed out"
    except Exception as e:
        return False, str(e)


def get_miranda_install_instructions() -> str:
    """Get installation instructions for miRanda"""
    return """
**Installing miRanda:**

**Option 1: Conda (Recommended)**
```bash
conda install -c bioconda miranda
```

**Option 2: From Source**
```bash
# Download from http://www.microrna.org/
wget http://cbio.mskcc.org/microrna_data/miRanda-aug2010.tar.gz
tar xzf miRanda-aug2010.tar.gz
cd miRanda-3.3a
./configure
make
sudo make install
```

**Option 3: Ubuntu/Debian**
```bash
sudo apt-get install miranda
```
"""


def run_miranda(
    mirna_fasta: str,
    target_fasta: str,
    score_threshold: float = 140.0,
    energy_threshold: float = -20.0,
    gap_open: float = -9.0,
    gap_extend: float = -4.0,
    scale: float = 4.0,
    strict: bool = False,
    quiet: bool = True,
    output_dir: Optional[Path] = None
) -> Dict:
    """
    Run miRanda target prediction

    Args:
        mirna_fasta: miRNA sequences in FASTA format (string)
        target_fasta: Target 3' UTR sequences in FASTA format (string)
        score_threshold: Minimum alignment score (default: 140)
        energy_threshold: Maximum free energy (default: -20 kcal/mol)
        gap_open: Gap open penalty (default: -9)
        gap_extend: Gap extend penalty (default: -4)
        scale: Scaling parameter (default: 4)
        strict: Use strict seed pairing (default: False)
        quiet: Suppress verbose output (default: True)
        output_dir: Directory for output files (uses temp if None)

    Returns:
        Dictionary with status, results DataFrame, and statistics
    """
    # Check if miRanda is installed
    is_installed, version = check_miranda_installed()
    if not is_installed:
        return {
            'status': 'error',
            'error': f"miRanda not installed: {version}",
            'install_instructions': get_miranda_install_instructions()
        }

    # Create temp directory if needed
    if output_dir is None:
        temp_dir = Path(tempfile.mkdtemp())
        cleanup_temp = True
    else:
        temp_dir = output_dir
        temp_dir.mkdir(parents=True, exist_ok=True)
        cleanup_temp = False

    try:
        # Write input files
        mirna_file = temp_dir / "mirnas.fa"
        target_file = temp_dir / "targets.fa"
        output_file = temp_dir / "miranda_output.txt"

        with open(mirna_file, 'w') as f:
            f.write(mirna_fasta)

        with open(target_file, 'w') as f:
            f.write(target_fasta)

        # Build miRanda command
        cmd = [
            'miranda',
            str(mirna_file),
            str(target_file),
            '-sc', str(score_threshold),
            '-en', str(energy_threshold),
            '-go', str(gap_open),
            '-ge', str(gap_extend),
            '-scale', str(scale),
            '-out', str(output_file)
        ]

        if strict:
            cmd.append('-strict')

        if quiet:
            cmd.append('-quiet')

        # Run miRanda
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600  # 1 hour timeout for large runs
        )

        if result.returncode != 0:
            return {
                'status': 'error',
                'error': f"miRanda failed: {result.stderr}",
                'command': ' '.join(cmd)
            }

        # Parse output
        hits = parse_miranda_output(output_file)

        # Convert to DataFrame
        if hits:
            results_df = pd.DataFrame([
                {
                    'miRNA': h.mirna_id,
                    'Target': h.target_id,
                    'Score': h.score,
                    'Energy': h.energy,
                    'miRNA_Start': h.query_start,
                    'miRNA_End': h.query_end,
                    'Target_Start': h.ref_start,
                    'Target_End': h.ref_end,
                    'Alignment_Length': h.alignment_length,
                    'Identity': h.identity,
                    'Similarity': h.similarity,
                    'miRNA_Alignment': h.query_seq,
                    'Match_Pattern': h.match_seq,
                    'Target_Alignment': h.ref_seq
                }
                for h in hits
            ])

            # Sort by score (descending) and energy (ascending/more negative)
            results_df = results_df.sort_values(
                ['Score', 'Energy'],
                ascending=[False, True]
            ).reset_index(drop=True)
        else:
            results_df = pd.DataFrame()

        # Calculate statistics
        stats = {
            'total_predictions': len(results_df),
            'unique_mirnas': results_df['miRNA'].nunique() if len(results_df) > 0 else 0,
            'unique_targets': results_df['Target'].nunique() if len(results_df) > 0 else 0,
            'mean_score': results_df['Score'].mean() if len(results_df) > 0 else 0,
            'mean_energy': results_df['Energy'].mean() if len(results_df) > 0 else 0,
            'version': version
        }

        return {
            'status': 'success',
            'results': results_df,
            'stats': stats,
            'parameters': {
                'score_threshold': score_threshold,
                'energy_threshold': energy_threshold,
                'strict': strict
            }
        }

    except subprocess.TimeoutExpired:
        return {
            'status': 'error',
            'error': 'miRanda prediction timed out (>1 hour)'
        }
    except Exception as e:
        return {
            'status': 'error',
            'error': str(e)
        }
    finally:
        if cleanup_temp and temp_dir.exists():
            shutil.rmtree(temp_dir, ignore_errors=True)


def parse_miranda_output(output_file: Path) -> List[MirandaHit]:
    """
    Parse miRanda output file

    Args:
        output_file: Path to miRanda output file

    Returns:
        List of MirandaHit objects
    """
    hits = []

    if not output_file.exists():
        return hits

    with open(output_file, 'r') as f:
        content = f.read()

    # miRanda output format:
    # >miRNA\tTarget\tScore\tEnergy\tQuery_Start\tQuery_End\tRef_Start\tRef_End\tAl_Len\tIdentity\tSimilarity\tQuery_Seq\tMatch_Seq\tRef_Seq
    # or blocks starting with ">>" for detailed alignments

    # Pattern for tab-delimited summary lines
    summary_pattern = re.compile(
        r'^>(\S+)\t(\S+)\t([\d.]+)\t([-\d.]+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t([\d.]+)%?\t([\d.]+)%?(?:\t(\S*)\t(\S*)\t(\S*))?',
        re.MULTILINE
    )

    for match in summary_pattern.finditer(content):
        try:
            hit = MirandaHit(
                mirna_id=match.group(1),
                target_id=match.group(2),
                score=float(match.group(3)),
                energy=float(match.group(4)),
                query_start=int(match.group(5)),
                query_end=int(match.group(6)),
                ref_start=int(match.group(7)),
                ref_end=int(match.group(8)),
                alignment_length=int(match.group(9)),
                identity=float(match.group(10)),
                similarity=float(match.group(11)),
                query_seq=match.group(12) or "",
                match_seq=match.group(13) or "",
                ref_seq=match.group(14) or ""
            )
            hits.append(hit)
        except (ValueError, IndexError):
            continue

    # If no summary lines found, try parsing detailed format
    if not hits:
        hits = parse_miranda_detailed_output(content)

    return hits


def parse_miranda_detailed_output(content: str) -> List[MirandaHit]:
    """
    Parse miRanda detailed alignment output format

    Args:
        content: Raw miRanda output content

    Returns:
        List of MirandaHit objects
    """
    hits = []

    # Pattern for alignment blocks
    block_pattern = re.compile(
        r'>>(\S+)\s+(\S+)\s+([\d.]+)\s+([-\d.]+).*?'
        r'Query:\s*(\d+)\s+(\S+)\s+(\d+).*?'
        r'(\|*\s*[:|]*\s*\|*).*?'
        r'Ref:\s*(\d+)\s+(\S+)\s+(\d+)',
        re.DOTALL
    )

    for match in block_pattern.finditer(content):
        try:
            query_seq = match.group(6)
            ref_seq = match.group(10)
            match_pattern = match.group(8).strip()

            hit = MirandaHit(
                mirna_id=match.group(1),
                target_id=match.group(2),
                score=float(match.group(3)),
                energy=float(match.group(4)),
                query_start=int(match.group(5)),
                query_end=int(match.group(7)),
                ref_start=int(match.group(9)),
                ref_end=int(match.group(11)),
                alignment_length=len(query_seq),
                identity=calculate_identity(query_seq, ref_seq),
                similarity=calculate_similarity(match_pattern),
                query_seq=query_seq,
                match_seq=match_pattern,
                ref_seq=ref_seq
            )
            hits.append(hit)
        except (ValueError, IndexError):
            continue

    return hits


def calculate_identity(seq1: str, seq2: str) -> float:
    """Calculate sequence identity percentage"""
    if not seq1 or not seq2:
        return 0.0
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return 100.0 * matches / max(len(seq1), len(seq2))


def calculate_similarity(match_pattern: str) -> float:
    """Calculate similarity from match pattern (| and :)"""
    if not match_pattern:
        return 0.0
    matches = match_pattern.count('|') + match_pattern.count(':')
    total = len(match_pattern.replace(' ', ''))
    return 100.0 * matches / total if total > 0 else 0.0


def filter_miranda_results(
    results: pd.DataFrame,
    min_score: float = 140.0,
    max_energy: float = -20.0,
    min_identity: float = 0.0,
    seed_strict: bool = False
) -> pd.DataFrame:
    """
    Filter miRanda results by various criteria

    Args:
        results: miRanda results DataFrame
        min_score: Minimum alignment score
        max_energy: Maximum free energy (more negative = stronger binding)
        min_identity: Minimum sequence identity percentage
        seed_strict: If True, filter for strict seed pairing (positions 2-8)

    Returns:
        Filtered DataFrame
    """
    if results.empty:
        return results

    filtered = results.copy()

    # Score filter
    if min_score > 0:
        filtered = filtered[filtered['Score'] >= min_score]

    # Energy filter (more negative = better)
    if max_energy < 0:
        filtered = filtered[filtered['Energy'] <= max_energy]

    # Identity filter
    if min_identity > 0:
        filtered = filtered[filtered['Identity'] >= min_identity]

    return filtered


def get_miranda_summary_stats(results: pd.DataFrame) -> Dict:
    """
    Generate summary statistics from miRanda results

    Args:
        results: miRanda results DataFrame

    Returns:
        Dictionary with summary statistics
    """
    if results.empty:
        return {
            'total_predictions': 0,
            'unique_mirnas': 0,
            'unique_targets': 0
        }

    # Targets per miRNA
    targets_per_mirna = results.groupby('miRNA')['Target'].nunique()

    # miRNAs per target
    mirnas_per_target = results.groupby('Target')['miRNA'].nunique()

    return {
        'total_predictions': len(results),
        'unique_mirnas': results['miRNA'].nunique(),
        'unique_targets': results['Target'].nunique(),
        'mean_targets_per_mirna': targets_per_mirna.mean(),
        'median_targets_per_mirna': targets_per_mirna.median(),
        'max_targets_per_mirna': targets_per_mirna.max(),
        'mean_mirnas_per_target': mirnas_per_target.mean(),
        'median_mirnas_per_target': mirnas_per_target.median(),
        'max_mirnas_per_target': mirnas_per_target.max(),
        'mean_score': results['Score'].mean(),
        'mean_energy': results['Energy'].mean(),
        'score_range': (results['Score'].min(), results['Score'].max()),
        'energy_range': (results['Energy'].min(), results['Energy'].max())
    }


def get_top_targets(
    results: pd.DataFrame,
    mirna_id: str,
    n: int = 10,
    sort_by: str = 'Score'
) -> pd.DataFrame:
    """
    Get top N targets for a specific miRNA

    Args:
        results: miRanda results DataFrame
        mirna_id: miRNA identifier
        n: Number of top targets to return
        sort_by: Column to sort by ('Score' or 'Energy')

    Returns:
        DataFrame with top targets
    """
    if results.empty:
        return pd.DataFrame()

    mirna_results = results[results['miRNA'] == mirna_id].copy()

    if sort_by == 'Energy':
        # More negative energy = better binding
        mirna_results = mirna_results.sort_values('Energy', ascending=True)
    else:
        # Higher score = better match
        mirna_results = mirna_results.sort_values('Score', ascending=False)

    return mirna_results.head(n)


def get_top_mirnas_for_target(
    results: pd.DataFrame,
    target_id: str,
    n: int = 10,
    sort_by: str = 'Score'
) -> pd.DataFrame:
    """
    Get top N miRNAs targeting a specific gene

    Args:
        results: miRanda results DataFrame
        target_id: Target gene identifier
        n: Number of top miRNAs to return
        sort_by: Column to sort by ('Score' or 'Energy')

    Returns:
        DataFrame with top miRNAs
    """
    if results.empty:
        return pd.DataFrame()

    target_results = results[results['Target'] == target_id].copy()

    if sort_by == 'Energy':
        target_results = target_results.sort_values('Energy', ascending=True)
    else:
        target_results = target_results.sort_values('Score', ascending=False)

    return target_results.head(n)
