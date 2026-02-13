"""
Error handling utilities for sRNAtlas
Provides actionable error messages with diagnostic commands
"""
import subprocess
import shutil
from pathlib import Path
from typing import Optional, Dict, List, Tuple
from dataclasses import dataclass


@dataclass
class DiagnosticError:
    """Structured error with diagnostic information"""
    message: str
    category: str
    diagnostic_commands: List[str]
    suggested_fix: str
    severity: str = "error"  # error, warning, info


class ErrorClassifier:
    """Classifies errors and provides actionable diagnostics"""

    @staticmethod
    def classify_bam_error(error_msg: str, bam_file: Optional[Path] = None) -> DiagnosticError:
        """Classify BAM file errors and provide diagnostics"""
        error_lower = error_msg.lower()

        if 'bgzf' in error_lower or 'truncated' in error_lower:
            return DiagnosticError(
                message="BAM file is corrupted or truncated",
                category="bam_corruption",
                diagnostic_commands=[
                    f"samtools quickcheck {bam_file}" if bam_file else "samtools quickcheck <your_bam_file>",
                    f"samtools view -H {bam_file} | head" if bam_file else "samtools view -H <your_bam_file> | head"
                ],
                suggested_fix="Delete the corrupted BAM file and re-run alignment. The file may have been interrupted during creation.",
                severity="error"
            )

        elif 'index' in error_lower or 'bai' in error_lower:
            return DiagnosticError(
                message="BAM index file (.bai) is missing or corrupted",
                category="bam_index_missing",
                diagnostic_commands=[
                    f"ls -la {bam_file}.bai" if bam_file else "ls -la <your_bam_file>.bai",
                    f"samtools index {bam_file}" if bam_file else "samtools index <your_bam_file>"
                ],
                suggested_fix="Run 'samtools index' to create the index, or re-run alignment which creates the index automatically.",
                severity="error"
            )

        elif 'not found' in error_lower or 'no such file' in error_lower:
            return DiagnosticError(
                message="BAM file not found",
                category="file_not_found",
                diagnostic_commands=[
                    f"ls -la {bam_file.parent if bam_file else '.'}" if bam_file else "ls -la <directory>",
                ],
                suggested_fix="Check that the file path is correct and the file exists. You may need to re-run alignment.",
                severity="error"
            )

        else:
            return DiagnosticError(
                message=f"BAM file error: {error_msg}",
                category="bam_unknown",
                diagnostic_commands=[
                    f"samtools quickcheck {bam_file}" if bam_file else "samtools quickcheck <your_bam_file>",
                    f"samtools flagstat {bam_file}" if bam_file else "samtools flagstat <your_bam_file>"
                ],
                suggested_fix="Try re-running alignment to regenerate the BAM file.",
                severity="error"
            )

    @staticmethod
    def classify_alignment_error(error_msg: str, index_prefix: Optional[str] = None) -> DiagnosticError:
        """Classify alignment errors and provide diagnostics"""
        error_lower = error_msg.lower()

        if 'could not locate' in error_lower or 'index' in error_lower:
            return DiagnosticError(
                message="Bowtie index not found",
                category="index_missing",
                diagnostic_commands=[
                    f"ls -la {index_prefix}*" if index_prefix else "ls -la <index_prefix>*",
                    "bowtie-inspect --summary <index_prefix>"
                ],
                suggested_fix="Build the Bowtie index first using 'bowtie-build reference.fa index_prefix', or select a valid reference in the Databases module.",
                severity="error"
            )

        elif '0% alignment' in error_lower or 'no alignments' in error_lower:
            return DiagnosticError(
                message="No reads aligned to reference",
                category="zero_alignment",
                diagnostic_commands=[
                    "head -4 <fastq_file>  # Check read format",
                    "bowtie-inspect --summary <index>  # Check reference"
                ],
                suggested_fix="Possible causes: (1) Wrong reference organism, (2) Adapters not trimmed, (3) Reference index corrupted. Try trimming adapters first or verify the reference.",
                severity="warning"
            )

        else:
            return DiagnosticError(
                message=f"Alignment error: {error_msg}",
                category="alignment_unknown",
                diagnostic_commands=[
                    "bowtie --version",
                    "samtools --version"
                ],
                suggested_fix="Check that Bowtie and Samtools are properly installed and the input files are valid.",
                severity="error"
            )

    @staticmethod
    def classify_de_error(error_msg: str, count_matrix=None, metadata=None) -> DiagnosticError:
        """Classify differential expression errors"""
        error_lower = error_msg.lower()

        if 'sample' in error_lower and ('mismatch' in error_lower or 'not found' in error_lower or 'differ' in error_lower):
            return DiagnosticError(
                message="Sample names in count matrix don't match metadata",
                category="sample_mismatch",
                diagnostic_commands=[
                    "# Check count matrix columns vs metadata sample IDs",
                    "# Count matrix samples: [list column names]",
                    "# Metadata samples: [list SampleID column]"
                ],
                suggested_fix="Ensure sample names in count matrix columns exactly match the SampleID column in metadata. Check for trailing spaces or case differences.",
                severity="error"
            )

        elif 'singular' in error_lower or 'rank deficient' in error_lower:
            return DiagnosticError(
                message="Design matrix is singular (not enough samples or confounded design)",
                category="design_singular",
                diagnostic_commands=[],
                suggested_fix="You need at least 2 samples per group for DE analysis. Check that your design factor has sufficient replicates in each condition.",
                severity="error"
            )

        elif 'negative' in error_lower and 'count' in error_lower:
            return DiagnosticError(
                message="Count matrix contains negative values",
                category="negative_counts",
                diagnostic_commands=[],
                suggested_fix="Count matrices must contain only non-negative integers. Check your count matrix for errors.",
                severity="error"
            )

        else:
            return DiagnosticError(
                message=f"DE analysis error: {error_msg}",
                category="de_unknown",
                diagnostic_commands=[
                    "# Check that pyDESeq2 is installed:",
                    "pip show pydeseq2"
                ],
                suggested_fix="Check that your count matrix and metadata are properly formatted.",
                severity="error"
            )

    @staticmethod
    def classify_tool_error(tool_name: str) -> DiagnosticError:
        """Classify missing tool errors"""
        install_commands = {
            'bowtie': [
                "conda install -c bioconda bowtie",
                "# Or on Ubuntu/Debian:",
                "sudo apt-get install bowtie"
            ],
            'samtools': [
                "conda install -c bioconda samtools",
                "# Or on Ubuntu/Debian:",
                "sudo apt-get install samtools"
            ],
            'cutadapt': [
                "pip install cutadapt",
                "# Or with conda:",
                "conda install -c bioconda cutadapt"
            ],
            'fastqc': [
                "conda install -c bioconda fastqc",
                "# Or on Ubuntu/Debian:",
                "sudo apt-get install fastqc"
            ]
        }

        commands = install_commands.get(tool_name.lower(), [f"# Install {tool_name} using your package manager"])

        return DiagnosticError(
            message=f"{tool_name} is not installed or not in PATH",
            category="tool_missing",
            diagnostic_commands=[
                f"which {tool_name}",
                f"{tool_name} --version"
            ],
            suggested_fix=f"Install {tool_name} using one of these methods:\n" + "\n".join(commands),
            severity="error"
        )


def check_required_tools() -> Dict[str, bool]:
    """Check if required external tools are available"""
    tools = ['bowtie', 'samtools', 'cutadapt']
    results = {}

    for tool in tools:
        results[tool] = shutil.which(tool) is not None

    return results


def get_tool_versions() -> Dict[str, str]:
    """Get versions of installed tools"""
    versions = {}

    tool_commands = {
        'bowtie': ['bowtie', '--version'],
        'samtools': ['samtools', '--version'],
        'cutadapt': ['cutadapt', '--version'],
        'python': ['python', '--version']
    }

    for tool, cmd in tool_commands.items():
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                # Extract first line which usually contains version
                output = result.stdout.strip() or result.stderr.strip()
                versions[tool] = output.split('\n')[0]
            else:
                versions[tool] = "not found"
        except (subprocess.TimeoutExpired, FileNotFoundError):
            versions[tool] = "not found"

    return versions


def format_error_for_streamlit(error: DiagnosticError) -> str:
    """Format a DiagnosticError for display in Streamlit"""
    output = f"**{error.message}**\n\n"

    if error.suggested_fix:
        output += f"**Suggested fix:** {error.suggested_fix}\n\n"

    if error.diagnostic_commands:
        output += "**Diagnostic commands:**\n```bash\n"
        output += "\n".join(error.diagnostic_commands)
        output += "\n```"

    return output


def validate_sample_names(count_matrix, metadata, sample_col: str = 'SampleID') -> Tuple[bool, str]:
    """
    Validate that sample names match between count matrix and metadata

    Returns:
        Tuple of (is_valid, error_message)
    """
    import pandas as pd

    if count_matrix is None or metadata is None:
        return False, "Count matrix or metadata is missing"

    # Get sample names from count matrix (columns)
    count_samples = set(count_matrix.columns)

    # Get sample names from metadata
    if sample_col not in metadata.columns:
        # Try to find a matching column
        possible_cols = [c for c in metadata.columns if 'sample' in c.lower()]
        if possible_cols:
            sample_col = possible_cols[0]
        else:
            return False, f"No '{sample_col}' column found in metadata. Available columns: {list(metadata.columns)}"

    metadata_samples = set(metadata[sample_col].astype(str))

    # Check for mismatches
    in_counts_not_meta = count_samples - metadata_samples
    in_meta_not_counts = metadata_samples - count_samples

    if in_counts_not_meta or in_meta_not_counts:
        msg = "Sample name mismatch detected:\n"
        if in_counts_not_meta:
            msg += f"- In count matrix but not metadata: {list(in_counts_not_meta)[:5]}\n"
        if in_meta_not_counts:
            msg += f"- In metadata but not count matrix: {list(in_meta_not_counts)[:5]}\n"
        return False, msg

    return True, ""
