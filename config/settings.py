"""
Configuration settings for sRNAtlas
"""
import os
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Optional

# Base paths
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data"
JOBS_DIR = BASE_DIR / "jobs"
REFERENCES_DIR = DATA_DIR / "references"

# Create directories if they don't exist
for dir_path in [DATA_DIR, JOBS_DIR, REFERENCES_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)

@dataclass
class AlignmentSettings:
    """Bowtie alignment parameters (optimized for small RNA)"""
    threads: int = 8

    # Mismatch settings (-v mode is best for small RNA)
    mismatches: int = 1  # -v parameter: allow 0, 1, 2, or 3 mismatches

    # Multi-mapping settings
    max_alignments: int = 10  # -k parameter: report up to k alignments
    suppress_multi: int = 0  # -m parameter: suppress reads with >m alignments (0=disabled)
    report_all: bool = False  # -a parameter: report all alignments (overrides -k)
    best_mode: bool = True  # --best: report best alignments first
    strata: bool = True  # --strata: only report alignments in best stratum

    # Output settings
    no_unal: bool = True  # --un: output unaligned reads to file (False = suppress)
    sam_output: bool = True  # -S: output in SAM format

    # Size selection (for small RNA)
    min_length: int = 18
    max_length: int = 35

    # Seed settings
    seed_length: int = 18  # -l: seed length (shorter = more sensitive but slower)
    seed_mismatches: int = 1  # -n: max mismatches in seed (only if using -n mode)

    @property
    def bowtie_args(self) -> List[str]:
        """Generate Bowtie command line arguments"""
        args = [
            f"-v {self.mismatches}",  # Use -v mode for exact mismatch control
            f"-p {self.threads}",
        ]

        if self.report_all:
            args.append("-a")
        else:
            args.append(f"-k {self.max_alignments}")

        if self.suppress_multi > 0:
            args.append(f"-m {self.suppress_multi}")

        if self.best_mode:
            args.append("--best")

        if self.strata:
            args.append("--strata")

        if self.sam_output:
            args.append("-S")

        return args


@dataclass
class QCSettings:
    """Quality control parameters"""
    min_quality: int = 20
    min_length: int = 18
    max_length: int = 200  # Increased to capture snoRNA/snRNA
    adapter_sequence: str = "TGGAATTCTCGGGTGCCAAGG"  # Illumina small RNA adapter
    trim_quality: int = 20

    # Size distribution bins
    size_bins: List[int] = field(default_factory=lambda: list(range(15, 201)))

    # Expected size ranges for different RNA types
    size_ranges: Dict[str, tuple] = field(default_factory=lambda: {
        'miRNA': (18, 25),
        'siRNA': (20, 24),
        'piRNA': (24, 32),
        'tRF/tsRNA': (14, 40),
        'rsRF': (15, 40),
        'snoRNA': (60, 300),
        'snRNA': (100, 300),
        'Y_RNA': (80, 120),
    })

    # Comprehensive small RNA type definitions
    rna_types: Dict[str, Dict] = field(default_factory=lambda: {
        'miRNA': {
            'name': 'MicroRNA',
            'size': (18, 25),
            'description': 'Gene expression regulators, 5\' U bias',
            'color': '#1f77b4'
        },
        'siRNA': {
            'name': 'Small interfering RNA',
            'size': (20, 24),
            'description': 'RNAi pathway, perfect complementarity',
            'color': '#ff7f0e'
        },
        'piRNA': {
            'name': 'PIWI-interacting RNA',
            'size': (24, 32),
            'description': 'Transposon silencing, germline',
            'color': '#2ca02c'
        },
        'tRF': {
            'name': 'tRNA-derived fragment',
            'size': (14, 40),
            'description': 'Stress response, gene regulation',
            'color': '#d62728'
        },
        'rsRF': {
            'name': 'rRNA-derived fragment',
            'size': (15, 40),
            'description': 'rRNA processing products',
            'color': '#9467bd'
        },
        'snoRNA': {
            'name': 'Small nucleolar RNA',
            'size': (60, 300),
            'description': 'rRNA/snRNA modification guides',
            'color': '#8c564b'
        },
        'snRNA': {
            'name': 'Small nuclear RNA',
            'size': (100, 300),
            'description': 'Splicing machinery (U1, U2, U4, U5, U6)',
            'color': '#e377c2'
        },
        'Y_RNA': {
            'name': 'Y RNA',
            'size': (80, 120),
            'description': 'DNA replication, RNA quality control',
            'color': '#7f7f7f'
        },
        'vtRNA': {
            'name': 'Vault RNA',
            'size': (80, 150),
            'description': 'Vault ribonucleoprotein complex',
            'color': '#bcbd22'
        },
        '7SL': {
            'name': 'Signal recognition particle RNA',
            'size': (280, 320),
            'description': 'Protein targeting to ER',
            'color': '#17becf'
        },
    })

    # Plant-specific small RNA types
    plant_rna_types: Dict[str, Dict] = field(default_factory=lambda: {
        'tasiRNA': {
            'name': 'Trans-acting siRNA',
            'size': (21, 21),
            'description': 'Plant-specific, phased production',
            'color': '#aec7e8'
        },
        'phasiRNA': {
            'name': 'Phased siRNA',
            'size': (21, 24),
            'description': 'Phased pattern, reproductive tissues',
            'color': '#ffbb78'
        },
        'natsiRNA': {
            'name': 'Natural antisense siRNA',
            'size': (21, 24),
            'description': 'From overlapping transcripts',
            'color': '#98df8a'
        },
        'hc-siRNA': {
            'name': 'Heterochromatic siRNA',
            'size': (24, 24),
            'description': 'DNA methylation, transposon silencing',
            'color': '#ff9896'
        },
    })


@dataclass
class CountingSettings:
    """Read counting parameters"""
    min_mapq: int = 0  # Allow multi-mappers (important for small RNAs)
    count_mode: str = "union"  # How to count overlapping features
    strand_specific: bool = False

    # Filtering
    min_counts_per_sample: int = 1
    min_samples_detected: int = 2
    min_cpm: float = 0.5  # Minimum CPM threshold for filtering


@dataclass
class DESeqSettings:
    """DESeq2 parameters"""
    # Filtering
    min_cpm: float = 0.5
    min_samples: int = 2

    # FDR thresholds
    padj_threshold: float = 0.05
    lfc_threshold: float = 0.585  # log2(1.5)

    # SVA settings
    use_sva: bool = True
    max_surrogate_variables: int = 5
    sv_correlation_threshold: float = 0.3  # Warn if SV correlates with biology

    # Normalization
    normalization_method: str = "median_ratio"  # DESeq2 default


@dataclass
class EnrichmentSettings:
    """GO/Pathway enrichment parameters"""
    # GO settings
    go_ontologies: List[str] = field(default_factory=lambda: ["BP", "MF", "CC"])
    min_gene_set_size: int = 10
    max_gene_set_size: int = 500

    # FDR
    enrichment_padj: float = 0.05

    # KEGG
    use_kegg: bool = True

    # Visualization
    top_terms_to_show: int = 20


@dataclass
class RNACentralSettings:
    """RNAcentral database settings"""
    base_url: str = "https://rnacentral.org"
    api_url: str = "https://rnacentral.org/api/v1"

    # Supported organisms (taxid -> name)
    organisms: Dict[int, str] = field(default_factory=lambda: {
        3880: "Medicago truncatula",
        3702: "Arabidopsis thaliana",
        4530: "Oryza sativa",
        4577: "Zea mays",
        3847: "Glycine max",
        9606: "Homo sapiens",
        10090: "Mus musculus",
        7227: "Drosophila melanogaster",
        6239: "Caenorhabditis elegans",
    })

    # RNA type categories for classification
    rna_categories: Dict[str, List[str]] = field(default_factory=lambda: {
        'miRNA': ['miRNA', 'microRNA', 'mir-', 'let-', 'MIR'],
        'tRNA': ['tRNA', 'transfer RNA', 'Mt-tRNA', 'mt_tRNA'],
        'tRF': ['tRF', 'tsRNA', 'tRNA-derived', 'tRNA fragment'],
        'rRNA': ['rRNA', 'ribosomal RNA', '5S', '5.8S', '18S', '28S'],
        'snoRNA': ['snoRNA', 'small nucleolar RNA', 'SNORD', 'SNORA', 'scaRNA'],
        'snRNA': ['snRNA', 'small nuclear RNA', 'U1', 'U2', 'U4', 'U5', 'U6', 'U11', 'U12'],
        'siRNA': ['siRNA', 'small interfering RNA'],
        'piRNA': ['piRNA', 'piwi-interacting RNA', 'PIWI'],
        'Y_RNA': ['Y_RNA', 'Y RNA', 'RNY1', 'RNY3', 'RNY4', 'RNY5'],
        'vtRNA': ['vtRNA', 'vault RNA', 'VTRNA'],
        '7SL': ['7SL', 'SRP RNA', 'RN7SL'],
        '7SK': ['7SK', 'RN7SK'],
    })


@dataclass
class PathSettings:
    """Path configuration settings"""
    base_dir: Path = field(default_factory=lambda: BASE_DIR)
    data_dir: Path = field(default_factory=lambda: DATA_DIR)
    jobs_dir: Path = field(default_factory=lambda: JOBS_DIR)
    reference_dir: Path = field(default_factory=lambda: REFERENCES_DIR)
    output_dir: Path = field(default_factory=lambda: BASE_DIR / "output")

    def __post_init__(self):
        # Create directories if they don't exist
        for attr in ['data_dir', 'jobs_dir', 'reference_dir', 'output_dir']:
            path = getattr(self, attr)
            if isinstance(path, Path):
                path.mkdir(parents=True, exist_ok=True)


# Global configuration instance
class Config:
    """Main configuration class"""

    def __init__(self):
        self.alignment = AlignmentSettings()
        self.qc = QCSettings()
        self.counting = CountingSettings()
        self.deseq = DESeqSettings()
        self.enrichment = EnrichmentSettings()
        self.rnacentral = RNACentralSettings()
        self.paths = PathSettings()

        # App settings
        self.app_name = "sRNAtlas"
        self.app_tagline = "Comprehensive Small RNA-seq Analysis Platform"
        self.version = "BETA"
        self.debug = os.getenv("DEBUG", "false").lower() == "true"

        # Job settings
        self.max_concurrent_jobs = 3
        self.job_timeout_hours = 24

        # File size limits (in MB)
        self.max_upload_size_mb = 5000
        self.max_reference_size_mb = 2000

    def to_dict(self) -> dict:
        """Export configuration as dictionary"""
        return {
            'alignment': self.alignment.__dict__,
            'qc': self.qc.__dict__,
            'counting': self.counting.__dict__,
            'deseq': self.deseq.__dict__,
            'enrichment': self.enrichment.__dict__,
        }

    def validate(self) -> list:
        """
        Validate all configuration settings.

        Returns:
            List of validation error strings (empty if valid)
        """
        errors = []

        # Validate alignment settings
        if hasattr(self, 'alignment'):
            a = self.alignment
            if hasattr(a, 'mismatches') and not (0 <= a.mismatches <= 3):
                errors.append(f"alignment.mismatches must be 0-3, got {a.mismatches}")
            if hasattr(a, 'max_alignments') and a.max_alignments < 1:
                errors.append(f"alignment.max_alignments must be >= 1, got {a.max_alignments}")

        # Validate QC settings
        if hasattr(self, 'qc'):
            q = self.qc
            if hasattr(q, 'min_length') and hasattr(q, 'max_length'):
                if q.min_length > q.max_length:
                    errors.append(f"qc.min_length ({q.min_length}) > qc.max_length ({q.max_length})")
            if hasattr(q, 'min_length') and q.min_length < 1:
                errors.append(f"qc.min_length must be >= 1, got {q.min_length}")

        # Validate DESeq settings
        if hasattr(self, 'deseq'):
            d = self.deseq
            if hasattr(d, 'padj_threshold') and not (0 < d.padj_threshold < 1):
                errors.append(f"deseq.padj_threshold must be between 0 and 1, got {d.padj_threshold}")
            if hasattr(d, 'lfc_threshold') and d.lfc_threshold < 0:
                errors.append(f"deseq.lfc_threshold must be >= 0, got {d.lfc_threshold}")

        # Validate counting settings
        if hasattr(self, 'counting'):
            c = self.counting
            if hasattr(c, 'min_cpm') and c.min_cpm < 0:
                errors.append(f"counting.min_cpm must be >= 0, got {c.min_cpm}")
            if hasattr(c, 'min_samples_detected') and c.min_samples_detected < 1:
                errors.append(f"counting.min_samples_detected must be >= 1, got {c.min_samples_detected}")

        # Validate enrichment settings
        if hasattr(self, 'enrichment'):
            e = self.enrichment
            if hasattr(e, 'enrichment_padj') and not (0 < e.enrichment_padj < 1):
                errors.append(f"enrichment.enrichment_padj must be between 0 and 1, got {e.enrichment_padj}")

        return errors


# Create global config instance
config = Config()
