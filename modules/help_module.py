"""
Help Module for sRNAtlas
Comprehensive documentation and usage guides for the Small RNA-seq Analysis Platform
"""
import streamlit as st
from pathlib import Path


def render_help_page():
    """Render the help page"""
    st.header("📚 Help & Documentation")

    tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs([
        "🚀 Quick Start",
        "📊 Pipeline Guide",
        "🔧 Module Reference",
        "📁 File Formats",
        "❓ FAQ",
        "🔧 Troubleshooting",
        "📥 Download Docs"
    ])

    with tab1:
        render_quick_start()

    with tab2:
        render_pipeline_guide()

    with tab3:
        render_module_reference()

    with tab4:
        render_file_formats()

    with tab5:
        render_faq()

    with tab6:
        render_troubleshooting()

    with tab7:
        render_download_docs()


def render_quick_start():
    """Render quick start guide"""
    st.subheader("🚀 Quick Start Guide")

    st.markdown("""
    Get your small RNA-seq analysis running in 10 minutes!

    ---

    ### Step 1: Create a Project
    1. Click **Project** in the sidebar
    2. Go to **Project Settings** tab
    3. Enter project name and select organism

    ### Step 2: Upload Your Data
    1. Stay in **Project** module
    2. Go to **Data Upload** tab
    3. Drag & drop your `.fastq.gz` files
    4. Click **Add Files to Project**

    ### Step 3: Quality Control
    1. Click **Quality Control** in sidebar
    2. Go to **Upload Data** tab
    3. Select **Use files from Project**
    4. Click **Run QC Analysis**
    5. Review results in **Read Statistics** and **Size Distribution** tabs

    ### Step 4: Trim Adapters
    1. Click **Trimming** in sidebar
    2. In **Input Files**, select "Use files from Project"
    3. In **Settings**, select your adapter preset
    4. Click **Run Trimming**

    ### Step 5: Set Up Reference
    1. Click **Databases** in sidebar
    2. In **miRBase** tab, select your organism
    3. Click **Download from miRBase**
    4. Click **Build Bowtie Index**

    ### Step 6: Align Reads
    1. Click **Alignment** in sidebar
    2. Verify reference index is selected
    3. In **Run Alignment**, select "Use trimmed files"
    4. Click **Start Alignment**

    ### Step 7: Post-Alignment QC
    1. Click **Post-Align QC** in sidebar
    2. Click **Run Post-Alignment QC**
    3. Review mapping statistics and RNA composition

    ### Step 8: Continue Analysis
    - **Counting**: Generate count matrix
    - **DE Analysis**: Differential expression
    - **Targets**: miRNA target prediction
    - **GO/Pathway**: Functional enrichment

    ---

    ### Quick Reference: Recommended Settings

    | Parameter | miRNA Analysis | Discovery Mode |
    |-----------|---------------|----------------|
    | Min length | 18 nt | 15 nt |
    | Max length | 25 nt | 40 nt |
    | Mismatches | 0-1 | 2 |
    | Multi-map | Suppress >10 | Keep all |

    ### Expected Results

    | Metric | Good | Acceptable | Check |
    |--------|------|------------|-------|
    | Trimming pass rate | >70% | 50-70% | <50% |
    | Alignment rate | >70% | 30-70% | <30% |
    | Unique mapping | >50% | 30-50% | <30% |
    """)


def render_pipeline_guide():
    """Render comprehensive pipeline guide"""
    st.subheader("📊 Analysis Pipeline Guide")

    st.markdown("""
    ### Pipeline Overview

    ```
    ┌─────────────────┐
    │   Raw FASTQ     │
    │   Files         │
    └────────┬────────┘
             │
             ▼
    ┌─────────────────┐
    │  Quality        │ ← Pre-alignment QC
    │  Control        │   (Read quality, size distribution)
    └────────┬────────┘
             │
             ▼
    ┌─────────────────┐
    │  Adapter        │ ← Remove adapters & filter by length
    │  Trimming       │   (Cutadapt)
    └────────┬────────┘
             │
             ▼
    ┌─────────────────┐
    │  Reference      │ ← miRBase, RNAcentral, or custom
    │  Database       │   (Build Bowtie index)
    └────────┬────────┘
             │
             ▼
    ┌─────────────────┐
    │  Alignment      │ ← Map reads to reference
    │  (Bowtie)       │   (Optimized for small RNA)
    └────────┬────────┘
             │
             ▼
    ┌─────────────────┐
    │  Post-Align     │ ← Mapping QC, strand bias,
    │  QC             │   RNA type composition
    └────────┬────────┘
             │
             ▼
    ┌─────────────────┐
    │  Read           │ ← Count reads per feature
    │  Counting       │   (miRNA, gene, etc.)
    └────────┬────────┘
             │
             ▼
    ┌─────────────────┐
    │  Differential   │ ← Compare conditions
    │  Expression     │   (Statistical analysis)
    └────────┬────────┘
             │
             ▼
    ┌─────────────────┐
    │  Functional     │ ← Target prediction,
    │  Analysis       │   GO/Pathway enrichment
    └─────────────────┘
    ```

    ---

    ### Supported RNA Types

    #### Core Small RNAs
    | RNA Type | Size Range | Description |
    |----------|-----------|-------------|
    | **miRNA** | 18-25 nt | MicroRNAs - gene expression regulators, 5' U bias |
    | **siRNA** | 20-24 nt | Small interfering RNAs - RNAi pathway |
    | **piRNA** | 24-32 nt | PIWI-interacting RNAs - transposon silencing |
    | **tRF/tsRNA** | 14-40 nt | tRNA-derived fragments - stress response |
    | **rsRF** | 15-40 nt | rRNA-derived fragments |

    #### Other Small RNAs
    | RNA Type | Size Range | Description |
    |----------|-----------|-------------|
    | **snoRNA** | 60-300 nt | Small nucleolar RNA - rRNA modification |
    | **snRNA** | 100-300 nt | Small nuclear RNA - splicing (U1, U2, etc.) |
    | **Y RNA** | 80-120 nt | DNA replication, RNA quality control |
    | **vtRNA** | 80-150 nt | Vault RNA |
    | **7SL/SRP** | ~300 nt | Signal recognition particle RNA |

    #### Plant-Specific Small RNAs
    | RNA Type | Size | Description |
    |----------|------|-------------|
    | **tasiRNA** | 21 nt | Trans-acting siRNA |
    | **phasiRNA** | 21/24 nt | Phased siRNA |
    | **natsiRNA** | 21-24 nt | Natural antisense siRNA |
    | **hc-siRNA** | 24 nt | Heterochromatic siRNA - DNA methylation |

    ---

    ### Key Analysis Parameters

    #### Trimming (Cutadapt)
    | Parameter | Default | Description |
    |-----------|---------|-------------|
    | 3' Adapter | (preset) | Adapter sequence to remove |
    | Min Length | 18 | Discard reads shorter than this |
    | Max Length | 35 | Discard reads longer than this |
    | Quality Cutoff | 20 | Trim low-quality bases |

    #### Alignment (Bowtie)
    | Parameter | Default | Description |
    |-----------|---------|-------------|
    | Mismatches (-v) | 1 | Allow 0-3 mismatches |
    | Max alignments (-k) | 10 | Report up to k alignments |
    | Best mode | On | Report best alignments first |
    | Strata | On | Only report best stratum |

    #### Differential Expression
    | Parameter | Default | Description |
    |-----------|---------|-------------|
    | FDR threshold | 0.05 | Significance cutoff |
    | Log2FC threshold | 0.585 | ~1.5-fold change |
    | Min count | 10 | Minimum read count |
    """)


def render_module_reference():
    """Render detailed module reference"""
    st.subheader("🔧 Module Reference")

    with st.expander("📁 Project Module", expanded=False):
        st.markdown("""
        ### Project Management

        **Purpose:** Organize your analysis session and data files.

        **Tabs:**
        - **Project Settings**: Name, organism, metadata
        - **Data Upload**: Upload FASTQ files
        - **Save/Load**: Save and restore projects

        **Session State Variables:**
        - `project_name`: Current project name
        - `project_fastq_files`: List of uploaded FASTQ paths
        - `project_fastq_names`: Original filenames
        - `organism`: Selected organism

        **Tips:**
        - Upload all samples at once
        - Files up to 2GB supported
        - Use compressed .fastq.gz files
        """)

    with st.expander("📊 Quality Control Module", expanded=False):
        st.markdown("""
        ### Quality Control Overview

        The QC module provides comprehensive analysis at multiple stages of your pipeline.

        **Tabs:**

        | Tab | Purpose |
        |-----|---------|
        | **Upload Data** | Select FASTQ files or count matrix |
        | **Pre-Alignment Stats** | Raw read quality metrics |
        | **Size Distribution** | Length histogram with RNA type annotations |
        | **Contamination Check** | Adapter content, poly-A/T, k-mers |
        | **Post-Alignment QC** | Mapping quality, strand bias, MAPQ |
        | **Pre vs Post Comparison** | Compare metrics before/after alignment |

        ---

        ### Pre-Alignment Metrics

        | Metric | Good | Warning |
        |--------|------|---------|
        | Mean Quality | >30 | <20 |
        | GC Content | 40-60% | <30% or >70% |
        | Adapter Content | <5% | >20% |

        **Size Distribution Interpretation:**
        - Peak at 21-23 nt → miRNA enrichment
        - Peak at 26-30 nt → piRNA enrichment
        - Broad distribution → possible degradation

        ---

        ### Post-Alignment Metrics

        | Metric | Good | Warning |
        |--------|------|---------|
        | Mapping Rate | >70% | <50% |
        | Unique Rate | >60% | <40% |
        | Mean MAPQ | >30 | <20 |
        | Strand Bias | <0.1 | >0.2 |

        **Key Post-Alignment Checks:**
        - Aligned read length distribution
        - Mapping quality (MAPQ) distribution
        - Unique vs multi-mapped reads
        - Forward/reverse strand balance
        """)

    with st.expander("✂️ Trimming Module", expanded=False):
        st.markdown("""
        ### Adapter Trimming

        **Purpose:** Remove adapter sequences and filter by quality/length.

        **Supported Platforms & Adapters:**

        | Platform | Library Kits |
        |----------|-------------|
        | **Illumina** | TruSeq Small RNA, NEBNext, QIAseq, Lexogen, Clontech SMARTer, Bioo NEXTflex |
        | **Ion Torrent** | Ion S5/Proton/Genexus Small RNA, A-adapter, P1-adapter |
        | **BGI/MGI** | DNBSEQ Small RNA, MGIEasy Small RNA |
        | **Nanopore** | cDNA-PCR, Direct RNA |
        | **PacBio** | SMRTbell |
        | **Legacy** | 454/Roche, SOLiD |
        | **Special** | Poly-A tail, Poly-T head, Custom |

        **Parameters:**

        | Parameter | Recommended | Notes |
        |-----------|-------------|-------|
        | Min length | 18 nt | Captures miRNAs |
        | Max length | 35 nt | For miRNA/siRNA (increase to 200 for snoRNA) |
        | Quality cutoff | 20 | Phred score |
        | Error rate | 0.1 (10%) | Adapter matching tolerance |

        **Troubleshooting 0% pass rate:**
        - Wrong adapter/platform selected
        - Adapters already trimmed in raw data
        - File format issues
        - Try "Custom" and enter adapter manually
        """)

    with st.expander("🗄️ Databases Module", expanded=False):
        st.markdown("""
        ### Reference Databases

        **Purpose:** Set up reference sequences for alignment.

        **Sources:**
        - **miRBase**: Mature miRNA sequences
        - **RNAcentral**: Comprehensive ncRNA database
        - **Custom**: Your own FASTA files

        **Index Building:**
        1. Download/upload reference FASTA
        2. Click "Build Bowtie Index"
        3. Wait for completion (1-5 minutes)

        **Organism Codes (miRBase):**
        - hsa = Human
        - mmu = Mouse
        - rno = Rat
        - dre = Zebrafish
        - dme = Drosophila
        - cel = C. elegans
        - ath = Arabidopsis
        """)

    with st.expander("🔗 Alignment Module", expanded=False):
        st.markdown("""
        ### Read Alignment

        **Purpose:** Map trimmed reads to reference using Bowtie.

        **Bowtie Parameters:**
        | Flag | Default | Description |
        |------|---------|-------------|
        | -v | 1 | Mismatches allowed (0-3) |
        | -k | 10 | Report up to k alignments |
        | --best | On | Report best alignments first |
        | --strata | On | Only best stratum |
        | -m | 0 | Suppress if >m hits (0=off) |

        **Output Files:**
        - `sample_sorted.bam`: Aligned reads
        - `sample_sorted.bam.bai`: BAM index
        - `sample_unaligned.fq`: Unaligned reads
        - `sample_alignment.log`: Statistics

        **Interpreting Alignment Rates:**
        - >70%: Excellent
        - 30-70%: Acceptable
        - <30%: Check reference/trimming
        """)

    with st.expander("🔬 Post-Alignment QC & Comparison", expanded=False):
        st.markdown("""
        ### Post-Alignment Quality Control

        **Purpose:** Assess alignment quality and compare to pre-alignment metrics.

        **Post-Alignment Metrics:**

        | Metric | Description | Threshold |
        |--------|-------------|-----------|
        | Mapping Rate | % reads aligned to reference | >70% good |
        | Unique Rate | % uniquely mapped (NH:i:1) | >60% good |
        | Multi-Mapped | Reads with multiple alignments | Monitor |
        | Mean MAPQ | Average mapping quality | >30 good |
        | Strand Bias | Balance of +/- strand | <0.1 good |

        ---

        ### Pre vs Post Comparison

        The comparison view shows how your data changes through the pipeline:

        **Read Count Flow:**
        - Input Reads → Aligned → Not Aligned → Suppressed

        **Length Changes:**
        - Compare raw read lengths to aligned read lengths
        - Large differences may indicate trimming or soft-clipping

        **Quality Indicators:**
        - Low alignment + good pre-QC → Check reference database
        - High suppression → Reduce -m parameter or use broader reference
        - Strand bias → May indicate library prep issues

        ---

        ### RNA Composition by Length

        | Read Length | RNA Type |
        |-------------|----------|
        | 18-25 nt | miRNA-like |
        | 20-24 nt | siRNA-like |
        | 24-32 nt | piRNA-like |
        | 28-40 nt | tRNA fragments |

        **5' Nucleotide Bias:**
        miRNAs typically show enrichment for Uracil (T in DNA) at the 5' end.
        """)

    with st.expander("📈 Counting Module", expanded=False):
        st.markdown("""
        ### Read Counting

        **Purpose:** Generate count matrices from BAM files.

        **Methods:**
        - Feature counting per reference sequence
        - Multi-mapper handling options:
          - Unique only
          - Fractional (divide by NH tag)
          - All alignments

        **Output:**
        - Count matrix (features × samples)
        - Raw counts for DE analysis
        """)

    with st.expander("🔬 DE Analysis Module", expanded=False):
        st.markdown("""
        ### Differential Expression Analysis

        **Purpose:** Identify differentially expressed small RNAs.

        **Methods:**
        - DESeq2 (recommended)
        - edgeR
        - Basic (fold-change + t-test)

        **Required Input:**
        - Count matrix
        - Sample metadata with condition column

        **Output:**
        - DE results table (log2FC, p-value, FDR)
        - MA plot
        - Volcano plot
        - Heatmap
        - PCA plot
        """)

    with st.expander("🎯 Target Prediction Module", expanded=False):
        st.markdown("""
        ### miRNA Target Prediction

        **Purpose:** Identify potential target genes for miRNAs.

        **Methods:**
        - **psRNATarget**: Plant-specific target prediction
        - **Local seed matching**: Sequence complementarity analysis
        - **TargetScan**: Animal miRNA targets (external)

        **Parameters:**
        | Parameter | Default | Description |
        |-----------|---------|-------------|
        | Seed region | 2-8 nt | Core binding region |
        | Max mismatches | 3 | Allowed in seed region |
        | Min score | 3.0 | Prediction confidence |

        **Output:**
        - Target gene list per miRNA
        - Binding site positions
        - Prediction scores
        - Target gene annotations
        """)

    with st.expander("🧬 GO/Pathway Module", expanded=False):
        st.markdown("""
        ### Functional Enrichment Analysis

        **Purpose:** Identify enriched biological functions and pathways.

        **Gene Ontology (GO) Analysis:**
        - **Biological Process (BP)**: Cellular functions
        - **Molecular Function (MF)**: Biochemical activities
        - **Cellular Component (CC)**: Subcellular locations

        **Pathway Analysis:**
        - **KEGG**: Metabolic and signaling pathways
        - **Reactome**: Curated pathway database

        **Input:**
        - DE miRNA list (or target genes)
        - Background gene set

        **Parameters:**
        | Parameter | Default | Description |
        |-----------|---------|-------------|
        | FDR threshold | 0.05 | Significance cutoff |
        | Min gene set | 10 | Minimum genes per term |
        | Max gene set | 500 | Maximum genes per term |

        **Output:**
        - Enriched GO terms table
        - KEGG pathway list
        - Dot plots and bar charts
        - Network visualization
        """)

    with st.expander("🔍 Novel miRNA Module", expanded=False):
        st.markdown("""
        ### Novel miRNA Discovery

        **Purpose:** Identify unannotated small RNAs with miRNA-like characteristics from unaligned reads.

        **Input:**
        - Unaligned reads from alignment step (FASTQ)
        - Reference genome (FASTA) for hairpin prediction

        **Discovery Criteria:**
        | Criterion | Default | Description |
        |-----------|---------|-------------|
        | Length | 20-24 nt | Typical miRNA size range |
        | Min reads | 100 | Minimum abundance threshold |
        | 5' U bias | Yes | Prefer sequences starting with U/T |
        | Hairpin MFE | -20 kcal/mol | Minimum folding energy |

        **Scoring Factors:**
        - Read abundance (higher = better)
        - 5' nucleotide preference (U/T = miRNA-like)
        - GC content (typical range: 40-60%)
        - Hairpin structure potential

        **Output:**
        - Candidate novel miRNA sequences
        - Expression counts
        - Characteristic scores
        - Downloadable CSV table
        """)

    with st.expander("🧫 isomiR Module", expanded=False):
        st.markdown("""
        ### isomiR Analysis

        **Purpose:** Detect and analyze miRNA isoforms (sequence variants).

        **What are isomiRs?**
        isomiRs are sequence variants of canonical miRNAs that differ by:
        - **5' variants**: Different start position (affects target specificity!)
        - **3' variants**: Different end position (most common)
        - **SNPs**: Internal nucleotide changes
        - **Non-templated additions (NTA)**: Added nucleotides (often A or U)

        **Input:**
        - Aligned BAM files from alignment step
        - Reference miRNA sequences

        **Detection Settings:**
        | Parameter | Default | Description |
        |-----------|---------|-------------|
        | Min reads | 10 | Minimum reads per isomiR |
        | Max 5' diff | 2 nt | Maximum 5' position change |
        | Max 3' diff | 4 nt | Maximum 3' position change |

        **Variant Types Detected:**
        - `canonical`: Matches reference exactly
        - `3p_addition`: Extended at 3' end
        - `3p_trimming`: Shortened at 3' end
        - `snp`: Single nucleotide polymorphism
        - `nta`: Non-templated addition (A/U)

        **Output:**
        - isomiR table with variant classifications
        - Variant type distribution (pie chart)
        - Per-miRNA isomiR profiles
        - Sequence alignment visualization
        """)

    with st.expander("⚡ Batch Module", expanded=False):
        st.markdown("""
        ### Batch Processing

        **Purpose:** Automate the full analysis pipeline for multiple samples.

        **Pipeline Steps:**
        1. Quality Control
        2. Adapter Trimming
        3. Alignment
        4. Post-Alignment QC
        5. Read Counting
        6. (Optional) DE Analysis

        **Configuration:**
        - Select which steps to run
        - Configure parameters for each step
        - Set up sample groups for DE analysis

        **Monitoring:**
        - Real-time progress tracking
        - Per-sample status updates
        - Error logging and reporting

        **Output:**
        - All intermediate files
        - Combined count matrix
        - Summary statistics
        - HTML batch report
        """)

    with st.expander("📋 Reports Module", expanded=False):
        st.markdown("""
        ### Report Generation

        **Purpose:** Generate comprehensive analysis reports and export results.

        **Report Types:**
        - **HTML Report**: Interactive, self-contained report
        - **PDF Report**: Print-ready summary
        - **Excel Export**: All results in spreadsheet format

        **Report Contents:**
        - Analysis parameters and settings
        - QC metrics and plots
        - Alignment statistics
        - DE analysis results
        - Enrichment findings
        - Methods section (for publications)

        **Figure Export:**
        - PNG: High-resolution images
        - SVG: Vector graphics (scalable)
        - HTML: Interactive Plotly figures

        **Project Export:**
        - Save complete project state
        - Export as ZIP archive
        - Share with collaborators
        """)

    with st.expander("⚙️ Settings Module", expanded=False):
        st.markdown("""
        ### Configuration Settings

        **Purpose:** Centralized configuration management with presets.

        **Preset Profiles:**
        | Preset | Use Case |
        |--------|----------|
        | miRNA Strict | Mature miRNA analysis (0 mismatches) |
        | miRNA Standard | Typical miRNA profiling |
        | Discovery | Novel small RNA discovery |
        | Plant miRNA | Plant-specific settings |

        **Configurable Parameters:**
        - Alignment settings (mismatches, multi-mapping)
        - Trimming parameters (adapters, length filters)
        - QC thresholds
        - DE analysis cutoffs
        - Enrichment settings

        **Features:**
        - Apply presets with one click
        - Export settings to JSON/YAML
        - Import settings from file
        - Generate methods text for publications
        - Reset to defaults

        **Methods Text Generator:**
        Automatically generates a methods paragraph describing your analysis parameters for inclusion in publications.
        """)


def render_file_formats():
    """Render file format documentation"""
    st.subheader("📁 File Formats")

    st.markdown("""
    ### Input Formats

    #### FASTQ Files
    Standard sequencing format with quality scores.
    ```
    @SEQ_ID
    GATTTGGGGTTCAAAGCAGTATCGATCAAA
    +
    !''*((((***+))%%%++)(%%%%).1***
    ```
    - Extensions: `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`
    - Gzip compression recommended

    ---

    #### FASTA Files (Reference)
    ```
    >hsa-miR-21-5p MIMAT0000076
    UAGCUUAUCAGACUGAUGUUGA
    >hsa-miR-155-5p MIMAT0000646
    UUAAUGCUAAUCGUGAUAGGGGU
    ```
    - Extensions: `.fasta`, `.fa`, `.fna`

    ---

    #### Count Matrix (CSV)
    ```csv
    ,Sample1,Sample2,Sample3,Sample4
    hsa-miR-21-5p,1500,1450,2800,2750
    hsa-miR-155-5p,200,180,650,620
    hsa-let-7a-5p,3500,3200,3100,2900
    ```
    - First column: Feature IDs
    - Other columns: Sample counts
    - Tab or comma separated

    ---

    #### Sample Metadata (CSV)
    ```csv
    sample,condition,batch,replicate
    Sample1,control,1,1
    Sample2,control,1,2
    Sample3,treatment,2,1
    Sample4,treatment,2,2
    ```
    - `sample`: Must match count matrix columns
    - `condition`: Required for DE analysis
    - Additional columns: Optional covariates

    ---

    ### Output Formats

    | Output | Format | Description |
    |--------|--------|-------------|
    | Trimmed reads | FASTQ.gz | Adapter-free, filtered reads |
    | Alignments | BAM | Binary alignment map |
    | Count matrix | CSV | Read counts per feature |
    | DE results | CSV/Excel | Fold changes, p-values |
    | Plots | PNG/HTML | Publication-ready figures |
    | Reports | HTML | Comprehensive analysis report |
    """)


def render_faq():
    """Render FAQ section"""
    st.subheader("❓ Frequently Asked Questions")

    with st.expander("What sequencing depth do I need?"):
        st.markdown("""
        **Recommendations:**
        - **miRNA profiling**: 5-10 million reads/sample
        - **Discovery**: 20+ million reads/sample
        - **Differential expression**: 10+ million reads/sample

        More depth improves detection of low-abundance RNAs.
        """)

    with st.expander("Can I use paired-end data?"):
        st.markdown("""
        Small RNA-seq is typically **single-end** because:
        - Small RNAs are shorter than read length
        - Paired-end provides no additional information

        If you have paired-end data, use only **Read 1**.
        """)

    with st.expander("Why use Bowtie instead of Bowtie2?"):
        st.markdown("""
        **Bowtie (v1)** is preferred for small RNA because:
        - Optimized for short reads (<50 bp)
        - Exact mismatch control with `-v` parameter
        - Better handling of multi-mappers

        **Bowtie2** is designed for longer reads and uses
        different alignment algorithms.
        """)

    with st.expander("How do I handle multi-mapping reads?"):
        st.markdown("""
        Multi-mapping is common in small RNA-seq because:
        - miRNA families have similar sequences
        - Paralogs and gene duplications

        **Strategies:**
        1. **Keep all** (-a): Report all alignments
        2. **Report best** (--best): Only best alignments
        3. **Limit reporting** (-k 10): Up to k alignments
        4. **Suppress** (-m 10): Discard if >m alignments

        Default: `-k 10 --best --strata`
        """)

    with st.expander("What does 5' U bias mean?"):
        st.markdown("""
        **5' Uracil (U) bias** is characteristic of miRNAs:
        - Argonaute proteins prefer U at the 5' end
        - ~70% of miRNAs start with U

        In sequencing data (DNA), this appears as **T enrichment**.

        **Interpretation:**
        - Strong T bias at 5' → Good miRNA library
        - No bias → May contain other RNA types
        """)

    with st.expander("How do I interpret alignment rates?"):
        st.markdown("""
        | Alignment Rate | Interpretation |
        |----------------|----------------|
        | >70% | Excellent - well-matched reference |
        | 50-70% | Good - typical for miRNA libraries |
        | 30-50% | Acceptable - may have other RNAs |
        | <30% | Check reference and trimming |

        **Low alignment causes:**
        - Wrong reference organism
        - Adapters not trimmed
        - Contamination (rRNA, degradation)
        - Novel/unannotated small RNAs
        """)

    with st.expander("What is the difference between miRBase and RNAcentral?"):
        st.markdown("""
        **miRBase:**
        - Specialized miRNA database
        - Curated mature and precursor miRNAs
        - ~2,600 human miRNAs
        - Best for miRNA-specific studies

        **RNAcentral:**
        - Comprehensive ncRNA database
        - Includes miRNA, tRNA, rRNA, snoRNA, etc.
        - Aggregates from multiple databases
        - Best for broad small RNA studies
        """)

    with st.expander("What are isomiRs and why do they matter?"):
        st.markdown("""
        **isomiRs** are sequence variants of canonical miRNAs:
        - **5' isomiRs**: Different start position - affects seed region and target specificity
        - **3' isomiRs**: Different end position - most common, less impact on targeting
        - **SNP variants**: Internal nucleotide changes
        - **NTA (non-templated additions)**: Added nucleotides, usually A or U

        **Why they matter:**
        - Can have different target genes than canonical miRNA
        - May have distinct biological functions
        - Important for accurate quantification
        - Can indicate RNA editing or processing changes
        """)

    with st.expander("How many replicates do I need?"):
        st.markdown("""
        **Minimum recommendations:**
        - **Exploratory studies**: 3 biological replicates per condition
        - **Publication-quality DE**: 4-6 replicates per condition
        - **Clinical studies**: 10+ samples per group

        **Why replicates matter:**
        - Capture biological variability
        - Improve statistical power
        - Reduce false positives
        - Required for meaningful p-values

        **Technical vs Biological replicates:**
        - Technical: Same sample, sequenced multiple times (less valuable)
        - Biological: Different individuals/samples (essential for DE)
        """)

    with st.expander("What normalization method should I use?"):
        st.markdown("""
        **DESeq2 Median Ratio (Default):**
        - Accounts for library size and composition
        - Robust to outliers
        - Best for most small RNA-seq experiments

        **CPM (Counts Per Million):**
        - Simple library size normalization
        - Good for visualization
        - Less robust for DE analysis

        **TPM/RPKM:**
        - Not recommended for small RNA-seq
        - Designed for mRNA with varying lengths
        - Small RNAs are similar length

        **When to use each:**
        | Method | Use Case |
        |--------|----------|
        | DESeq2 | Differential expression |
        | CPM | Visualization, filtering |
        | Raw counts | Input to DE tools |
        """)

    with st.expander("How do I choose significance thresholds?"):
        st.markdown("""
        **Common thresholds:**
        | Parameter | Strict | Standard | Relaxed |
        |-----------|--------|----------|---------|
        | FDR (padj) | 0.01 | 0.05 | 0.1 |
        | log2FC | 1.0 | 0.585 | 0.0 |
        | Fold Change | 2x | 1.5x | Any |

        **Considerations:**
        - More stringent = fewer false positives, may miss real changes
        - More relaxed = more candidates, higher false positive rate
        - Adjust based on validation capacity
        - Consider effect size (log2FC) not just p-value
        """)

    with st.expander("Can I analyze degraded RNA samples?"):
        st.markdown("""
        **Small RNAs are relatively stable**, but degradation can affect results:

        **Signs of degradation:**
        - Broad size distribution (no clear peaks)
        - High proportion of very short reads
        - Low alignment rates
        - Inconsistent results across replicates

        **Recommendations:**
        - Check RIN/RQN values before sequencing
        - Use appropriate size selection during library prep
        - Consider excluding severely degraded samples
        - Note: tRNA fragments increase with degradation
        """)

    with st.expander("What is the seed region and why is it important?"):
        st.markdown("""
        **Seed region**: Nucleotides 2-8 from the 5' end of the miRNA

        **Why it's critical:**
        - Primary determinant of target recognition
        - Perfect complementarity usually required
        - ~16,000 possible seed sequences (4^7)
        - miRNAs with same seed = same seed family

        **Target prediction:**
        - Seed match is the minimum requirement
        - Additional 3' pairing improves binding
        - 5' isomiRs change the seed = different targets!

        ```
        Position:  1 2 3 4 5 6 7 8 9 ...
        miRNA:     U G A G G U A G U ...
                     ─────────────
                        SEED REGION
        ```
        """)

    with st.expander("How do I handle batch effects?"):
        st.markdown("""
        **Batch effects** are technical variations between experimental batches:

        **Prevention (best approach):**
        - Randomize samples across batches
        - Process all samples together when possible
        - Use same reagent lots

        **Detection:**
        - PCA plot colored by batch
        - Hierarchical clustering
        - Check if samples cluster by batch vs condition

        **Correction methods:**
        - Include batch as covariate in DE model
        - ComBat (for severe batch effects)
        - SVA (Surrogate Variable Analysis)

        **In sRNAtlas:**
        - Add 'batch' column to sample metadata
        - Include in DE analysis design
        """)

    with st.expander("What file formats does sRNAtlas support?"):
        st.markdown("""
        **Input formats:**
        | Format | Extension | Description |
        |--------|-----------|-------------|
        | FASTQ | .fastq, .fq, .fastq.gz | Raw sequencing reads |
        | FASTA | .fasta, .fa | Reference sequences |
        | BAM | .bam | Aligned reads |
        | CSV | .csv | Count matrices, metadata |

        **Output formats:**
        | Format | Description |
        |--------|-------------|
        | FASTQ.gz | Trimmed reads |
        | BAM/BAI | Alignments + index |
        | CSV/Excel | Count matrices, DE results |
        | PNG/SVG | Publication figures |
        | HTML | Interactive reports |
        """)

    with st.expander("Can I use sRNAtlas for non-model organisms?"):
        st.markdown("""
        **Yes!** sRNAtlas supports any organism with:

        **Option 1: Use related species**
        - Align to closest available reference
        - Expect lower alignment rates
        - Good for conserved miRNAs

        **Option 2: Custom reference**
        - Upload your own FASTA file
        - Use predicted miRNAs from your genome
        - Combine with miRBase for related species

        **Option 3: Discovery mode**
        - Use broader reference (RNAcentral)
        - Enable novel miRNA discovery
        - Identify unannotated small RNAs

        **Tips:**
        - Start with relaxed alignment parameters
        - Check size distribution for expected peaks
        - Validate candidates experimentally
        """)

    with st.expander("How long does analysis take?"):
        st.markdown("""
        **Typical processing times** (8-core machine, 10M reads/sample):

        | Step | Time per Sample |
        |------|-----------------|
        | QC | 30 seconds |
        | Trimming | 1-2 minutes |
        | Index building | 1-5 minutes (once) |
        | Alignment | 2-5 minutes |
        | Counting | 1 minute |
        | DE Analysis | 1-2 minutes |

        **Full pipeline**: ~10-15 minutes per sample

        **Factors affecting speed:**
        - Read count (more reads = longer)
        - Reference size
        - CPU cores available
        - Storage speed (SSD vs HDD)
        """)

    with st.expander("How do I cite sRNAtlas?"):
        st.markdown("""
        If you use sRNAtlas in your research, please cite:

        ```
        sRNAtlas: A Comprehensive Platform for Small RNA-seq Analysis
        [Version BETA]
        https://github.com/YOUR_USERNAME/sRNAtlas
        ```

        Also cite the underlying tools:
        - **Bowtie**: Langmead et al., Genome Biology 2009
        - **Cutadapt**: Martin, EMBnet.journal 2011
        - **DESeq2/pyDESeq2**: Love et al., Genome Biology 2014
        - **miRBase**: Kozomara et al., NAR 2019
        """)


def render_troubleshooting():
    """Render troubleshooting section"""
    st.subheader("🔧 Troubleshooting Guide")

    st.markdown("### Common Issues and Solutions")

    with st.expander("🔴 0 reads after trimming", expanded=True):
        st.markdown("""
        **Symptoms:** Trimming shows 0 reads processed, 0% pass rate

        **Causes & Solutions:**

        1. **Wrong adapter sequence**
           - Try different adapter presets
           - Check your library prep kit documentation
           - Use "Custom" and enter the correct sequence

        2. **Files already trimmed**
           - Skip trimming step
           - Go directly to alignment

        3. **File format issues**
           - Verify files are valid FASTQ
           - Check for corruption during upload
           - Try re-uploading files

        4. **Path/permission issues**
           - Check debug information in Results tab
           - Verify temp directory is writable
        """)

    with st.expander("🔴 0% alignment rate"):
        st.markdown("""
        **Symptoms:** Alignment completes but no reads align

        **Causes & Solutions:**

        1. **Wrong reference database**
           - Verify organism matches your samples
           - Check miRBase organism code

        2. **Reference index not built**
           - Go to Databases module
           - Click "Build Bowtie Index"
           - Wait for completion

        3. **Reads still have adapters**
           - Check trimming results
           - Verify adapter removal was successful

        4. **Mismatches too strict**
           - Increase allowed mismatches (-v 2)
           - Use --best mode

        5. **Bowtie command syntax issue**
           - Check alignment logs for errors
        """)

    with st.expander("🟡 Low alignment rate (<30%)"):
        st.markdown("""
        **Possible Causes:**

        1. **Contamination**
           - High rRNA content
           - Degraded samples
           - Run Contamination Check in QC module

        2. **Novel small RNAs**
           - Not in reference database
           - Consider custom database

        3. **Organism mismatch**
           - Partial alignment to related species
           - Verify correct organism selected

        4. **Technical issues**
           - Adapter dimers
           - Poor library quality
        """)

    with st.expander("🟡 Memory errors"):
        st.markdown("""
        **Solutions:**

        1. **Process fewer samples at once**
           - Upload in batches of 5-10

        2. **Use compressed files**
           - .fastq.gz uses less memory than .fastq

        3. **Increase system resources**
           - Close other applications
           - Add more RAM

        4. **Reduce analysis scope**
           - Subsample reads for initial QC
           - Use smaller reference database
        """)

    with st.expander("🟡 Tool not found errors"):
        st.markdown("""
        **"bowtie: command not found"**
        ```bash
        # Check if installed
        which bowtie

        # Install via conda
        conda install -c bioconda bowtie

        # Or via homebrew (macOS)
        brew install bowtie
        ```

        **"samtools: command not found"**
        ```bash
        conda install -c bioconda samtools
        ```

        **"cutadapt: command not found"**
        ```bash
        pip install cutadapt
        ```

        **Verify PATH:**
        ```bash
        echo $PATH
        # Add tool directory to PATH if needed
        export PATH=$PATH:/path/to/tools
        ```
        """)

    with st.expander("🟡 pysam installation issues"):
        st.markdown("""
        **Error:** "No module named 'pysam'"

        **Solutions:**

        ```bash
        # Standard install
        pip install pysam

        # If that fails (macOS)
        brew install htslib
        pip install pysam

        # Using conda (most reliable)
        conda install -c bioconda pysam
        ```
        """)


def render_download_docs():
    """Render documentation download section"""
    st.subheader("📥 Download Documentation")

    st.markdown("""
    Download the full documentation for offline reference.
    """)

    # Check if docs exist
    docs_dir = Path(__file__).parent.parent / "docs"

    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown("### 📘 User Guide")
        st.markdown("Complete documentation covering all features")
        user_guide = docs_dir / "USER_GUIDE.md"
        if user_guide.exists():
            content = user_guide.read_text()
            st.download_button(
                "📥 Download User Guide",
                content,
                file_name="USER_GUIDE.md",
                mime="text/markdown"
            )
        else:
            st.info("User guide not found")

    with col2:
        st.markdown("### 🚀 Quick Start")
        st.markdown("Get running in 10 minutes")
        quick_start = docs_dir / "QUICK_START.md"
        if quick_start.exists():
            content = quick_start.read_text()
            st.download_button(
                "📥 Download Quick Start",
                content,
                file_name="QUICK_START.md",
                mime="text/markdown"
            )
        else:
            st.info("Quick start not found")

    with col3:
        st.markdown("### 💻 Installation")
        st.markdown("Detailed setup instructions")
        install_guide = docs_dir / "INSTALLATION.md"
        if install_guide.exists():
            content = install_guide.read_text()
            st.download_button(
                "📥 Download Installation Guide",
                content,
                file_name="INSTALLATION.md",
                mime="text/markdown"
            )
        else:
            st.info("Installation guide not found")

    st.divider()

    st.markdown("### 📚 All Documentation")

    # Create combined documentation
    all_docs = """# sRNAtlas - Complete Documentation

Comprehensive Small RNA-seq Analysis Platform

"""
    for doc_name in ["QUICK_START.md", "USER_GUIDE.md", "INSTALLATION.md"]:
        doc_path = docs_dir / doc_name
        if doc_path.exists():
            all_docs += f"\n\n---\n\n# {doc_name.replace('.md', '').replace('_', ' ')}\n\n"
            all_docs += doc_path.read_text()

    st.download_button(
        "📥 Download All Documentation (Combined)",
        all_docs,
        file_name="sRNAtlas_Documentation.md",
        mime="text/markdown",
        type="primary"
    )

    st.divider()

    st.markdown("""
    ### 🔗 Additional Resources

    - **GitHub Repository**: Source code and issue tracking
    - **Video Tutorials**: Coming soon
    - **Example Datasets**: Available in the repository

    ### 📧 Contact

    For questions, bug reports, or feature requests:
    - Open an issue on GitHub
    - Email: support@example.com
    """)
