# sRNAtlas User Guide

Comprehensive Small RNA-seq Analysis Platform

## Table of Contents

1. [Introduction](#introduction)
2. [Getting Started](#getting-started)
3. [Analysis Pipeline](#analysis-pipeline)
4. [Module Reference](#module-reference)
5. [Best Practices](#best-practices)
6. [Troubleshooting](#troubleshooting)

---

## Introduction

sRNAtlas is a web-based platform for comprehensive small RNA sequencing (sRNA-seq) analysis. It provides an intuitive interface for processing raw sequencing data through to differential expression and functional enrichment analysis.

### Key Features

- **Multi-organism support**: Plants and animals via miRBase and RNAcentral
- **Complete pipeline**: From raw FASTQ to functional enrichment
- **Adapter trimming**: Cutadapt integration with common presets
- **Optimized alignment**: Bowtie v1 for short read mapping
- **Robust statistics**: pyDESeq2 for differential expression
- **Target prediction**: psRNATarget (plants), miRanda (animals)
- **Batch processing**: Automate full pipeline runs
- **Project management**: Save and restore analysis sessions
- **QC Scorecard**: Traffic-light quality assessment with outlier detection
- **Provenance tracking**: Full reproducibility with YAML/JSON export
- **isomiR analysis**: Differential usage and arm switching detection
- **Performance caching**: Fast repeat analyses with Streamlit caching

### Supported RNA Types

| RNA Type | Size Range | Description |
|----------|-----------|-------------|
| miRNA | 18-25 nt | MicroRNAs - gene expression regulators |
| siRNA | 20-24 nt | Small interfering RNAs |
| piRNA | 24-32 nt | PIWI-interacting RNAs |
| tRNA fragments | 28-36 nt | Transfer RNA-derived fragments |

---

## Getting Started

### System Requirements

- Python 3.9+
- 8+ GB RAM recommended
- External tools: Bowtie, Samtools, Cutadapt

### Quick Start

1. **Create Project**: Set up project name and organism
2. **Upload Data**: Add FASTQ files
3. **Quality Control**: Assess read quality
4. **Trim Adapters**: Remove adapter sequences
5. **Build Reference**: Download and index reference database
6. **Align Reads**: Map to reference with Bowtie
7. **Post-Alignment QC**: Verify alignment quality
8. **Count Reads**: Generate count matrix
9. **DE Analysis**: Identify differential expression
10. **Functional Analysis**: GO/pathway enrichment

---

## Analysis Pipeline

### Pipeline Overview

```
Raw FASTQ → QC → Trimming → Alignment → Post-Align QC → Counting → DE → Enrichment
```

### Step 1: Quality Control

Pre-alignment quality assessment:
- Total read counts
- Read length distribution
- Quality score profiles
- GC content analysis
- Adapter content detection
- Contamination screening

**Key metrics to check:**
- Mean quality score > 30
- GC content 40-60%
- Clear size peak (21-23 nt for miRNA)

### Step 2: Adapter Trimming

Remove sequencing adapters and filter reads:
- Supported adapters: Illumina TruSeq, NEBNext, QIAseq
- Custom adapter sequences supported
- Length filtering (default: 18-35 nt)
- Quality trimming

**Recommended settings:**
- Min length: 18 nt (captures mature miRNAs)
- Max length: 35 nt (includes piRNAs)
- Quality cutoff: 20

### Step 3: Reference Database

Set up alignment reference:
- **miRBase**: Curated miRNA sequences
- **RNAcentral**: Comprehensive ncRNA database
- **Custom FASTA**: Your own sequences

After downloading, build Bowtie index before alignment.

### Step 4: Alignment

Map reads using Bowtie (optimized for small RNA):
- Uses `-v` mode for exact mismatch control
- Best-stratum alignment for multi-mappers
- Configurable multi-mapping handling

**Default parameters:**
- Mismatches: 1
- Report up to 10 alignments
- Best mode + Strata enabled

### Step 5: Post-Alignment QC

Assess alignment quality:
- Mapping rate (target: >50%)
- Unique vs multi-mapped reads
- Strand distribution
- 5' nucleotide bias (U enrichment indicates miRNA)
- RNA type composition by length

### Step 6: Read Counting

Generate count matrices:
- Counts per reference sequence
- Multi-mapper handling options
- Export as CSV for downstream analysis

### Step 7: Differential Expression

Compare conditions with pyDESeq2:
- Normalized counts
- Log2 fold changes
- Statistical testing with FDR correction
- Visualization: MA plot, volcano plot, heatmap

### Step 8: Functional Enrichment

Identify enriched pathways:
- Gene Ontology (BP, MF, CC)
- KEGG pathways
- miRNA target prediction

---

## Module Reference

### Project Module
Organize your analysis: Create/name projects, upload FASTQ files, set organism and metadata, save/load project sessions.

### Quality Control Module
Tabs: Upload Data, Read Statistics, Size Distribution, Contamination Check

### Trimming Module
Adapter removal with Cutadapt: Preset adapter sequences, custom adapter input, length and quality filtering.

### Databases Module
Reference management: Download from miRBase/RNAcentral, upload custom FASTA, build Bowtie indexes.

### Alignment Module
Bowtie alignment: Configure mismatch parameters, multi-mapping options, output BAM files.

### Post-Alignment QC Module
Alignment quality: Mapping statistics, length distributions, 5' nucleotide bias, RNA composition.

### Counting Module
Read quantification: Per-feature counting, count matrix export.

**Multi-mapper Counting Modes:**
- **All alignments**: Count every alignment (default)
- **Unique only**: Count only reads with NH=1
- **Fractional**: Weight each alignment by 1/n (where n = number of alignments)
- **Primary only**: Count only primary alignments (SAM flag)

### DE Analysis Module
Differential expression: Sample metadata input, pyDESeq2 analysis, visualization.

### Target Prediction Module
miRNA targets with multiple algorithms:

**For Plants:**
- **psRNATarget** (API): Well-established plant target prediction
- **Local Seed Matching**: Fast algorithm based on seed complementarity

**For Animals:**
- **miRanda**: Thermodynamics-based algorithm with score/energy thresholds
- **Local Seed Matching**: Fallback when miRanda not installed

**Parameters:**
- Score threshold (miRanda): Minimum alignment score (default: 140)
- Energy threshold (miRanda): Maximum free energy (default: -20 kcal/mol)
- Strict seed pairing: Require perfect seed complementarity

### GO/Pathway Module
Functional enrichment: GO term analysis, KEGG pathway enrichment.

### Batch Module
Pipeline automation: Configure full pipeline, process multiple samples.

### Reports Module
Generate reports: HTML summaries, figure export, ZIP archives.

**Provenance Tracking Tab:**
- Records all pipeline parameters and tool versions
- Tracks file checksums (MD5) for reproducibility
- Exports as YAML or JSON format
- Shows Python package versions
- Automatically tracks each analysis step

### isomiR Module (Advanced)

**Differential Usage Tab:**
Compare isomiR ratios between conditions to identify processing changes.
- Define sample groups (from metadata or manually)
- Statistical testing with FDR correction
- Visualization of top differential isomiRs

**Arm Switching Tab:**
Detect changes in 5p/3p arm dominance between conditions.
- Automatic detection of 5p/3p miRNA pairs
- T-test on log-transformed ratios
- Identifies statistically significant arm switches

### QC Scorecard

Traffic-light quality assessment for each sample:

| Status | Meaning |
|--------|---------|
| ✅ OK | Metric within acceptable range |
| ⚠️ WARNING | Metric approaching threshold |
| ❌ CRITICAL | Metric outside acceptable range |

**Default Thresholds:**
| Metric | Warning | Critical |
|--------|---------|----------|
| Total Reads | < 1M | < 100K |
| Mean Quality | < 28 | < 20 |
| Alignment Rate | < 50% | < 20% |
| rRNA % | > 10% | > 30% |

**Multi-sample Outlier Detection:**
- MAD (Median Absolute Deviation) based detection
- PCA clustering for batch effect visualization
- Overlay distribution plots

---

## Best Practices

### Quality Thresholds

| Metric | Good | Acceptable | Review |
|--------|------|------------|--------|
| Trim pass rate | >70% | 50-70% | <50% |
| Alignment rate | >70% | 30-70% | <30% |
| Unique mapping | >50% | 30-50% | <30% |

---

## Troubleshooting

### 0 reads after trimming
- Wrong adapter sequence selected
- Files already trimmed
- Try different adapter presets

### 0% alignment rate
- Wrong reference organism
- Index not built
- Verify organism selection and build Bowtie index

### Low alignment rate (<30%)
- rRNA contamination
- Novel small RNAs not in reference
- Run contamination check

---

## Citation

If you use sRNAtlas in your research, please cite:

```
sRNAtlas: A Comprehensive Platform for Small RNA-seq Analysis
```

---

## FAQ

### Q: Why is my analysis slow the first time but fast on repeat?
**A:** sRNAtlas uses Streamlit caching (`@st.cache_data`) to store results of heavy computations. The first run calculates everything; subsequent runs with the same inputs use cached results.

### Q: How do I interpret the QC Scorecard colors?
**A:**
- ✅ Green (OK): Your sample passes quality thresholds
- ⚠️ Yellow (WARNING): Approaching thresholds, review recommended
- ❌ Red (CRITICAL): Below acceptable thresholds, investigate before proceeding

### Q: What's the difference between multi-mapper counting modes?
**A:**
- **All**: Every alignment counted (may inflate counts for multi-mappers)
- **Unique**: Only reads mapping to one location (most conservative)
- **Fractional**: Each alignment gets 1/n weight (balanced approach)
- **Primary**: Only the best alignment per read

### Q: How do I export my analysis for reproducibility?
**A:** Go to Reports → Provenance tab → Download as YAML or JSON. This includes all parameters, tool versions, and file checksums.

### Q: miRanda is not installed. How do I add it?
**A:** Install via Conda: `conda install -c bioconda miranda`. The tool will auto-detect miRanda when available.

### Q: What is arm switching?
**A:** When the dominant mature miRNA from a precursor changes between conditions (e.g., 5p-dominant in control, 3p-dominant in treatment). This can indicate altered miRNA processing.

---

*sRNAtlas v1.4.0*
