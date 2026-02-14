# Tutorial 2: Your First Small RNA-seq Analysis

A step-by-step walkthrough of the complete sRNAtlas pipeline — from raw FASTQ files to differential expression results.

---

## Overview

|                         |                                                                                                                                         |
| ----------------------- | --------------------------------------------------------------------------------------------------------------------------------------- |
| **Learning objectives** | Run the full sRNAtlas pipeline end-to-end; interpret QC metrics, alignment statistics, and DE results; export publication-ready figures |
| **Prerequisites**       | sRNAtlas installed and running ([Tutorial 1: Installation](01-installation.md))                                                         |
| **Estimated time**      | 20–30 minutes                                                                                                                           |
| **Difficulty**          | Beginner                                                                                                                                |

---

## About the Tutorial Dataset

This tutorial uses a small RNA-seq experiment comparing two conditions (**control** vs. **treatment**), with 2 biological replicates each. The dataset is a downsampled subset of a published *Arabidopsis thaliana* small RNA-seq study, chosen for fast processing times while retaining real biological signal.

| Sample | Condition | Reads |
|--------|-----------|-------|
| `ctrl_rep1.fastq.gz` | control | ~50,000 |
| `ctrl_rep2.fastq.gz` | control | ~50,000 |
| `treat_rep1.fastq.gz` | treatment | ~50,000 |
| `treat_rep2.fastq.gz` | treatment | ~50,000 |

Library preparation: **Illumina TruSeq Small RNA** (3' adapter: `TGGAATTCTCGGGTGCCAAGG`).

### Download the Test Data

Download the tutorial dataset from Zenodo:

```bash
# Option 1: wget
wget https://zenodo.org/records/XXXXXXX/files/sRNAtlas_tutorial_data.tar.gz
tar -xzf sRNAtlas_tutorial_data.tar.gz

# Option 2: curl
curl -OL https://zenodo.org/records/XXXXXXX/files/sRNAtlas_tutorial_data.tar.gz
tar -xzf sRNAtlas_tutorial_data.tar.gz
```

> **Note**: Replace `XXXXXXX` with the actual Zenodo record ID once the dataset is deposited. See [Appendix A: Preparing Your Own Test Dataset](#appendix-a-preparing-your-own-test-dataset) for instructions on creating a test dataset from any published study.

The archive contains:

```
tutorial_data/
├── ctrl_rep1.fastq.gz
├── ctrl_rep2.fastq.gz
├── treat_rep1.fastq.gz
├── treat_rep2.fastq.gz
└── metadata.csv
```

The `metadata.csv` file:

```csv
sample,condition,batch
ctrl_rep1,control,1
ctrl_rep2,control,1
treat_rep1,treatment,1
treat_rep2,treatment,1
```

---

## Step 1: Launch sRNAtlas

Start the application from your terminal:

```bash
conda activate srnatlas
cd sRNAtlas
streamlit run app/main.py
```

Your browser should open automatically to `http://localhost:8501`. You will see the **Home** page with the sRNAtlas logo and navigation sidebar on the left.

> **Troubleshooting**: If the browser doesn't open, navigate to `http://localhost:8501` manually. If port 8501 is in use, Streamlit will suggest an alternative port in the terminal output.

---

## Step 2: Create a Project (1 min)

A project in sRNAtlas organizes all your files, settings, and results into a single session that can be saved and restored later.

1. Click **Project** in the left sidebar
2. In the **Project Settings** tab:
   - **Project name**: Enter `tutorial_ath_srna` (or any descriptive name)
   - **Organism**: Select `Arabidopsis thaliana`
3. Click **Create Project**

A green success message confirms the project was created. All subsequent analysis results will be stored within this project.

> **Tip**: Use a naming convention that includes organism and date, such as `srna-ath_20260212`. This matches the naming convention you'll see in saved `.srna` project files.

---

## Step 3: Upload Data (2 min)

1. Stay in the **Project** module
2. Switch to the **Data Upload** tab
3. **Drag and drop** (or click to browse) the four `.fastq.gz` files:
   - `ctrl_rep1.fastq.gz`
   - `ctrl_rep2.fastq.gz`
   - `treat_rep1.fastq.gz`
   - `treat_rep2.fastq.gz`
4. Click **Add Files to Project**

You should see all four files listed with their sizes. Files up to 2 GB are supported; gzip compression is recommended.

> **Tip**: If you're working on a remote server, the Streamlit upload widget has a default file size limit. For very large files, you can instead place them directly in the project's data directory and sRNAtlas will detect them.

---

## Step 4: Quality Control (2 min)

Pre-alignment QC helps you understand your data *before* any processing. This step evaluates read quality, length distribution, and potential contamination.

1. Click **Quality Control** in the sidebar
2. In the **Upload Data** tab, select **Use files from Project**
3. Click **Run QC Analysis**

Once complete, explore the result tabs:

### 4a. Read Statistics

This tab shows a summary table for each sample:

| Metric | What to expect | What to check |
|--------|---------------|---------------|
| Total reads | ~50,000 (for tutorial data) | Very low count may indicate a failed library |
| Mean quality | >28 (Phred) | Values <20 suggest degraded quality |
| Mean length | Depends on adapter trimming status | Raw reads will be the full sequencing length (e.g., 50 nt) |
| GC content | ~50% (varies by organism) | Extreme values may indicate contamination |

### 4b. Size Distribution

This is the **most informative QC plot for small RNA-seq**. Before trimming, most reads will appear at the full sequencing read length (e.g., 50 nt), because adapters haven't been removed yet. After trimming (Step 5), you will re-examine this distribution and expect to see characteristic peaks.

### 4c. Contamination Check

Evaluates adapter content and rRNA contamination. For a properly prepared small RNA library:

- **High adapter content is expected and normal** — unlike in standard RNA-seq, small RNA inserts (18–25 nt) are shorter than the read length (50 nt), so adapter sequence is read through. This is a *good sign*.
- **rRNA content** should be low (<10%). High rRNA suggests ribosomal RNA depletion was insufficient.

> ✅ **Checkpoint**: Do all four samples have similar read counts and quality scores? If one sample looks very different, take note — it may behave as an outlier later in DE analysis.

---

## Step 5: Adapter Trimming (2 min)

Adapter trimming is critical for small RNA-seq. Because mature miRNAs (~22 nt) are shorter than the sequencing read length, each read contains the insert followed by adapter sequence. Cutadapt removes this adapter and filters reads by length.

1. Click **Trimming** in the sidebar
2. In the **Input Files** tab, select **Use files from Project**
3. In the **Trimming Settings** tab, configure:

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| **Adapter preset** | `Illumina TruSeq Small RNA (RA3)` | Matches our library prep; adapter = `TGGAATTCTCGGGTGCCAAGG` |
| **Min length** | `18` | Shortest expected miRNA; anything shorter is likely degradation |
| **Max length** | `35` | Captures miRNAs (18–25 nt) and piRNAs/siRNAs (up to ~32 nt) |
| **Quality cutoff** | `20` | Trim bases with Phred quality <20 from the 3' end |
| **Error rate** | `0.1` | Allow 10% mismatches in adapter matching (Cutadapt default) |

4. Click **Run Trimming**

### Interpreting Trimming Results

After trimming completes, review the statistics table:

| Metric | Good | Acceptable | Investigate |
|--------|------|------------|-------------|
| **Reads passing filters** | >70% | 50–70% | <50% |
| **Reads with adapter** | >80% | 60–80% | <60% |
| **Reads too short (<18 nt)** | <15% | 15–30% | >30% |
| **Reads too long (>35 nt)** | <10% | 10–20% | >20% |

> ⚠️ **Common problem**: If "Reads with adapter" is very low (<20%), you likely selected the **wrong adapter preset**. Go back and try a different preset, or check your library preparation protocol documentation.

> ⚠️ **Common problem**: If "Reads too short" is very high (>50%), this may indicate severe degradation or adapter dimers. Check your gel images from library prep if available.

> ✅ **Checkpoint**: You should see >70% of reads passing filters. If using the tutorial dataset, expect approximately 75–85% pass rate.

---

## Step 6: Set Up Reference Database (2 min)

Before aligning reads, you need a reference database of known small RNA sequences. sRNAtlas supports miRBase (miRNAs) and RNAcentral (broader ncRNA types).

1. Click **Databases** in the sidebar
2. In the **miRBase** tab:
   - **Organism**: Select `Arabidopsis thaliana` (species code: `ath`)
   - **Sequence type**: Select `mature` (for mature miRNA sequences — the biologically active form)
   - Click **Download from miRBase**
3. Wait for the download to complete (typically a few seconds)
4. You should see a confirmation message showing the number of sequences downloaded (e.g., "Downloaded 428 mature miRNA sequences for ath")

### Build the Bowtie Index

After downloading, you need to build a Bowtie index for efficient alignment:

1. In the same **Databases** module, click **Build Bowtie Index**
2. Wait for index building to complete (1–2 minutes)
3. A success message confirms the index is ready

The index is built using `bowtie-build` and creates the `.ebwt` files that Bowtie needs for rapid alignment.

> **Why Bowtie v1 (not Bowtie2)?** Bowtie v1 is specifically designed for short reads with few mismatches — exactly the characteristics of small RNA-seq data. Bowtie2 is optimized for longer reads and uses a different alignment algorithm that is less suitable for 18–25 nt sequences.

> **Advanced**: For a more comprehensive annotation, you can also download sequences from **RNAcentral** (in the **RNAcentral** tab). This provides tRNAs, rRNAs, snoRNAs, and other ncRNA types beyond miRNAs. For this tutorial, miRBase alone is sufficient.

---

## Step 7: Align Reads (2 min)

Alignment maps your trimmed reads to the reference database to identify which known miRNAs are present in each sample.

1. Click **Alignment** in the sidebar
2. Verify the **Reference index** shows the miRBase index you just built
3. In the **Run Alignment** tab:
   - Select **Use trimmed files** (from Step 5)
   - Review the default alignment parameters:

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `-v 1` | 1 mismatch | Allow 0 or 1 mismatch across the entire read. Accounts for SNPs and sequencing errors |
| `-k 10` | Report up to 10 | Report up to 10 alignments per read (for multi-mappers) |
| `--best` | Enabled | Report alignments in best-first order |
| `--strata` | Enabled | Only report alignments in the best mismatch stratum (e.g., if a 0-mismatch alignment exists, don't report 1-mismatch ones) |

4. Click **Start Alignment**

### Interpreting Alignment Results

After alignment, review the per-sample statistics:

| Metric | Good | Acceptable | Investigate |
|--------|------|------------|-------------|
| **Overall alignment rate** | >70% | 30–70% | <30% |
| **Unique mapping rate** | >50% | 30–50% | <30% |

> ⚠️ **Low alignment rate?** The most common causes are:
> - **Wrong reference**: You selected the wrong organism in miRBase
> - **Trimming issue**: Adapters were not fully removed (re-check Step 5)
> - **Contamination**: Many reads are rRNA or other non-miRNA species (try aligning to RNAcentral)
> - **Novel biology**: Many reads are bona fide small RNAs not yet in miRBase (explore the Novel miRNA module later)

> **What happens to unaligned reads?** sRNAtlas automatically saves reads that did not align to the reference. These can be used later in the **Novel miRNA Discovery** module (Tutorial 6) to identify unannotated small RNAs.

> ✅ **Checkpoint**: Each sample should show a similar alignment rate. Large differences between samples may indicate technical problems with specific libraries.

---

## Step 8: Post-Alignment QC (2 min)

Post-alignment QC evaluates the quality of your alignments and provides insight into the RNA composition of your libraries.

1. Click **Post-Align QC** in the sidebar
2. Click **Run Post-Alignment QC**

### Key Metrics to Review

#### Mapping Statistics
A summary of mapped vs. unmapped reads, unique vs. multi-mapped reads, and overall alignment rates — consolidating the information from Step 7.

#### Strand Bias
For miRNA-seq aligned to mature sequences, you expect **predominantly forward-strand** alignments (strand bias close to 0 or 1, depending on the metric). A roughly equal forward/reverse split would suggest non-specific mapping.

#### 5' Nucleotide Composition
This is a signature metric for miRNA biology:
- **Plant miRNAs**: Expect a strong **5' U (uridine)** bias, especially for 21-nt reads. This reflects the preference of AGO1 for 5'-U miRNAs.
- **Animal miRNAs**: Also show a 5' U bias for most miRNA families.
- If you see 5' A enrichment in 24-nt reads (in plants), these are likely **heterochromatic siRNAs** loaded into AGO4.

#### Read Length Distribution (Post-Alignment)
Now that reads are aligned to known miRNAs, the length distribution should be much more informative than pre-trimming:
- **Peak at 21 nt**: miRNAs (the dominant class in most sRNA-seq experiments)
- **Peak at 24 nt** (plants): hc-siRNAs (if you aligned to a broader reference that includes siRNAs)

#### RNA Type Composition
A pie chart or bar chart showing the relative abundance of each miRNA across samples. The most abundant miRNAs (e.g., ath-miR156, ath-miR159, ath-miR166 in *Arabidopsis*) should dominate.

> ✅ **Checkpoint**: All samples should show consistent 5' nucleotide profiles and length distributions. If one sample deviates, flag it as a potential outlier.

---

## Step 9: Read Counting (1 min)

The counting module quantifies how many reads align to each reference feature (miRNA), producing a **count matrix** that is the input for differential expression analysis.

1. Click **Counting** in the sidebar
2. Select the BAM files from the alignment step (they should auto-populate)
3. Click **Generate Count Matrix**

### Understanding the Count Matrix

The output is a table with miRNA names as rows and samples as columns:

| | ctrl_rep1 | ctrl_rep2 | treat_rep1 | treat_rep2 |
|---|---|---|---|---|
| ath-miR156a-5p | 1,250 | 1,180 | 580 | 620 |
| ath-miR159a | 3,400 | 3,250 | 3,100 | 3,350 |
| ath-miR166a-3p | 2,800 | 2,650 | 2,900 | 2,750 |
| ath-miR168a-5p | 890 | 920 | 1,850 | 1,780 |
| ... | ... | ... | ... | ... |

> **Note**: These numbers are illustrative. Your actual counts will depend on the tutorial dataset.

Key observations to check:
- **Total counts per sample** should be roughly similar (within 2-fold). Large differences suggest library size issues.
- **Number of detected miRNAs** (rows with >0 counts) should be consistent across samples.

The count matrix is automatically saved in your project and will be the input for the next step.

---

## Step 10: Differential Expression Analysis (5 min)

This is the central analysis step — identifying which miRNAs are significantly up- or down-regulated between conditions. sRNAtlas uses **pyDESeq2**, a Python implementation of the widely-used DESeq2 algorithm.

1. Click **DE Analysis** in the sidebar

### 10a. Set Up the Comparison

1. **Upload or confirm the count matrix**: The count matrix from Step 9 should auto-populate
2. **Upload sample metadata**: Upload the `metadata.csv` file, or enter the metadata manually:
   
   | SampleID | condition |
   |----------|-----------|
   | ctrl_rep1 | control |
   | ctrl_rep2 | control |
   | treat_rep1 | treatment |
   | treat_rep2 | treatment |

3. **Design factor**: Select `condition` as the main factor
4. **Set up comparison**: 
   - **Numerator** (test): `treatment`
   - **Denominator** (reference): `control`
   
   This means positive log2FC = higher in treatment, negative log2FC = higher in control.

5. **Filtering settings**:
   - **Minimum count filter**: `10` (remove miRNAs with fewer than 10 total reads across all samples — these are too lowly expressed for reliable statistics)
   - **FDR threshold**: `0.05`
   - **log2FC threshold**: `0.585` (equivalent to a 1.5-fold change)

6. Click **Run DE Analysis**

### 10b. Interpreting the Results

After the analysis completes, explore each visualization tab:

#### PCA Plot
Principal Component Analysis shows the overall relationship between samples:
- **Replicates should cluster together** (control samples close to each other, treatment samples close to each other)
- **Conditions should separate** along PC1 (the axis explaining the most variance)
- If a replicate clusters with the wrong group, it may be mislabeled or an outlier

#### Volcano Plot
The volcano plot displays statistical significance (-log10 adjusted p-value) against biological effect size (log2 fold change):
- **Upper-left quadrant**: Significantly down-regulated in treatment (high significance, negative fold change)
- **Upper-right quadrant**: Significantly up-regulated in treatment (high significance, positive fold change)
- **Colored points**: Pass both the FDR and fold-change thresholds
- **Gray points**: Not significant

#### MA Plot
Mean expression (average across all samples) vs. log2 fold change:
- Significant miRNAs are highlighted
- Useful for spotting fold-change compression at low expression levels
- Well-behaved data shows a symmetric cloud centered at log2FC = 0

#### Heatmap
Hierarchical clustering of the top differentially expressed miRNAs:
- Rows = miRNAs, columns = samples
- Color intensity = normalized expression (z-score)
- Samples should cluster by condition (control together, treatment together)
- miRNAs with similar expression patterns cluster together, potentially revealing co-regulated groups

#### DE Results Table
A downloadable table with all results:

| Column | Description |
|--------|-------------|
| `baseMean` | Mean normalized count across all samples |
| `log2FoldChange` | Effect size (positive = up in treatment) |
| `lfcSE` | Standard error of the log2FC estimate |
| `stat` | Wald test statistic |
| `pvalue` | Raw p-value |
| `padj` | Adjusted p-value (Benjamini-Hochberg FDR) |

**How to read the results**: A miRNA with `padj < 0.05` and `|log2FoldChange| > 0.585` is considered significantly differentially expressed. For example:

| miRNA | baseMean | log2FC | padj | Interpretation |
|-------|----------|--------|------|---------------|
| ath-miR156a-5p | 907 | -1.15 | 0.003 | Significantly **down** in treatment |
| ath-miR168a-5p | 1,360 | 1.02 | 0.008 | Significantly **up** in treatment |
| ath-miR159a | 3,275 | -0.06 | 0.82 | **Not significant** (no change) |

### 10c. Download Results

Click the download buttons to export:
- **CSV table**: Full DE results for all miRNAs
- **Figures**: PNG, SVG, PDF, or HTML versions of each plot

> **Publication tip**: The interactive HTML figures are excellent for supplementary materials. For main-text figures, use the SVG or PDF exports (requires `kaleido`: `pip install kaleido`).

> ✅ **Checkpoint**: You should see at least a few significantly DE miRNAs in the tutorial dataset. If you see zero significant results, check that: (a) your metadata correctly matches sample names, (b) you haven't accidentally swapped the condition labels, (c) there are enough reads per miRNA after filtering.

---

## Step 11: Export a Report (1 min)

Generate a comprehensive HTML report summarizing the entire analysis.

1. Click **Reports** in the sidebar
2. Select which sections to include:
   - ✅ QC summary
   - ✅ Trimming statistics
   - ✅ Alignment results
   - ✅ Count matrix summary
   - ✅ DE analysis results
   - ✅ Figures
3. Click **Generate Report**
4. Click **Download** to save the HTML report

The report is a self-contained HTML file that includes all metrics, tables, and figures — suitable for sharing with collaborators or archiving.

---

## Step 12: Save Your Project (1 min)

Save your complete analysis session so you can return to it later.

1. Click **Project** in the sidebar
2. Switch to the **Save/Load** tab (or equivalent)
3. Click **Save Project**
4. The project is saved as a `.srna` file (JSON format) in the `output/projects/` directory

The saved project includes:
- Count matrices and annotations
- DE analysis results
- Enrichment results
- All analysis settings
- Sample metadata

To restore a project later, use the **Load Project** option in the same tab.

---

## Summary

Congratulations! You have completed a full small RNA-seq analysis using sRNAtlas. Here is what you accomplished:

| Step | Module | Output |
|------|--------|--------|
| 1 | **Project** | Created organized project |
| 2 | **Project** | Uploaded 4 FASTQ files |
| 3 | **Quality Control** | Assessed raw read quality |
| 4 | **Trimming** | Removed adapters, filtered by length |
| 5 | **Databases** | Downloaded miRBase, built Bowtie index |
| 6 | **Alignment** | Mapped reads to miRNA reference |
| 7 | **Post-Align QC** | Verified alignment quality |
| 8 | **Counting** | Generated count matrix |
| 9 | **DE Analysis** | Identified differentially expressed miRNAs |
| 10 | **Reports** | Exported HTML report |
| 11 | **Project** | Saved the analysis session |

---

## What's Next?

Now that you've completed the core pipeline, explore the advanced modules:

| Tutorial | Module | What you'll learn |
|----------|--------|-------------------|
| [Tutorial 3: Quality Control Deep Dive](03-quality-control.md) | QC | Interpret size distributions, diagnose bad libraries |
| [Tutorial 4: Databases & Alignment](04-databases-alignment.md) | Databases + Alignment | RNAcentral, custom references, multi-database strategies |
| [Tutorial 5: Differential Expression](05-differential-expression.md) | DE Analysis | Multi-group comparisons, batch effects, advanced visualization |
| [Tutorial 6: Novel miRNA & isomiR](06-novel-mirna-isomir.md) | Novel miRNA + isomiR | Discover unannotated miRNAs, analyze sequence variants |
| [Tutorial 7: Targets & Enrichment](07-targets-enrichment.md) | Targets + GO/Pathway | Predict miRNA targets, run GO and KEGG enrichment |
| [Tutorial 8: Batch & Reproducibility](08-batch-reproducibility.md) | Batch + Reports | Automate full pipelines, manage projects |

---

## Quick Reference Card

### Pipeline at a Glance

```
FASTQ → QC → Trimming → Database → Alignment → Post-Align QC → Counting → DE → Report
```

### Recommended Default Parameters

| Parameter | miRNA-focused | Discovery mode |
|-----------|--------------|----------------|
| Min length | 18 nt | 15 nt |
| Max length | 25 nt | 40 nt |
| Adapter | Illumina TruSeq RA3 | Check your kit |
| Mismatches (Bowtie -v) | 0–1 | 2 |
| Multi-mapper (-k) | 10 | 10 |
| DE FDR threshold | 0.05 | 0.05 |
| DE log2FC threshold | 0.585 (1.5×) | 1.0 (2×) |

### Quality Benchmarks

| Metric | ✅ Good | ⚠️ Acceptable | ❌ Investigate |
|--------|---------|---------------|---------------|
| Trimming pass rate | >70% | 50–70% | <50% |
| Adapter found | >80% | 60–80% | <60% |
| Alignment rate | >70% | 30–70% | <30% |
| Unique mapping | >50% | 30–50% | <30% |

---

## Appendix A: Preparing Your Own Test Dataset

If you want to create a tutorial dataset from a public study instead of using the pre-bundled data:

### 1. Find a Dataset

Search the [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) for small RNA-seq data:
- Search term: `"small RNA-seq" AND "Arabidopsis thaliana"` (or your organism)
- Filter by: Library Strategy = miRNA-Seq
- Look for studies with ≥2 conditions and ≥2 replicates per condition

### 2. Download with SRA Toolkit

```bash
# Install SRA Toolkit
conda install -c bioconda sra-tools

# Download (replace SRR numbers with your chosen accessions)
fastq-dump --gzip SRR1234567
fastq-dump --gzip SRR1234568
fastq-dump --gzip SRR1234569
fastq-dump --gzip SRR1234570
```

### 3. Downsample for Tutorial Speed

```bash
# Downsample to 50,000 reads using seqtk
conda install -c bioconda seqtk

seqtk sample -s 42 SRR1234567.fastq.gz 50000 | gzip > ctrl_rep1.fastq.gz
seqtk sample -s 42 SRR1234568.fastq.gz 50000 | gzip > ctrl_rep2.fastq.gz
seqtk sample -s 42 SRR1234569.fastq.gz 50000 | gzip > treat_rep1.fastq.gz
seqtk sample -s 42 SRR1234570.fastq.gz 50000 | gzip > treat_rep2.fastq.gz
```

### 4. Create Metadata

```bash
cat > metadata.csv << EOF
sample,condition,batch
ctrl_rep1,control,1
ctrl_rep2,control,1
treat_rep1,treatment,1
treat_rep2,treatment,1
EOF
```

---

## Appendix B: Troubleshooting

| Problem | Likely cause | Solution |
|---------|-------------|----------|
| App won't start | Port conflict or missing dependency | Check `pip install -r requirements.txt`; try `streamlit run app/main.py --server.port 8502` |
| 0 reads after trimming | Wrong adapter selected | Try a different adapter preset; check your library prep kit documentation |
| 0% alignment rate | Wrong organism or un-trimmed adapters | Verify the organism code matches your data; re-run trimming |
| No BAM files found | Alignment didn't complete | Check the Alignment module for error messages; verify Bowtie is installed (`bowtie --version`) |
| Count matrix is empty | BAM files have no mapped reads | Re-run alignment with different parameters (try `-v 2` for more mismatches) |
| No DE results | Sample names don't match metadata | Ensure `SampleID` in metadata exactly matches the sample names in the count matrix columns |
| `bowtie: not found` | Not in PATH | Run `conda install -c bioconda bowtie` or add Bowtie to your PATH |
| `pysam` import error | Missing htslib dependency | Run `conda install -c bioconda pysam` (conda handles the C dependencies better than pip) |
| Memory error with large files | Insufficient RAM | Process fewer samples at once, or increase available memory |
| Figures not exporting as PNG/SVG/PDF | kaleido not installed | Run `pip install kaleido` |

---

## Appendix C: Glossary

| Term | Definition |
|------|-----------|
| **Adapter** | Short DNA sequence ligated during library preparation; must be removed bioinformatically |
| **BAM** | Binary Alignment Map — compressed format for storing aligned reads |
| **Bowtie index** | Pre-built data structure enabling rapid alignment of reads to a reference |
| **Count matrix** | Table of read counts per feature (miRNA) per sample; input for DE analysis |
| **DE (Differential Expression)** | Statistical identification of features with significantly different abundance between conditions |
| **FDR** | False Discovery Rate — the expected proportion of false positives among significant results (Benjamini-Hochberg correction) |
| **FASTQ** | Text-based format for storing sequencing reads with quality scores |
| **log2FC** | Log2 fold change — a measure of effect size. log2FC of 1 = 2-fold increase; log2FC of -1 = 2-fold decrease |
| **miRBase** | The primary database of published miRNA sequences and annotations |
| **Multi-mapper** | A read that aligns equally well to multiple reference sequences |
| **Normalized counts** | Read counts adjusted for library size differences between samples |
| **padj** | Adjusted p-value (FDR-corrected); the primary significance metric in DE analysis |
| **Phred score** | Logarithmic quality score for each sequenced base. Q20 = 99% accuracy; Q30 = 99.9% accuracy |
| **pyDESeq2** | Python implementation of the DESeq2 statistical framework for count-based DE analysis |
| **RNAcentral** | Comprehensive database of non-coding RNA sequences from multiple source databases |

---

*Tutorial written for sRNAtlas v0.1.0 (BETA). Last updated: 2026-02-12.*
