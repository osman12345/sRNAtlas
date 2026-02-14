# sRNAtlas: an interactive platform for end-to-end small RNA sequencing analysis

**Osman Radwan**^1,\*

^1 Centre for European Masters Programmes (CEMP), Riga, Latvia

\* Corresponding author. E-mail: [corresponding email]

**Associate Editor:** [Editor name]

Received on XXXXX; revised on XXXXX; accepted on XXXXX

---

## Abstract

**Summary:** Small RNA sequencing (sRNA-seq) is the primary method for profiling microRNAs, siRNAs, piRNAs, and other non-coding RNA fragments, yet existing analysis tools either require substantial command-line expertise or cover only a subset of the analytical workflow. Here we present sRNAtlas, an open-source, browser-based platform that provides a complete, interactive pipeline from raw FASTQ files to differential expression results, target prediction, and functional enrichment. Built with Streamlit and integrating established tools—Cutadapt, Bowtie, SAMtools, and pyDESeq2—sRNAtlas offers adapter trimming with multi-platform presets, automated reference management for miRBase and RNAcentral, organism-aware quality control with diagnostic size-distribution interpretation, novel miRNA discovery, isomiR detection, miRNA target prediction via psRNATarget, and GO/KEGG pathway enrichment via g:Profiler. All parameters are transparent and adjustable, and full projects can be saved, reloaded, and batch-processed. sRNAtlas is designed to make rigorous small RNA analysis accessible to bench scientists while remaining fully configurable for bioinformaticians.

**Availability and implementation:** sRNAtlas is freely available at https://github.com/osman12345/sRNAtlas. It runs on Linux, macOS, and Windows via Conda, pip, or Docker. Documentation and tutorials are provided at the repository.

**Contact:** [corresponding email]

**Supplementary information:** Supplementary data are available at *Bioinformatics* online.

---

## 1 Introduction

Small non-coding RNAs (sRNAs) regulate gene expression, chromatin state, and genome defence across eukaryotes (Borges and Martienssen, 2015). Small RNA sequencing (sRNA-seq) is the standard method for their genome-wide profiling, but the analysis presents unique computational challenges. Inserts of 18–32 nt are shorter than typical read lengths, making adapter trimming critical; read-length distributions carry direct biological information about RNA class composition; and multi-mapping is pervasive due to conserved miRNA families (Axtell, 2013; Kozomara *et al.*, 2019).

Several tools address parts of this workflow. sRNAbench/sRNAtoolbox provides a comprehensive web server but requires uploading data to an external host and is limited to pre-defined organisms (Aparicio-Puerta *et al.*, 2022). The nf-core/smrnaseq Nextflow pipeline offers reproducible processing but demands familiarity with Nextflow and the command line (Ewels *et al.*, 2020). miRDeep2 focuses on novel miRNA prediction but does not provide trimming, QC, or differential expression (Friedländer *et al.*, 2012). sRNAflow integrates several steps into a Snakemake pipeline but lacks an interactive interface and functional enrichment modules (Zanini *et al.*, 2024). Galaxy-based workflows offer interactivity but require server administration or reliance on public instances with upload constraints (Afgan *et al.*, 2018).

A persistent gap exists for researchers who need a *complete*, *interactive*, *locally installable* platform that covers the entire sRNA-seq workflow—from raw reads to biological interpretation—without requiring command-line skills or external data uploads. sRNAtlas addresses this gap.

## 2 Implementation

sRNAtlas is implemented in Python and deployed as a Streamlit web application. The interface is organized into 14 interconnected modules (Fig. 1) that guide users through a logical analytical workflow while allowing entry at any stage with pre-existing data (e.g., count matrices). All external tool invocations, intermediate files, and parameter choices are logged for reproducibility.

### 2.1 Adapter trimming and quality control

The **Trimming** module wraps Cutadapt (Martin, 2011) and provides pre-configured adapter presets for Illumina TruSeq Small RNA, NEBNext, QIAseq, Ion Torrent, BGI/MGI, Nanopore, and Takara SMARTer kits. Users select a preset or supply a custom sequence; minimum/maximum length filters (default: 18–35 nt), quality thresholds (default: Phred ≥ 20), and other Cutadapt parameters are exposed and adjustable.

The **Quality Control** module performs FASTQ-level analysis by sampling up to 100 000 reads per file, computing per-read length, quality, and GC content distributions, and running contamination scanning for adapter remnants (12-bp signatures for nine common platforms), poly-A/T content, and overrepresented 6-mers. Read-length distributions are overlaid with shaded RNA-type size ranges (miRNA: 18–25 nt; siRNA: 20–24 nt; piRNA: 24–32 nt; tRF: 14–40 nt; snoRNA: 60–300 nt; snRNA: 100–300 nt), with automatic peak detection and organism-aware interpretation. Plant-specific RNA types—tasiRNA, phasiRNA, natsiRNA, and hc-siRNA—are included in the size-range definitions.

### 2.2 Reference database management

The **Databases** module automates retrieval and indexing of reference sequences. For miRBase (Kozomara *et al.*, 2019), mature or hairpin FASTA files are fetched and filtered by three-letter organism codes (e.g., `ath`, `mtr`, `hsa`), with 39 supported organisms. For RNAcentral (The RNAcentral Consortium, 2021), species-specific sequences for selected RNA types (miRNA, tRNA, rRNA, snoRNA, snRNA, piRNA, lncRNA) are downloaded. Users can also supply custom FASTA files. All references are indexed with `bowtie-build` (Langmead *et al.*, 2009) within the interface; resulting indices are stored and automatically populated in downstream modules.

### 2.3 Alignment

The **Alignment** module invokes Bowtie v1 (Langmead *et al.*, 2009) with parameters optimized for short ncRNAs: `-v 1` (one mismatch across the full read), `-k 10` (up to 10 alignments per read), and `--best --strata` (report only the best-scoring stratum). These defaults reflect the short read length, the prevalence of near-identical miRNA family members, and the importance of capturing multi-mappers for accurate quantification. All parameters—mismatches, maximum alignments, multi-mapper suppression threshold, strand handling—are adjustable. SAM output is converted to sorted, indexed BAM via SAMtools (Danecek *et al.*, 2021). Unaligned reads are written to a separate FASTQ file for downstream novel miRNA discovery or hierarchical alignment to secondary references.

### 2.4 Post-alignment quality control

The **Post-Align QC** module parses BAM files to compute mapping rate, unique vs. multi-mapped read fractions, strand bias, MAPQ distribution, and 5ʹ-nucleotide composition. The 5ʹ-nucleotide profile serves as a biological validation: in plants, 5ʹ-U enrichment among 21-nt reads confirms AGO1-loaded miRNAs, while 5ʹ-A enrichment among 24-nt reads indicates AGO4-loaded hc-siRNAs (Mi *et al.*, 2008). Post-alignment size distributions are compared to pre-alignment profiles to quantify the fraction of each size class captured by the reference.

### 2.5 Read counting

The **Counting** module extracts per-feature read counts from BAM files, generates a sample × feature count matrix, and annotates features with RNA type metadata from the reference FASTA headers. Low-abundance filtering (configurable minimum counts per sample, minimum samples detected, and CPM threshold) is applied before export.

### 2.6 Differential expression analysis

The **DE Analysis** module implements differential expression testing via pyDESeq2 (Muzellec *et al.*, 2023), providing a pure-Python reimplementation of the DESeq2 model (Love *et al.*, 2014) that eliminates the R dependency and enables deployment on Streamlit Cloud. Users configure the design factor, define pairwise comparisons (reference vs. test), and set FDR and log₂ fold-change thresholds (defaults: adjusted *P* < 0.05, |log₂FC| > 0.585). A global FDR correction option adjusts *P*-values across all comparisons simultaneously. Results are presented as interactive volcano plots, MA plots, heatmaps of top DE features, and PCA of normalized counts. All plots are generated with Plotly and can be exported as publication-quality PNG/SVG files.

### 2.7 Novel miRNA discovery and isomiR analysis

The **Novel miRNA** module uses unaligned reads and a genome reference to identify candidate novel miRNA loci based on structural and expression criteria. The **isomiR** module detects sequence variants of known miRNAs—5ʹ and 3ʹ length variants, internal SNPs, and non-templated additions—by comparing aligned reads to canonical mature sequences from miRBase.

### 2.8 Target prediction and functional enrichment

The **Target Prediction** module supports two modes: (i) remote submission to the psRNATarget server (Dai *et al.*, 2018), which uses an established scoring scheme for plant miRNA–target complementarity, and (ii) a local seed-matching algorithm for rapid preliminary screening. Differentially expressed miRNAs can be fed directly from the DE results.

The **GO/Pathway Enrichment** module calls g:Profiler (Kolberg *et al.*, 2023) to perform Gene Ontology (BP, MF, CC) and KEGG pathway enrichment on predicted target gene lists, with configurable background gene sets and FDR thresholds.

### 2.9 Batch processing, project management, and reporting

The **Batch** module enables automated execution of the full pipeline or individual stages across all samples, with job queuing and progress tracking. The **Project** module manages sample metadata, experimental design, and file associations; entire analysis sessions—including all intermediate results—can be saved to disk and reloaded. The **Reports** module generates HTML summary reports containing QC metrics, DE results, enrichment tables, and all figures, and supports bulk export of data tables and visualizations.

## 3 Comparison with existing tools

Table 1 compares sRNAtlas with five widely used sRNA-seq analysis platforms across key functional and usability dimensions.

**Table 1.** Feature comparison of sRNA-seq analysis platforms

| Feature | sRNAtlas | sRNAtoolbox | nf-core/smrnaseq | miRDeep2 | sRNAflow | Galaxy |
|---|---|---|---|---|---|---|
| Interface | Web GUI (local) | Web GUI (remote) | CLI | CLI | CLI | Web GUI (remote/local) |
| Local installation | ✓ | — | ✓ | ✓ | ✓ | Complex |
| No data upload required | ✓ | — | ✓ | ✓ | ✓ | Depends |
| Adapter trimming | ✓ (Cutadapt, 9 presets) | ✓ | ✓ | — | ✓ | ✓ |
| QC with size interpretation | ✓ (organism-aware) | Partial | ✓ | — | Partial | ✓ |
| miRBase integration | ✓ (auto-fetch) | ✓ | ✓ | ✓ | ✓ | Manual |
| RNAcentral integration | ✓ (auto-fetch) | — | — | — | — | Manual |
| Alignment | Bowtie v1 | Bowtie v1 | Bowtie v1 | Bowtie v1 | Bowtie v1/v2 | Configurable |
| Differential expression | ✓ (pyDESeq2) | ✓ (edgeR) | — | — | ✓ (DESeq2) | ✓ |
| Novel miRNA discovery | ✓ | ✓ | ✓ (miRDeep2) | ✓ | ✓ | ✓ |
| isomiR detection | ✓ | ✓ | — | — | — | — |
| Target prediction | ✓ (psRNATarget) | ✓ | — | — | — | ✓ |
| GO/Pathway enrichment | ✓ (g:Profiler) | Partial | — | — | — | ✓ |
| Batch processing | ✓ | ✓ | ✓ | — | ✓ | ✓ |
| Project save/reload | ✓ | — | — | — | — | ✓ |
| Plant sRNA types | ✓ (tasiRNA, phasiRNA, hc-siRNA) | Partial | — | — | — | — |
| R dependency | None | Yes | Yes | Yes | Yes | Depends |
| Docker support | ✓ | — | ✓ | — | ✓ | ✓ |

Key differentiators of sRNAtlas include: (i) a complete end-to-end workflow from FASTQ to enrichment in a single local application with no R dependency; (ii) integrated RNAcentral support for non-miRNA small RNA classes (tRFs, rsRFs, snoRNA fragments); (iii) organism-aware QC diagnostics with plant-specific small RNA type definitions; (iv) transparent, adjustable parameters at every stage; and (v) project persistence with full session save/reload.

## 4 Conclusion

sRNAtlas provides researchers with a comprehensive, locally deployable platform for small RNA-seq analysis that bridges the gap between command-line pipelines and limited web servers. By integrating established algorithms into a transparent, interactive workflow—and extending support to non-miRNA small RNA classes, plant-specific biology, and functional interpretation—sRNAtlas aims to make rigorous sRNA-seq analysis accessible to the broader life sciences community.

## Acknowledgements

The author thanks the Streamlit, pyDESeq2, Cutadapt, Bowtie, SAMtools, miRBase, RNAcentral, psRNATarget, and g:Profiler development teams for making their tools freely available.

## Funding

This work has been supported by [funding source].

*Conflict of Interest:* none declared.

---

## References

Afgan,E. *et al.* (2018) The Galaxy platform for accessible, reproducible and collaborative biomedical analyses: 2018 update. *Nucleic Acids Res.*, **46**, W537–W544.

Aparicio-Puerta,E. *et al.* (2022) sRNAbench and sRNAtoolbox 2022 update: accurate miRNA and sRNA profiling for model and non-model organisms. *Nucleic Acids Res.*, **50**, W710–W717.

Axtell,M.J. (2013) Classification and comparison of small RNAs from plants. *Annu. Rev. Plant Biol.*, **64**, 137–159.

Borges,F. and Martienssen,R.A. (2015) The expanding world of small RNAs in plants. *Nat. Rev. Mol. Cell Biol.*, **16**, 727–741.

Dai,X. *et al.* (2018) psRNATarget: a plant small RNA target analysis server (2017 release). *Nucleic Acids Res.*, **46**, W49–W54.

Danecek,P. *et al.* (2021) Twelve years of SAMtools and BCFtools. *GigaScience*, **10**, giab008.

Ewels,P.A. *et al.* (2020) The nf-core framework for community-curated bioinformatics pipelines. *Nat. Biotechnol.*, **38**, 276–278.

Friedländer,M.R. *et al.* (2012) miRDeep2 accurately identifies known and hundreds of novel microRNA genes in seven animal clades. *Nucleic Acids Res.*, **40**, 37–52.

Kolberg,L. *et al.* (2023) g:Profiler—interoperable web service for functional enrichment analysis and gene identifier mapping (2023 update). *Nucleic Acids Res.*, **51**, W431–W437.

Kozomara,A. *et al.* (2019) miRBase: from microRNA sequences to function. *Nucleic Acids Res.*, **47**, D155–D162.

Langmead,B. *et al.* (2009) Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. *Genome Biol.*, **10**, R25.

Love,M.I. *et al.* (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol.*, **15**, 550.

Martin,M. (2011) Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal*, **17**, 10–12.

Mi,S. *et al.* (2008) Sorting of small RNAs into Arabidopsis argonaute complexes is directed by the 5ʹ terminal nucleotide. *Cell*, **133**, 116–127.

Muzellec,B. *et al.* (2023) PyDESeq2: a python package for bulk RNA-seq differential expression analysis. *Bioinformatics*, **39**, btad547.

The RNAcentral Consortium (2021) RNAcentral 2021: secondary structure integration, improved sequence search and new member databases. *Nucleic Acids Res.*, **49**, D212–D220.

Zanini,S. *et al.* (2024) sRNAflow: a tool for the analysis of small RNA-seq data. *Bioinform. Adv.*, **4**, vbae003.

---

## Figure Legend

**Fig. 1. sRNAtlas architecture and workflow.** The platform consists of 14 interconnected modules organized into four analytical phases. **Phase 1 — Preprocessing** (blue): Project setup, quality control with organism-aware size-distribution diagnostics, adapter trimming via Cutadapt with multi-platform presets, and reference database management integrating miRBase and RNAcentral. **Phase 2 — Mapping and Quantification** (green): Bowtie v1 alignment with small-RNA-optimized parameters, post-alignment QC including 5ʹ-nucleotide and strand-bias analysis, and count matrix generation with RNA-type annotation. **Phase 3 — Statistical Analysis** (orange): Differential expression via pyDESeq2 with interactive volcano/MA/PCA plots, novel miRNA discovery from unaligned reads, and isomiR variant detection. **Phase 4 — Biological Interpretation** (purple): miRNA target prediction via psRNATarget, GO and KEGG pathway enrichment via g:Profiler, batch automation, project persistence, and HTML report generation. Data flows sequentially through phases (solid arrows) but users may enter at any stage with pre-existing data (dashed arrows). All modules share a persistent session state, enabling iterative refinement without re-running upstream steps. The interface is implemented in Streamlit with Plotly-based interactive visualizations.
