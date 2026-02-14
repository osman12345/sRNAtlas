# Tutorial 8 – Advanced Analysis: Hierarchical Alignment, isomiRs, and Novel miRNA Discovery

In this tutorial you will explore three specialized modules that extend the core sRNA-seq workflow:

1. **Hierarchical alignment** – sequential alignment to multiple references (miRBase → RNAcentral → genome) to classify reads by RNA type.
2. **isomiR detection** – identifying sequence variants of known miRNAs.
3. **Novel miRNA discovery** – screening unaligned reads for candidate novel miRNA loci.

These analyses are optional extensions to the basic pipeline (Tutorials 1–5) and are most relevant when:

- Your organism has complex small RNA biogenesis (plants with multiple siRNA classes, animals with tRNA fragments).
- You want to characterize miRNA isoform diversity.
- You suspect your samples contain unannotated miRNAs.

***

## 8.1 Hierarchical alignment for multi-ncRNA classification

### 8.1.1 Motivation

In a typical sRNA-seq library, reads derive from multiple RNA classes: miRNAs, siRNAs, tRNA fragments (tRFs), rRNA fragments, snoRNAs, and degradation products. A single alignment to miRBase captures only the miRNA fraction, leaving the rest "unaligned."

**Hierarchical alignment** classifies reads progressively:

1. Align to **miRBase** → capture miRNAs.
2. Take unaligned reads → align to **RNAcentral tRNA** → capture tRFs.
3. Take remaining unaligned reads → align to **RNAcentral rRNA** → capture rRNA contaminants.
4. Take final unaligned reads → align to **genome** → capture novel loci or other genomic ncRNAs.

This strategy, analogous to SPORTS1.1's sequential annotation (Shi *et al.*, 2018), provides:

- A comprehensive accounting of where reads originate.
- QC insights (e.g., high rRNA fraction indicates poor rRNA depletion).
- Input for specialized analyses (e.g., tRF-specific differential expression).

### 8.1.2 Workflow setup

#### Step 1: Prepare references in the Databases module

Before starting hierarchical alignment, ensure you have:

1. **Primary reference (miRBase)**: Already indexed (see Tutorial 2).
2. **Secondary references (RNAcentral)**:
   - Go to **Databases → RNAcentral**.
   - Select your organism (e.g., *Medicago truncatula*).
   - Check the RNA types to download:
     - `tRNA` (for tRNA fragments)
     - `rRNA` (for ribosomal RNA contaminants)
     - Optionally: `snoRNA`, `snRNA`, `piRNA`, `lncRNA`.
   - Click **"Download sequences"**.
   - After download completes, click **"Build index"** for each RNA type.

3. **Genome reference (optional)**:
   - If you plan to align final unaligned reads to the genome, upload a genome FASTA (or chromosome subset) in **Databases → Custom References**.
   - Index it with Bowtie.

#### Step 2: Run the first alignment (miRBase)

1. Go to **Alignment** module.
2. Select your trimmed FASTQ files (output from Tutorial 3).
3. Choose the **miRBase reference** for your organism.
4. Set alignment parameters (defaults: `-v 1 -k 10 --best --strata`).
5. **Enable "Write unaligned reads"** – this is critical for hierarchical alignment.
   - Check the box **"Save unaligned reads to FASTQ"**.
   - Unaligned reads will be saved as `<sample>_unaligned.fastq.gz` in the output directory.
6. Click **"Run alignment"**.

After completion, you will have:

- `<sample>.bam` – reads that aligned to miRBase.
- `<sample>_unaligned.fastq.gz` – reads that did not align.

#### Step 3: Run secondary alignment (tRNA)

1. In the **Alignment** module, load the **unaligned FASTQ files** from Step 2.
   - Use the file uploader or browse to `output/alignment/<sample>_unaligned.fastq.gz`.
2. Select the **RNAcentral tRNA reference**.
3. Keep alignment parameters consistent (same `-v`, `-k` settings).
4. Again, **enable "Write unaligned reads"**.
5. Click **"Run alignment"**.

This produces:

- `<sample>_tRNA.bam` – tRNA-aligned reads.
- `<sample>_tRNA_unaligned.fastq.gz` – still-unaligned reads.

#### Step 4: Run tertiary alignment (rRNA)

Repeat Step 3, but:

- Input: `<sample>_tRNA_unaligned.fastq.gz`.
- Reference: **RNAcentral rRNA**.
- Output: `<sample>_rRNA.bam` and `<sample>_rRNA_unaligned.fastq.gz`.

#### Step 5 (optional): Genome alignment

If you want to capture novel miRNAs or unannotated genomic small RNAs:

- Input: `<sample>_rRNA_unaligned.fastq.gz`.
- Reference: **Genome index**.
- Output: `<sample>_genome.bam` and `<sample>_genome_unaligned.fastq.gz`.

At this point, `<sample>_genome_unaligned.fastq.gz` contains reads that did not map to any of your references—these may be:

- Adapter dimers or quality artifacts (if they survived trimming).
- Reads from contaminating organisms.
- Highly degraded RNAs.

### 8.1.3 Quantifying RNA-type composition

After completing the hierarchical alignment:

1. Go to **Counting** module.
2. Load **all BAM files** for each sample:
   - `sample1.bam` (miRBase)
   - `sample1_tRNA.bam`
   - `sample1_rRNA.bam`
   - `sample1_genome.bam`
3. Generate count matrices separately for each reference, or use the **"Combined counting"** option if exposed in the UI (this creates a single matrix with features labeled by RNA type).

4. Go to **Post-Align QC** → **RNA Type Composition** (if available), or use the summary table in the Counting module.

Typical output:

| Sample | miRNA | tRNA | rRNA | Genome | Unaligned | Total |
|---|---|---|---|---|---|---|
| Control_rep1 | 8,234,567 (68%) | 1,456,789 (12%) | 345,678 (3%) | 567,890 (5%) | 1,456,123 (12%) | 12,000,000 |
| Drought_rep1 | 7,890,123 (65%) | 1,678,901 (14%) | 234,567 (2%) | 678,901 (6%) | 1,567,890 (13%) | 12,000,000 |

**Interpretation:**

- **High miRNA fraction (60–80%)**: Expected for well-prepared small RNA libraries.
- **High tRNA fraction (>20%)**: May indicate tRNA fragment accumulation under stress (biologically interesting) or tRNA contamination (technical artifact).
- **High rRNA fraction (>10%)**: Poor rRNA depletion; consider filtering these reads before DE analysis.
- **High genome fraction (>15%)**: May indicate novel small RNA loci or genomic repeat-derived siRNAs.

### 8.1.4 Exporting hierarchical alignment results

For publication or downstream analysis:

1. Export the **RNA-type composition table** as CSV.
2. Generate a **stacked bar chart** showing the proportion of each RNA type per sample (available in **Reports** module or via custom plotting).
3. Include the alignment parameters and reference versions in the **Methods** section:

   > "Hierarchical alignment was performed using Bowtie v1.3.1 with the following sequential references: (i) miRBase v22.1 (mature miRNAs, *Medicago truncatula*); (ii) RNAcentral v21 (tRNA and rRNA sequences, taxid 3880); (iii) *M. truncatula* genome v5.0 (Mt5.0). Alignment parameters were: -v 1 -k 10 --best --strata. Reads were classified as the first reference to which they aligned; subsequent references received only previously unaligned reads."

***

## 8.2 isomiR detection and analysis

### 8.2.1 What are isomiRs?

isomiRs (isoforms of miRNAs) are sequence variants of canonical miRNAs that arise from:

- **5′ or 3′ trimming/addition** – length variants due to imprecise Drosha/Dicer cleavage or post-transcriptional modifications.
- **Internal sequence changes** – SNPs, A-to-I editing.
- **Non-templated additions (NTA)** – uridylation or adenylation at the 3′ end.

isomiRs can have distinct biological functions (different target sets) or represent processing noise. Characterizing isomiR diversity is important for:

- Validating miRNA annotation (canonical vs. variant abundance).
- Identifying functional isoform switches (e.g., 5′ shift changes seed sequence).
- Understanding regulatory complexity.

### 8.2.2 Running isomiR detection

The **isomiR Analysis** module compares aligned reads to canonical miRNA sequences and classifies variants.

#### Step 1: Input data

1. Go to **isomiR Analysis** module.
2. Select **aligned BAM files** from the miRBase alignment (not the hierarchical secondary alignments).
3. Select the **miRBase reference FASTA** used for alignment (the module needs the canonical sequences for comparison).

#### Step 2: Configure detection parameters

- **Minimum read count**: Only analyze isomiRs with ≥ *N* reads (default: 10). This filters out sequencing errors.
- **Maximum edit distance**: Allow up to *N* differences from the canonical sequence (default: 3). This includes length differences + internal mismatches.
- **Classification categories**:
  - 5′ variant (trimmed or extended)
  - 3′ variant (trimmed or extended)
  - Internal SNP
  - Non-templated addition (NTA)
  - Canonical (exact match)

#### Step 3: Run analysis

Click **"Detect isomiRs"**. The module will:

1. Extract all reads aligned to each miRNA feature.
2. Compare each read sequence to the canonical mature sequence from miRBase.
3. Classify the variant type.
4. Count the abundance of each isoform.

#### Step 4: Inspect results

The results table includes:

| miRNA | Variant_type | Variant_sequence | Read_count | Fraction_of_miRNA | Sample |
|---|---|---|---|---|---|
| mtr-miR160a | canonical | UGCCUGGCUCCCUGUAUGCCA | 12,345 | 0.68 | Control_rep1 |
| mtr-miR160a | 3p_trimming | UGCCUGGCUCCCUGUAUGCC | 3,456 | 0.19 | Control_rep1 |
| mtr-miR160a | 3p_addition | UGCCUGGCUCCCUGUAUGCCAA | 1,234 | 0.07 | Control_rep1 |
| mtr-miR160a | nta_uridylation | UGCCUGGCUCCCUGUAUGCCAU | 890 | 0.05 | Control_rep1 |

**Key metrics:**

- **Canonical fraction**: Proportion of reads matching the reference sequence exactly. High fractions (>70%) indicate well-annotated miRNAs.
- **Dominant isoform**: The most abundant variant. If a non-canonical isoform is more abundant than the canonical, this may warrant further investigation.

#### Step 5: Visualizations

The module provides:

1. **Stacked bar chart**: Proportion of each variant type per miRNA.
2. **Length distribution**: Histogram of read lengths for each miRNA (shows trimming/extension patterns).
3. **Logo plot** (optional): Sequence logo showing nucleotide variability at each position.

#### Step 6: Export and interpretation

Export the isomiR table as CSV for:

- **Differential isoform analysis**: Compare isoform ratios between conditions (e.g., does drought shift the 5′ variant distribution?).
- **Target prediction refinement**: Re-run target prediction using the dominant isoform sequence if it differs from the canonical.
- **Annotation curation**: Submit novel dominant isoforms to miRBase.

**Example interpretation:**

> "We detected 1,234 isomiR variants across 89 miRNA families in *Medicago truncatula* drought-stressed leaves. The canonical sequence accounted for 68 ± 12% of reads per miRNA (mean ± SD). For mtr-miR160a, 3′ trimming variants represented 19% of reads, with the −1 nt isoform (UGCCUGGCUCCCUGUAUGCC) being the dominant variant. This 3′ shortening does not affect the seed region and is likely a processing artifact rather than a functional isoform."

***

## 8.3 Novel miRNA discovery

### 8.3.1 Rationale

Novel miRNA discovery aims to identify unannotated miRNAs—particularly relevant for:

- Non-model organisms with incomplete miRBase annotation.
- Tissue- or condition-specific miRNAs expressed only under certain circumstances.
- Species-specific miRNAs not conserved across kingdoms.

**Important caveat**: sRNAtlas provides a *candidate screening* mode, not a full *de novo* locus annotation pipeline. For rigorous plant miRNA discovery, **ShortStack** (Axtell, 2013; Johnson *et al.*, 2016) remains the gold standard, as it incorporates:

- Strand-specific clustering.
- Phased siRNA detection (21-nt or 24-nt phasing).
- Stringent hairpin structure validation via RNAfold.
- DICER-LIKE cleavage precision metrics.

sRNAtlas's **Novel miRNA** module is designed to:

- Flag high-confidence candidates quickly.
- Provide an entry point for users unfamiliar with command-line tools.
- Integrate with the rest of the sRNAtlas workflow (DE, targets, enrichment for novel candidates).

Candidates identified here should be validated experimentally or cross-checked with ShortStack.

### 8.3.2 Workflow

#### Step 1: Prepare unaligned reads

You need reads that did **not** align to known miRNAs. Sources:

1. **From hierarchical alignment** (Section 8.1): Use the final unaligned FASTQ after miRBase → tRNA → rRNA alignment.
2. **From standard alignment**: Use the `<sample>_unaligned.fastq.gz` output from the Alignment module with "Save unaligned reads" enabled.

These reads may contain:

- Novel miRNAs.
- Unannotated small RNAs (siRNAs, novel tRFs, etc.).
- Low-quality reads or adapter artifacts.

#### Step 2: Configure discovery parameters

Go to **Novel miRNA Discovery** module.

**Input section:**

- Upload the unaligned FASTQ files, or select them from the session if they were generated in the current project.
- Upload or select a **genome reference FASTA** (required for hairpin prediction if implemented, or for genomic context mapping).

**Parameter settings:**

1. **Size selection**:
   - Minimum length: 18 nt (default)
   - Maximum length: 24 nt (default)
   - Rationale: Canonical miRNAs are typically 20–24 nt. Adjust if targeting other small RNA classes.

2. **Read abundance filter**:
   - Minimum read count: 100 (default)
   - Rationale: True miRNAs are usually abundant; low-count sequences are more likely to be degradation products or errors.

3. **5′ nucleotide bias** (optional but recommended for plants):
   - Check **"Prefer 5′ Uracil (U/T)"**.
   - Rationale: Most plant miRNAs start with U (AGO1 sorting signal). This is less strict for animals.

4. **Hairpin structure validation** (if implemented):
   - Check **"Validate hairpin potential"**.
   - Minimum free energy (MFE) threshold: −20 kcal/mol (default).
   - Rationale: miRNA precursors fold into stable stem-loop structures. More negative MFE = more stable.
   - **Current limitation** (as of Feb 2026): The sRNAtlas UI captures this setting, but RNAfold integration is not yet complete. A disclaimer should appear in the interface. Candidates are currently filtered by sequence criteria only (length, count, 5′ bias).

#### Step 3: Run discovery

Click **"Start Discovery"**. The module will:

1. Extract unique sequences from unaligned reads.
2. Filter by length and count thresholds.
3. Optionally apply 5′ nucleotide filter.
4. (If hairpin validation is enabled and implemented): Extract genomic flanking regions, predict secondary structure via RNAfold, and filter by MFE.
5. Score candidates based on:
   - Read abundance (higher = better).
   - 5′ U bias (if enabled).
   - GC content (miRNAs typically have 40–60% GC).
   - (If hairpin validation works): Structure stability.

#### Step 4: Inspect candidates

The results table shows:

| Candidate_ID | Sequence | Length | Read_count | 5′_nt | GC% | Score | Hairpin_MFE |
|---|---|---|---|---|---|---|---|
| Novel_001 | UGCCUGGCUCCCUGUAUGCCA | 21 | 4,567 | U | 52.4 | 0.87 | −28.3 |
| Novel_002 | UAGAAGCUGCCAGCAUGAUCU | 21 | 3,456 | U | 47.6 | 0.82 | −22.1 |
| Novel_003 | AAGGCUUCGUCGUCGUCGUC | 20 | 2,345 | A | 55.0 | 0.65 | −18.5 |

**Interpretation:**

- **Novel_001 and Novel_002**: High score, start with U, abundant, stable hairpin (if MFE available). **High-confidence candidates**.
- **Novel_003**: Starts with A (lower confidence for plant miRNAs), weaker hairpin. **Medium-confidence candidate** or may be a siRNA.

#### Step 5: Validation and downstream analysis

**Experimental validation:**

- Northern blot to confirm expression.
- Small RNA sequencing of precursor-overexpressing lines.
- qRT-PCR to detect the precursor transcript.

**Computational validation:**

1. **BLAST the candidate sequence** against miRBase to ensure it's not a known miRNA from a related species (may indicate incomplete annotation for your species).
2. **Map the candidate to the genome** (if not already done) to identify the locus and extract flanking sequence.
3. **Fold the precursor** using a standalone RNAfold run (ViennaRNA package) to manually verify the hairpin structure:

   ```bash
   echo ">candidate_001_precursor" > candidate_001.fa
   echo "AGCUGCCUGGCUCCCUGUAUGCCAAGCUAGCUGCAUGCUAGCU..." >> candidate_001.fa
   RNAfold candidate_001.fa
   ```

   Look for:
   - A single, unbranched hairpin.
   - The mature sequence located in one arm of the hairpin.
   - No large internal loops or bulges near the mature region.

4. **Run ShortStack** on your original BAM files (aligned to genome) to get a rigorous *de novo* annotation:

   ```bash
   ShortStack --readfile alignment.bam \
              --genomefile genome.fa \
              --outdir shortstack_output \
              --mincov 100 \
              --dicermin 20 --dicermax 24
   ```

   Compare sRNAtlas candidates with ShortStack-annotated MIRNA loci. ShortStack may identify the same locus with additional confidence metrics (Dicer call, strand ratio, phasing).

#### Step 6: Integrate into the main analysis

If you validate a novel miRNA:

1. Add its sequence to your **custom reference FASTA**.
2. Re-run alignment and counting to include the novel miRNA in the count matrix.
3. Perform DE analysis to see if the novel miRNA is differentially expressed.
4. Run target prediction (psRNATarget or local seed matching) to identify potential targets.
5. Include in the enrichment analysis if targets are identified.

***

## 8.4 Combined workflow example: Medicago drought stress

Here's how you might combine all three advanced modules in a real project:

**Experimental setup:**

- Organism: *Medicago truncatula*
- Samples: 4 control, 4 drought-stressed leaf tissue
- Goal: Identify miRNAs, other small RNAs, and novel candidates responding to drought.

**Pipeline:**

1. **Tutorials 1–5**: Trim, QC, align to miRBase, count, DE analysis.
   - Result: 45 differentially expressed miRNAs.

2. **Tutorial 8.1 (Hierarchical alignment)**:
   - Take unaligned reads from Step 1.
   - Align to RNAcentral tRNA → identify tRNA fragments.
   - Align remaining reads to RNAcentral rRNA → remove contaminants.
   - Result: 12% of reads are tRFs (up-regulated under drought), 2% rRNA.

3. **Tutorial 8.2 (isomiR analysis)**:
   - Analyze miRBase-aligned reads for isoform diversity.
   - Result: For mtr-miR398, the 5′ +1 nt isoform is the dominant form under drought (shifts the seed sequence → different targets).

4. **Tutorial 8.3 (Novel miRNA discovery)**:
   - Use the final unaligned reads from Step 2 (post-rRNA removal).
   - Apply 5′ U filter, minimum 100 reads, length 20–24 nt.
   - Result: 8 high-confidence candidates. BLAST confirms 3 are known miRNAs from *Medicago sativa* (not in *M. truncatula* miRBase). The other 5 are truly novel.
   - ShortStack validation: 4 of the 5 novel candidates confirmed as MIRNA loci.

5. **Downstream (Tutorial 6)**:
   - Add the 4 validated novel miRNAs to the custom reference.
   - Re-count and re-run DE → 2 of the novel miRNAs are drought-responsive.
   - Predict targets for the novel miRNAs using psRNATarget.
   - Enrichment: Targets of novel-miRNA-001 are enriched for "response to water deprivation" (GO:0009414).

**Result**: A comprehensive dissection of the drought-responsive small RNAome, including known miRNAs, isomiRs, tRNA fragments, and novel drought-specific miRNAs.

***

## 8.5 Best practices and considerations

### 8.5.1 When to use hierarchical alignment

- **Always**, if you care about RNA-type composition (QC, contamination assessment).
- **Essential** for studies of tRNA fragments, rRNA-derived fragments, or multi-class small RNA biology.
- **Less critical** if you only care about miRNA differential expression and your libraries are high-quality (>80% miRNA).

### 8.5.2 When to analyze isomiRs

- When studying **miRNA processing** or **post-transcriptional modifications**.
- When the **dominant isoform differs from the canonical** (important for target prediction).
- When you need to **validate miRNA annotations** (is the miRBase sequence the actual predominant form in your tissue?).

### 8.5.3 When to run novel miRNA discovery

- **Non-model organisms** with incomplete miRBase annotation.
- **Tissue- or stress-specific studies** where novel miRNAs may be expressed.
- **After checking miRBase** – if your organism already has comprehensive annotation (e.g., human, mouse, *Arabidopsis*), the yield of true novel miRNAs will be low.

### 8.5.4 Computational resource considerations

- **Hierarchical alignment**: Multiplies alignment time by the number of references (3–4× longer than single alignment). Plan accordingly for large datasets.
- **isomiR analysis**: Memory-intensive if BAM files are large (many samples, deep sequencing). Subsample or process samples individually if memory is limited.
- **Novel discovery**: Requires a genome reference (can be large). If RNAfold is integrated, structure prediction is CPU-intensive (~1 minute per candidate per sample on a single core). Use multi-threading if available.

***

## 8.6 Summary

This tutorial demonstrated:

- **Hierarchical alignment** to classify reads by RNA type (miRNA → tRNA → rRNA → genome), providing comprehensive RNA-type composition analysis.
- **isomiR detection** to characterize sequence variants of known miRNAs and identify dominant isoforms that may differ from miRBase annotations.
- **Novel miRNA discovery** to screen unaligned reads for candidate unannotated miRNAs, with guidelines for validation using ShortStack and experimental methods.

These advanced analyses complement the core pipeline and enable:

- More detailed biological interpretation.
- Discovery of non-miRNA small RNA classes.
- Identification of species- or condition-specific regulatory RNAs.

**Next tutorial:**

- **Tutorial 9** will cover **Visualization and Reporting**: generating publication-quality figures, HTML reports, and exporting data for external tools.
