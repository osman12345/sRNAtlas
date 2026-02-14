# Tutorial 3: Quality Control Deep Dive

How to read every QC metric in sRNAtlas, diagnose problematic libraries, and make informed decisions before committing to downstream analysis.

---

## Overview

| | |
|---|---|
| **Learning objectives** | Interpret size distributions for different organisms and RNA types; diagnose adapter contamination, degradation, and library failures; compare pre- vs. post-trimming QC; use post-alignment QC to validate biological expectations |
| **Prerequisites** | [Tutorial 2: Your First Analysis](02-first-analysis.md) completed; familiarity with the sRNAtlas interface |
| **Estimated time** | 30–40 minutes (reading + hands-on) |
| **Difficulty** | Intermediate |

---

## Why QC Matters More for Small RNA-seq

Quality control in small RNA-seq is fundamentally different from standard RNA-seq. In a typical mRNA-seq experiment, reads are 75–150 nt and map to exons; QC focuses on duplication rates, gene body coverage, and 3' bias. In small RNA-seq, the insert itself is only 18–32 nt — shorter than the sequencing read length — which creates unique challenges:

- **Adapter read-through is expected**, not a failure. Every read contains the small RNA insert followed by adapter sequence.
- **Read length distribution is the single most informative QC metric.** It reveals which RNA classes dominate your library before any annotation.
- **Small variations in trimming parameters can discard real biology or retain artifacts.** A 1-nt shift in the minimum length filter changes which miRNAs survive.
- **Contamination signatures are size-dependent.** rRNA fragments, adapter dimers, and degradation products each occupy characteristic size ranges.

This tutorial teaches you to read every QC output in sRNAtlas like a diagnostic — understanding not just *what* the numbers are, but *what they mean* and *what to do about them*.

---

## Part 1: Pre-Trimming QC — Raw Read Assessment

### 1.1 Running the QC Module

1. Click **Quality Control** in the sidebar
2. In the input section, select **Use files from Project**
3. Click **Run QC Analysis**

sRNAtlas samples up to 100,000 reads per file (configurable) for fast analysis. The results appear in three sections: **Read Statistics Summary**, **Size Distribution**, and **Contamination Check**.

### 1.2 Read Statistics Summary

The summary table shows per-sample metrics:

| Metric | What it measures | Healthy range | Red flag |
|--------|-----------------|---------------|----------|
| **Total Reads** | Sequencing depth | >1M for real experiments; ~50k for tutorial data | <500k may lack statistical power for DE |
| **Mean Length** | Average read length before trimming | Equal to your sequencing cycle length (e.g., 50 nt for SE50) | If much shorter than expected, reads may be pre-trimmed or truncated |
| **Mean Quality** | Average Phred score across all bases | >28 | <20 indicates serious quality problems |
| **GC %** | Guanine + cytosine content | 45–55% for most organisms | >60% or <35% suggests contamination or extreme bias |

The four metric cards at the top (Total Reads, Avg Quality, Avg GC%, Avg Length) give you an instant cross-sample overview. **The most important check at this stage is consistency** — all samples from the same experiment should have similar values. If one sample has dramatically different read counts or quality, flag it immediately.

> **What "Mean Length" tells you before trimming**: If your sequencing was SE50 (single-end, 50 cycles), all raw reads should be exactly 50 nt. A mean length of 50.0 confirms the data is untrimmed. If you see a mean length of ~25 nt, your data was already adapter-trimmed by the sequencing facility — skip to Step 5 (Trimming) or proceed directly to alignment.

### 1.3 Size Distribution — Before Trimming

Before adapter removal, the size distribution plot is **not yet informative about RNA biology**. You will typically see:

- **A single spike at the sequencing read length** (e.g., 50 nt) — this is normal for untrimmed data
- Or a **broad distribution centered on the read length** if quality trimming was applied upstream

This plot becomes biologically meaningful only after trimming (see Part 2). At this stage, its main purpose is confirming whether your data has been pre-processed.

### 1.4 Contamination Check

Click **Run Contamination Analysis** to scan for adapter sequences, poly-A/T content, and overrepresented k-mers.

#### Adapter Content

sRNAtlas searches for 12-bp signature sequences from multiple platforms:

| Platform | Adapter name | 12-bp signature |
|----------|-------------|-----------------|
| Illumina | Universal | `AGATCGGAAGAG` |
| Illumina | Small RNA (RA3) | `TGGAATTCTCGG` |
| Illumina | Nextera | `CTGTCTCTTATA` |
| Ion Torrent | Small RNA | `ATCACCGACTGC` |
| BGI/MGI | Standard | `AGTCGGAGGCCA` |
| Nanopore | cDNA | `ACTTGCCTGTCG` |
| QIAseq | miRNA | `AACTGTAGGCAC` |

**Interpretation guide:**

| Adapter content level | Meaning | Action |
|----------------------|---------|--------|
| **>80% Illumina Small RNA** | Expected and healthy — confirms small RNA inserts are shorter than read length | Proceed to trimming with the matched adapter preset |
| **>50% Illumina Small RNA** | Normal — some reads may be longer ncRNAs | Proceed to trimming |
| **<20% for any adapter** | Unusual — reads may already be trimmed, or the wrong adapter was used | Check if data is pre-trimmed; verify library prep protocol |
| **Multiple adapters detected** | Possible index hopping or mixed libraries | Investigate demultiplexing; check for sample swap |

> ⚠️ **Critical insight**: High adapter content in raw small RNA-seq data is a *good sign*, not a problem. It confirms that your inserts are genuine small RNAs (18–25 nt) that are shorter than the read length. This is the opposite of standard RNA-seq, where adapter content indicates a problem (short inserts or adapter dimers).

#### Poly-A/T Content

| Poly content | Meaning |
|-------------|---------|
| **Poly-A <1%** | Normal |
| **Poly-A 1–5%** | May indicate poly-A tailing in library prep (e.g., Takara SMARTer, Diagenode CATS kits) |
| **Poly-A >5%** | Possible mRNA contamination or specific library prep artifact |
| **Poly-T >5%** | Unusual — investigate library prep |

#### Overrepresented 6-mers

sRNAtlas reports the top overrepresented 6-mers (>0.5% of all sampled k-mers). Common benign hits include:

- `AGATCG` — Illumina adapter signature
- `TGGAAT` — Illumina small RNA adapter (RA3) start
- `AAAAAA` — Poly-A (if present)
- Organism-specific sequences (e.g., highly abundant miRNAs like miR-166 in plants)

**Suspicious k-mers** to investigate: sequences that don't match known adapters *and* don't correspond to abundant miRNAs may indicate contamination from another organism or a spike-in control.

---

## Part 2: Post-Trimming QC — The Critical Assessment

After running the **Trimming** module (Tutorial 2, Step 5), return to the **Quality Control** module and re-run the analysis on your **trimmed files**. This is where the real diagnostic power lies.

### 2.1 Size Distribution — After Trimming

The post-trimming size distribution is **the single most important QC plot in small RNA-seq**. It reveals the RNA composition of your library at a glance, before any alignment or annotation.

sRNAtlas overlays colored shaded regions for expected RNA type size ranges on the histogram:

| Color region | RNA type | Size range | What to expect |
|-------------|----------|------------|----------------|
| Blue | miRNA | 18–25 nt | Dominant peak at 21–22 nt in most experiments |
| Orange | siRNA | 20–24 nt | Overlaps with miRNA range |
| Green | piRNA | 24–32 nt | Prominent in animal germline tissues |
| Red | tRF | 28–40 nt | Stress-responsive; variable |
| Purple | snoRNA | 60–90 nt | Only if library prep captures longer RNAs |
| Brown | snRNA | 100–150 nt | Rare in standard sRNA-seq |

sRNAtlas also performs automatic **peak detection** — it identifies the top 5 most abundant read lengths and reports their counts and percentages. A message interprets the dominant peak:

- Peak at 20–24 nt → "✅ Peak consistent with **miRNA**"
- Peak at 26–32 nt → "ℹ️ Peak consistent with **piRNA**"
- Peak at >40 nt → "⚠️ Peak length suggests reads may not be trimmed or are not small RNA"

### 2.2 Reading Size Distributions by Organism

Different organisms produce characteristic size distribution profiles. Learning to recognize these patterns is one of the most valuable QC skills in small RNA biology.

#### Animal Tissues (Human, Mouse, Drosophila)

```
Expected pattern:
    │
    │      ██
    │     ████
    │    ██████
    │   ████████
    │  ██████████
    └──────────────────
    18  21  24  28  32  nt

Dominant peak: 22 nt (miRNAs)
Secondary peak: ~30 nt (piRNAs, in germline tissues only)
```

- **Somatic tissues** (liver, brain, blood): Single sharp peak at 22 nt. miRNAs dominate (>70% of mapped reads). Minimal signal at other sizes.
- **Germline tissues** (testis, ovary): Bimodal distribution — a 22-nt miRNA peak *and* a broad 26–30 nt piRNA peak. In mouse testis, piRNAs can exceed miRNAs in abundance.
- **Extracellular vesicles / biofluids** (serum, plasma): Broad distribution with peaks at 22 nt (miRNAs) and 30–33 nt (tRNA halves). Higher adapter dimer contamination is common due to low input RNA.

#### Plant Tissues (Arabidopsis, Medicago, Rice)

```
Expected pattern:
    │
    │  ██            ██
    │ ████          ████
    │██████        ██████
    │████████    ████████
    └──────────────────────
    18  21  24  27  30  nt

Dominant peaks: 21 nt AND 24 nt
```

- **Two characteristic peaks**: 21 nt (miRNAs + ta-siRNAs) and 24 nt (heterochromatic siRNAs / hc-siRNAs). This bimodal pattern is a hallmark of plant sRNA-seq and reflects two distinct silencing pathways:
  - 21 nt: DCL1-processed miRNAs (loaded into AGO1) and DCL4-processed ta-siRNAs
  - 24 nt: DCL3-processed hc-siRNAs (loaded into AGO4, guide RNA-directed DNA methylation)
- **Relative ratio matters**: In leaves, the 24-nt peak often exceeds the 21-nt peak. In flowers/inflorescences, the 21-nt peak is typically stronger. In reproductive tissues, 21-nt and 24-nt phasiRNAs can be very abundant.
- **22-nt peak** (if present): Secondary siRNAs triggered by 22-nt miRNAs (e.g., miR173, miR472 in *Arabidopsis*). These are biologically meaningful.

#### Fungi, Bacteria, Archaea

- **Fungi** (e.g., *Neurospora*): Peak at 21–25 nt (milRNAs and siRNAs if RNAi is active). Many fungi lack RNAi entirely.
- **Bacteria**: No canonical small RNA-seq profile. Expect tRNA fragments and rRNA fragments to dominate.

### 2.3 Diagnostic Patterns — What Went Wrong?

Not every library produces a clean size distribution. Here are the most common failure patterns and their diagnoses:

#### Pattern A: Spike at Very Short Reads (15–17 nt)

```
    │██
    │██
    │██ ██
    │██████
    │████████
    └──────────────
    15  18  21  24  nt
```

**Diagnosis**: **Adapter dimers.** When adapters ligate to each other without an insert, the resulting reads after trimming are very short (the few nucleotides between the adapter junctions). Adapter dimer contamination is one of the most common artifacts in sRNA-seq. A recent study showed it affects ~28% of publicly available sRNA-seq datasets and correlates with batch effects and sequencing failure.

**What to do**:
- If your min length filter is set to 18 nt (the sRNAtlas default), these should already be filtered out
- If the spike represents >30% of total reads, the library has low complexity — consider whether enough reads remain after filtering
- For future experiments: optimize the insert:adapter ratio during library prep, or use gel size-selection to remove dimers

#### Pattern B: Flat Distribution, No Clear Peak

```
    │
    │████████████████
    │████████████████
    │████████████████
    └──────────────────
    18  21  24  28  32  nt
```

**Diagnosis**: **RNA degradation or non-specific library.** Degraded RNA produces fragments across all size ranges. Alternatively, the library prep may have captured random degradation products rather than genuine small RNAs.

**What to do**:
- Check the RNA Integrity Number (RIN) from your Bioanalyzer/TapeStation results (if available). Low RIN (<7) confirms degradation.
- Proceed with analysis but expect high background noise and lower sensitivity for detecting differentially expressed miRNAs.
- For future experiments: use higher quality RNA input and verify small RNA enrichment on a Bioanalyzer Small RNA chip.

#### Pattern C: Dominant Peak at 28–33 nt, Weak miRNA Peak

```
    │
    │              ████
    │   █         ██████
    │  ███       ████████
    │ █████     ██████████
    └──────────────────────
    18  21  24  28  32  nt
```

**Diagnosis**: **tRNA fragment / tRNA half dominance.** tRNA-derived fragments (tRFs) and tRNA halves (tiRNAs) are 28–36 nt and can overwhelm miRNA signal, especially in samples under stress, from biofluids, or when the library was prepared from total RNA without size selection.

**What to do**:
- This is not necessarily a problem — tRFs are biologically real and increasingly studied.
- If your goal is miRNA analysis, you can computationally filter tRF-derived reads after alignment (align to a tRNA reference first, then remove those reads).
- If you aligned to miRBase only, these reads will appear as "unaligned" — they are not lost.

#### Pattern D: Single Spike at the Sequencing Read Length

```
    │                         ██
    │                         ██
    │                         ██
    │                         ██
    └──────────────────────────
    18  21  24  28  40  50  nt
```

**Diagnosis**: **Trimming failed or wrong adapter.** Reads retained their full length, meaning the adapter was not found and removed.

**What to do**:
- Return to the Trimming module and try a different adapter preset
- If you don't know which adapter was used, try the top 3 most common:
  1. `Illumina TruSeq Small RNA (RA3)` — `TGGAATTCTCGGGTGCCAAGG`
  2. `NEBNext Small RNA` — `AGATCGGAAGAGCACACGTCT`
  3. `QIAseq miRNA` — `AACTGTAGGCACCATCAAT`
- Check your library preparation protocol documentation or contact the sequencing provider

#### Pattern E: Peak at 24 nt Only (in Plants)

```
    │
    │            ████
    │           ██████
    │          ████████
    │        ████████████
    └──────────────────────
    18  21  24  28  32  nt
```

**Diagnosis**: **hc-siRNA-dominated library** or tissue with high repeat/transposon content. In plants, certain tissues (e.g., endosperm, pollen) or genotypes produce very high levels of 24-nt hc-siRNAs.

**What to do**:
- This is biologically valid — 24-nt hc-siRNAs are involved in RNA-directed DNA methylation (RdDM).
- If your goal is miRNA analysis, be aware that miRNAs will be a smaller fraction of total reads. You may need deeper sequencing.
- Consider whether this matches your experimental expectations (tissue type, mutant genotype, stress condition).

### 2.4 The "Reads by Expected RNA Type" Table

After the size distribution plot, sRNAtlas displays a categorization table:

| RNA Type | Count | Percentage |
|----------|-------|------------|
| miRNA (18–25 nt) | 34,200 | 68.4% |
| siRNA (20–24 nt) | 28,100 | 56.2% |
| piRNA (24–32 nt) | 8,900 | 17.8% |
| tRF/tsRNA (14–40 nt) | 41,500 | 83.0% |
| Short (<14 nt) | 200 | 0.4% |
| Long (>50 nt) | 1,100 | 2.2% |

> ⚠️ **Important**: These size ranges overlap. A 22-nt read is counted in *both* the miRNA and siRNA categories. A 24-nt read is counted in miRNA, siRNA, *and* piRNA ranges. The percentages will sum to more than 100%. **This is by design** — size alone cannot definitively assign an RNA type. Actual annotation (alignment to databases) is required for definitive classification. This table shows *potential* composition based on size.

---

## Part 3: Post-Alignment QC — Validating Your Mapping

After alignment (Tutorial 2, Steps 7–8), the **Post-Align QC** module provides mapping-level diagnostics.

1. Click **Post-Align QC** in the sidebar
2. Click **Run Post-Alignment QC**

### 3.1 Mapping Statistics

| Metric | Definition | Healthy range | Concern |
|--------|-----------|---------------|---------|
| **Mapping rate** | % reads aligned to reference | >70% (miRBase) | <30% (wrong reference or untrimmed reads) |
| **Unique mapping** | % mapped reads with a single best hit | >50% | <30% (many multi-family miRNAs, or broad reference) |
| **Multi-mapped** | Reads aligning to ≥2 reference sequences equally well | <50% | >70% (reference contains redundant sequences) |
| **Forward strand** | Reads mapping to the sense strand | ~100% for mature miRNA ref | ~50% suggests non-specific mapping |
| **Reverse strand** | Reads mapping to antisense | ~0% for mature miRNA ref | High reverse mapping to mature sequences is suspicious |

**Why do multi-mappers occur in miRNA-seq?**

Many miRNA genes exist as families with nearly identical mature sequences. For example, in *Arabidopsis*:
- ath-miR166**a** through ath-miR166**g** all produce the same or nearly identical 21-nt mature sequence
- A read matching this sequence maps equally well to all 7 loci

sRNAtlas handles this with `-k 10` (report up to 10 alignments) and `--best --strata` (only the best-scoring alignments). Multi-mapper rates of 30–50% are normal for miRNA analysis.

### 3.2 Strand Bias

When aligning to **mature miRNA sequences** from miRBase, you expect essentially 100% forward-strand mapping. The mature sequences in miRBase are already in the 5'→3' orientation of the functional miRNA, so trimmed reads should align in the same direction.

| Strand ratio (fwd:rev) | Interpretation |
|------------------------|----------------|
| >95% forward | Normal for mature miRNA reference |
| ~50:50 | Suspicious — may indicate mapping to a genomic reference rather than mature sequences, or non-specific alignments |
| >95% reverse | Data may be reverse-complemented; check FASTQ orientation |

If you aligned to **hairpin sequences** instead of mature, expect both forward and reverse reads (the hairpin contains both the miRNA and miRNA* strands).

### 3.3 5' Nucleotide Composition

The 5' nucleotide of a small RNA determines which Argonaute protein loads it. This is a powerful biological validation:

#### Animals (Human, Mouse)

| 5' nucleotide | Expected enrichment | Loaded into | Biology |
|---------------|-------------------|-------------|---------|
| **U (T in DNA)** | ~50–60% of miRNAs | AGO1, AGO2, AGO3 | Most canonical miRNAs start with U |
| **A** | ~20–25% | AGO2, AGO4 | Some miRNA families; piRNAs show 1U and 10A (ping-pong signature) |
| **C** | ~15% | AGO1 | Less common |
| **G** | ~10% | AGO4 | Least common for miRNAs |

#### Plants (Arabidopsis, Medicago, Rice)

| 5' nucleotide | Expected enrichment | Loaded into | Biology |
|---------------|-------------------|-------------|---------|
| **U** | Dominant for 21-nt reads | **AGO1** | miRNAs — post-transcriptional gene silencing |
| **A** | Dominant for 24-nt reads | **AGO4** | hc-siRNAs — RNA-directed DNA methylation |
| **C** | Enriched in 22-nt reads | **AGO5** | Some specialized small RNAs |
| **G** | Minor | Various | Least common |

**Diagnostic value**: If you see strong 5'-U enrichment in 21-nt reads and 5'-A enrichment in 24-nt reads in a plant sample, this is strong evidence that your library has captured genuine, biologically active small RNAs.

> **Tip**: In sRNAtlas, the 5' nucleotide is computed from the `query_sequence` of each aligned read. For reverse-strand alignments, the 5' end is taken from the last base of the stored sequence (sRNAtlas handles the strand orientation automatically in the Post-Align QC module).

### 3.4 Read Length Distribution (Post-Alignment)

The post-alignment length distribution shows the sizes of reads that *successfully mapped* to your reference. Comparing this to the pre-alignment (post-trimming) distribution tells you which size classes are being captured:

| Scenario | Pre-alignment peak | Post-alignment peak | Interpretation |
|----------|-------------------|--------------------|----|
| Normal miRNA analysis | 21 nt | 21 nt | Most reads in the miRNA range map successfully |
| Plant with miRBase-only ref | 21 + 24 nt | 21 nt only | 24-nt hc-siRNAs don't map to miRBase (expected — miRBase contains miRNAs, not siRNAs) |
| Low mapping rate | 21 nt | 21 nt (few reads) | Many reads in the miRNA size range don't match known miRNAs — candidates for novel miRNA discovery |

### 3.5 RNA Type Composition

If your reference includes multiple RNA types (e.g., from RNAcentral), sRNAtlas displays a pie chart of reads mapping to each category. A typical healthy profile:

| RNA type | Expected fraction | Concern threshold |
|----------|------------------|-------------------|
| miRNA | 50–80% | <30% in a miRNA-focused experiment |
| tRNA/tRF | 5–20% | >50% indicates tRNA contamination |
| rRNA/rsRF | 1–10% | >20% indicates rRNA depletion failure |
| snoRNA | 1–5% | — |
| Other/unknown | 5–15% | >50% suggests annotation gaps |

---

## Part 4: Cross-Sample Comparison — Spotting Outliers

One of the most important QC tasks is comparing metrics **across all samples** in your experiment. sRNAtlas displays per-sample statistics side by side in the summary table.

### What to Compare

| Metric | Acceptable variation | Flag if |
|--------|---------------------|---------|
| Total reads | Within 3-fold | One sample has 10× fewer reads than others |
| Trimming pass rate | Within 10% | One sample is 30% lower than the rest |
| Alignment rate | Within 15% | One sample maps at 20% while others map at 70% |
| Size distribution shape | Same peak positions | One sample has a completely different peak pattern |
| GC content | Within 5% | One sample is >10% different from others |

### Decision Framework: Keep, Flag, or Remove?

```
For each sample, ask:

1. Is the total read count adequate?
   ├─ YES → continue
   └─ NO (<500k) → FLAG: may lack statistical power

2. Is the size distribution consistent with other samples?
   ├─ YES → continue
   └─ NO → Is there a biological explanation (different tissue, genotype)?
           ├─ YES → KEEP but note in analysis
           └─ NO → FLAG as potential technical outlier

3. Is the alignment rate consistent with other samples?
   ├─ YES → continue
   └─ NO (>20% deviation) → Investigate:
           ├─ Contamination? → Run against rRNA/tRNA ref
           ├─ Sample swap? → Check metadata
           └─ Library failure? → REMOVE if confirmed

4. Does PCA (in DE module) cluster this sample with its group?
   ├─ YES → KEEP
   └─ NO → REMOVE from DE analysis, document reason
```

> **Golden rule**: Never silently remove outlier samples. Always document the QC evidence that led to exclusion, and report both the full and filtered analyses if the exclusion changes results substantially.

---

## Part 5: QC Checklist — Before Proceeding to DE Analysis

Before moving on to differential expression, verify every item:

### ✅ Pre-Analysis Checklist

- [ ] **All samples have sufficient reads** (>500k for real experiments; >20k for tutorial data)
- [ ] **Adapter was correctly identified** (>60% adapter content in raw reads, or high trimming pass rate)
- [ ] **Post-trimming size distribution shows expected peaks** (21 nt for miRNA; 21+24 nt for plants)
- [ ] **No adapter dimer spike** at 15–17 nt (or if present, it represents <30% of reads)
- [ ] **Alignment rate is consistent** across samples (within 15%)
- [ ] **5' nucleotide composition matches biology** (U-bias for miRNAs)
- [ ] **Strand bias is as expected** (~100% forward for mature miRNA reference)
- [ ] **No sample is a clear outlier** on all metrics simultaneously
- [ ] **rRNA contamination is low** (<20% of mapped reads)

If any item fails, address it before proceeding. The Troubleshooting table below maps each failure to its resolution.

---

## Part 6: Troubleshooting Reference

### Size Distribution Problems

| Symptom | Diagnosis | Resolution |
|---------|-----------|------------|
| All reads at 50 nt (full read length) | Adapter not trimmed | Try a different adapter preset in Trimming module |
| Spike at 15–17 nt | Adapter dimers | Increase min length to 18 (default); if >30% of reads, library may have low complexity |
| Flat distribution, no peak | RNA degradation | Check RNA quality (RIN); proceed with caution |
| Peak at 28–33 nt only | tRNA fragment dominance | Filter tRFs computationally; or accept if tRFs are your target |
| Only 24-nt peak in plants | hc-siRNA-dominated | Biologically valid; miRNA signal may be low |
| Peak at 40+ nt | Captured longer ncRNAs or degradation | Adjust max length filter; check library prep protocol |

### Alignment Problems

| Symptom | Diagnosis | Resolution |
|---------|-----------|------------|
| 0% alignment | Wrong reference or untrimmed reads | Verify organism; confirm trimming completed |
| <30% alignment | Partial mismatch (subspecies, strain) | Try `-v 2` (allow 2 mismatches); try RNAcentral in addition to miRBase |
| >90% alignment + >70% multi-mapped | Redundant reference sequences | Normal for miRBase (family members); consider collapsing counts by miRNA family |
| High reverse-strand mapping | Incorrect reference orientation or hairpin reference | Verify you used mature (not hairpin) sequences; check FASTQ read orientation |
| Alignment rate varies >20% between samples | Sample-specific contamination or quality | Investigate the outlier sample's QC metrics individually |

### Contamination Problems

| Symptom | Diagnosis | Resolution |
|---------|-----------|------------|
| rRNA >20% | Ribosome depletion failure | Filter rRNA reads computationally (align to rRNA first, discard mapped) |
| Multiple adapter types detected | Mixed libraries or index hopping | Re-demultiplex; check for sample contamination |
| High poly-A content | mRNA contamination or kit-specific (SMARTer, CATS) | If using poly-A–based kit, this is expected; otherwise investigate |
| Overrepresented k-mers from another organism | Cross-contamination | Align to potential contaminant genomes to quantify; consider excluding |

---

## Summary

Quality control in small RNA-seq is not a checkbox — it's a diagnostic process. The key takeaways:

1. **Size distribution is your most powerful tool.** Learn to recognize the characteristic patterns for your organism and tissue type.
2. **High adapter content is normal.** Don't panic when >80% of raw reads contain adapter — that's exactly what healthy small RNA data looks like.
3. **Consistency across samples matters more than absolute values.** An alignment rate of 40% is fine if all samples are at 40%. One sample at 40% while others are at 80% is a problem.
4. **Address QC failures before DE analysis.** Garbage in, garbage out — no statistical method can rescue a failed library.
5. **Document everything.** Record which QC checks each sample passed, and justify any sample exclusions with specific metrics.

---

## What's Next?

| Tutorial | Topic |
|----------|-------|
| [Tutorial 4: Databases & Alignment Strategy](04-databases-alignment.md) | Reference selection, Bowtie parameters, multi-database workflows |
| [Tutorial 5: Differential Expression](05-differential-expression.md) | From count matrix to publication-ready DE results |
| [Tutorial 6: Novel miRNA & isomiR Analysis](06-novel-mirna-isomir.md) | Discover unannotated small RNAs and sequence variants |

---

*Tutorial written for sRNAtlas v0.1.0 (BETA). Last updated: 2026-02-12.*
