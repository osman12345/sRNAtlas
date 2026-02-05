# sRNAtlas Quick Start Guide

Get your small RNA-seq analysis running in 10 minutes!

---

## Before You Begin

Ensure you have:
- FASTQ files from small RNA-seq
- Know your adapter sequence (or library prep kit)
- Know your organism

---

## Step 1: Create a Project (1 min)

1. Click **Project** in the sidebar
2. Go to **Project Settings** tab
3. Enter your project name
4. Select organism (optional)
5. Click **Create Project**

---

## Step 2: Upload Your Data (2 min)

1. Stay in **Project** module
2. Go to **Data Upload** tab
3. Drag & drop your `.fastq.gz` files
4. Click **Add Files to Project**

**Tip:** Files up to 2GB supported. Use gzip compression.

---

## Step 3: Quality Control (2 min)

1. Click **Quality Control** in sidebar
2. Go to **Upload Data** tab
3. Select **Use files from Project**
4. Click **Run QC Analysis**

**Check these tabs:**
- **Read Statistics**: Total reads, quality scores
- **Size Distribution**: Look for peak at 21-23 nt (miRNA)
- **Contamination Check**: Adapter and rRNA content

---

## Step 4: Trim Adapters (2 min)

1. Click **Trimming** in sidebar
2. In **Input Files** tab, select "Use files from Project"
3. In **Trimming Settings** tab:
   - Select your adapter preset (e.g., Illumina TruSeq)
   - Set min length: 18 nt
   - Set max length: 35 nt
4. Click **Run Trimming**

**Expected:** >70% reads passing filters

---

## Step 5: Set Up Reference (2 min)

1. Click **Databases** in sidebar
2. In **miRBase** tab:
   - Select your organism (e.g., `hsa` for human)
   - Click **Download from miRBase**
3. After download completes, click **Build Bowtie Index**
4. Wait for index to build (1-2 min)

---

## Step 6: Align Reads (2 min)

1. Click **Alignment** in sidebar
2. Verify reference index is selected (from step 5)
3. In **Run Alignment** tab:
   - Select "Use trimmed files from previous step"
   - Keep default settings
4. Click **Start Alignment**

**Expected:** >50% alignment rate

---

## Step 7: Post-Alignment QC (1 min)

1. Click **Post-Align QC** in sidebar
2. Click **Run Post-Alignment QC**

**Check:**
- Mapping statistics
- 5' nucleotide bias (T/U enrichment = good miRNA library)
- RNA composition

---

## Next Steps

Continue your analysis:
- **Counting**: Generate count matrix for all samples
- **DE Analysis**: Compare conditions (need replicates)
- **Targets**: Predict miRNA targets
- **GO/Pathway**: Functional enrichment

---

## Quick Reference

### Recommended Settings for miRNA

| Parameter | Value |
|-----------|-------|
| Min read length | 18 nt |
| Max read length | 25 nt |
| Mismatches | 0-1 |
| Multi-map limit | 10 |

### Expected Results

| Metric | Good | Check |
|--------|------|-------|
| Trimming pass | >70% | <50% |
| Alignment rate | >70% | <30% |
| 5' U bias | >50% | <30% |

---

## Common Issues

**0 reads after trimming?**
- Try different adapter preset
- Check if files are already trimmed

**0% alignment?**
- Build the Bowtie index first
- Verify organism matches your samples

**Low alignment (<30%)?**
- Check for rRNA contamination
- Verify adapters were trimmed

---

## Need More Help?

- Click **Help** in sidebar for full documentation
- Check the Troubleshooting section
- Review the FAQ

---

*sRNAtlas - Comprehensive Small RNA-seq Analysis*
