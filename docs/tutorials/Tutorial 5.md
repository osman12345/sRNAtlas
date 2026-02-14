# Tutorial 5: Differential Expression Analysis

How to configure contrasts, run pyDESeq2, and interpret volcano plots, heatmaps, and PCA inside sRNAtlas.

---

## Overview

| | |
|---|---|
| **Learning objectives** | Set up DE design and contrasts; run pyDESeq2; interpret volcano plots, MA plots, heatmaps, and PCA; export results for downstream work |
| **Prerequisites** | Tutorials 2–4 completed; a count matrix generated in **Counting**; sample metadata defined in **Project** |
| **Estimated time** | 30–45 minutes |
| **Difficulty** | Intermediate |

---

## Part 1: Before You Start

Differential expression (DE) in sRNAtlas uses **pyDESeq2**, a Python implementation of the DESeq2 method, so all the usual assumptions apply:

- Counts must be **raw integers** (no CPM/TPM/logs).
- There should be at least **2–3 biological replicates per group** for reliable inference.
- Large global shifts (e.g., many miRNAs up in one condition) are handled by DESeq2’s **median-ratio normalization**.

In this tutorial, we assume:

- You have run **Counting** and have a count matrix in `st.session_state.count_matrix`.
- You have defined your experimental groups in **Project → Metadata** (e.g. `Condition` column with `Control` vs `Treatment`).

---

## Part 2: Setting Up DE Analysis

### 2.1 Open the DE Analysis module

1. Click **DE Analysis** in the sidebar.
2. You should see three main sections:
   - **Setup**
   - **Run Analysis**
   - **Results & Visualizations**

If sRNAtlas cannot find a count matrix or metadata, it will show a warning. Return to **Project** and **Counting** to fix that before continuing.

### 2.2 Choose design factor and metadata

In the **Setup** panel:

1. **Select count matrix source**
   - If you followed previous tutorials, choose **“Use current project count matrix”**.
   - Optionally, you can upload an external count matrix (CSV) instead.

2. **Select metadata**
   - If you defined metadata in **Project**, choose **“Use project metadata”**.
   - Make sure the **SampleID** column matches the column names in your count matrix.

3. **Design factor**
   - Pick the main categorical variable you want to test (e.g. `Condition`).
   - This is the factor that will enter the DESeq2 design formula.

> **Tip:** For simple two-group comparisons (Control vs Treatment), keep the design as a single factor. For more complex designs (e.g. `~ Batch + Condition`), consider starting simple in the GUI and documenting any extra modeling in the paper.

---

## Part 3: Defining Comparisons

### 3.1 Single comparison (two groups)

In **Setup → Comparisons**:

1. Click **“Add comparison”**.
2. Choose:
   - **Factor:** usually same as design factor (`Condition`).
   - **Reference level:** e.g. `Control`.
   - **Test level:** e.g. `Treatment`.
3. Give the comparison a **name**, for example: `Treatment_vs_Control`.

This sets up a contrast equivalent to:

\[
\text{log}_2(\text{Treatment} / \text{Control})
\]

### 3.2 Multiple pairwise comparisons

You can add multiple comparisons, for example in a three-condition design:

- `TreatmentA_vs_Control`
- `TreatmentB_vs_Control`
- `TreatmentB_vs_TreatmentA`

Add each one using the **Add comparison** button.

### 3.3 Global FDR option

In **Setup → Multiple testing**:

- Toggle **“Global FDR across comparisons”** if you want a joint FDR correction across *all* tests in *all* contrasts.
- Otherwise, each comparison is adjusted independently (DESeq2 standard).

> **Recommendation:** For a small number of planned comparisons (e.g., 2–3), per-contrast FDR is usually fine. For many exploratory contrasts, consider enabling global FDR.

---

## Part 4: Running pyDESeq2

### 4.1 Launch analysis

1. Go to the **Run Analysis** section.
2. Check that:
   - The **count matrix** preview looks correct (rows = features, columns = samples).
   - Metadata looks correct (sample IDs, group labels).
3. Choose:
   - **FDR threshold** (*alpha*), default 0.05.
   - Optional minimum base mean / counts if the UI exposes it.

4. Click **“🚀 Run DE Analysis”**.

The app will:

- Prepare counts and metadata.
- Create a `DeseqDataSet` object.
- Run `dds.deseq2()` for each comparison.
- Extract results into a table with:
  - log2FoldChange
  - baseMean
  - pvalue
  - padj (or global_FDR if global FDR was selected)

### 4.2 Checking for errors

If something goes wrong, common issues include:

- **Sample name mismatch**: no overlap between count columns and metadata SampleID.
- **Too few samples per group**: pyDESeq2 may complain or produce unstable estimates.
- **All zeros** for a feature across samples: DESeq2 automatically filters those.

Fix the underlying issue (usually in **Project → Metadata** or **Counting**) and rerun.

---

## Part 5: Inspecting DE Results

Once the run completes, sRNAtlas stores results in `st.session_state.de_results`.

### 5.1 Results table

In the **Results** tab:

1. Choose a **comparison** from the dropdown (e.g. `Treatment_vs_Control`).
2. You’ll see a table with columns such as:
   - `baseMean`
   - `log2FoldChange`
   - `lfcSE`
   - `stat`
   - `pvalue`
   - `padj` (or `global_FDR`)

3. Use the filters:
   - FDR threshold (e.g. 0.05).
   - |log2FC| threshold (e.g. > 1).

4. You can **download** the full table as CSV for external analysis.

> **Interpretation:** For miRNAs, even modest log2FC (~0.6, i.e. 1.5-fold) can be biologically meaningful, especially if the miRNA is abundant.

---

## Part 6: Visualizations

Click **Visualizations** in the DE module to explore:

### 6.1 Volcano plot

1. Select **“Volcano Plot”**.
2. Pick the comparison.
3. Adjust:
   - **LFC Threshold** (horizontal cutoffs).
   - **FDR Threshold** (vertical cutoff).

The plot shows:

- x-axis: log2FoldChange.
- y-axis: −log10(p-value).
- Points coloured by significance (up, down, not significant).

Use the **download buttons** to export PNG/SVG for figures.

### 6.2 MA plot

1. Select **“MA Plot”**.
2. Choose a comparison.

The MA plot shows:

- x-axis: mean expression (baseMean).
- y-axis: log2FoldChange.

DE features stand out from the cloud. Low-count features tend to show more scatter.

### 6.3 PCA

1. Choose **“PCA”**.
2. Select a metadata column to colour by (e.g. `Condition`).

The PCA is computed on **normalized counts**:

- Samples should cluster by biological condition if the signal is strong.
- If they cluster by batch or library prep instead, you likely have technical confounding.

### 6.4 Sample correlation

1. Choose **“Sample Correlation”**.
2. Select correlation method (Spearman or Pearson).

You get a correlation heatmap of samples:

- Tight blocks along the diagonal indicate well‑correlated replicates.
- Outliers appear as low‑correlation rows/columns.

### 6.5 Heatmap of top DE features

1. Choose **“Heatmap”**.
2. Select a comparison.
3. Set:
   - FDR threshold.
   - Number of top features (e.g. 50–100).

The module:

- Selects significant features.
- Takes the top N by FDR.
- Z‑scores them across samples.
- Plots a hierarchical cluster heatmap.

---

## Part 7: Practical Example

Suppose you have:

- 3 controls (C1–C3).
- 3 treatments (T1–T3).
- One comparison: `Treatment_vs_Control` on factor `Condition`.

A typical workflow:

1. Run DE with FDR 0.05, no global FDR.
2. In **Results**, filter for:
   - `padj < 0.05`
   - `|log2FoldChange| > 0.585` (~1.5‑fold).
3. Inspect the **volcano plot**:
   - Identify a cluster of strongly upregulated miRNAs (log2FC ~ 2–3).
4. Use **PCA** coloured by `Condition`:
   - Confirm that C1–C3 cluster together and T1–T3 cluster separately.
5. Use the **Heatmap**:
   - Visualize the top 50 miRNAs to see if they separate the two groups coherently.
6. Export:
   - DE table CSV.
   - Volcano and heatmap PNG for figures.

---

## Part 8: Connecting DE to Targets and Enrichment

Tutorials 6 and 7 will show how to:

- Take your list of **significant miRNAs** and run **target prediction** (psRNATarget for plants, TargetScan/miRanda for animals if you integrate them).
- Use **GO & Pathway Enrichment** to interpret target gene lists.

In sRNAtlas, these modules can read DE results directly, so you don’t have to re‑upload anything.

---

## Checklist

Before moving on:

- [ ] Design factor and levels correctly defined in **Project**.
- [ ] Comparisons specified (ref vs test) and named clearly.
- [ ] DE run completed without errors.
- [ ] Volcano, MA, PCA, and heatmap inspected.
- [ ] DE tables exported for record‑keeping.

---

*Save this file as `docs/tutorials/05-differential-expression.md` in your repository.*
