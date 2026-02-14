# Tutorial 9 – Visualization and Reporting: Publication-Quality Figures and HTML Reports

In this final tutorial, you will learn how to generate publication-ready figures, customize visualizations, and export comprehensive HTML reports that document your entire analysis. This tutorial assumes you have completed at least Tutorials 1–5 (the core pipeline through differential expression).

***

## 9.1 Overview of visualization capabilities

sRNAtlas provides interactive visualizations at every stage of the analysis, all generated using **Plotly** for high-quality, interactive graphics. Key visualization types include:

- **QC plots**: Length distributions, quality profiles, GC content, contamination heatmaps
- **Alignment plots**: Mapping rate bar charts, 5′-nucleotide composition, MAPQ distributions
- **DE plots**: Volcano plots, MA plots, heatmaps, PCA/clustering
- **Enrichment plots**: Dot plots, bar charts, network diagrams (if implemented)
- **Custom plots**: User-defined visualizations via the plotting interface

All plots can be:

- **Viewed interactively** in the browser (zoom, pan, hover for details)
- **Exported as static images** (PNG, SVG, PDF) for publications
- **Exported as HTML** (preserves interactivity for supplementary materials)
- **Customized** (colors, labels, axes, themes)

***

## 9.2 Customizing plots for publication

### 9.2.1 Accessing the plot customization interface

When a plot is displayed in any module (QC, Alignment, DE, Enrichment), you will see a **"Customize plot"** expander below the figure. Click it to reveal customization options.

**General customization options available for most plots:**

1. **Dimensions**:
   - Width (pixels): 600–1200 (default: 800)
   - Height (pixels): 400–800 (default: 600)
   - Aspect ratio lock: Maintain proportions when resizing

2. **Title and labels**:
   - Plot title (e.g., "Volcano Plot: Drought vs. Control")
   - X-axis label (e.g., "log₂ Fold Change")
   - Y-axis label (e.g., "−log₁₀ Adjusted P-value")
   - Font size for title (12–24 pt)
   - Font size for axis labels (10–18 pt)

3. **Color scheme**:
   - **Single-color plots** (bar charts, line plots): Choose from palette or enter hex code
   - **Multi-category plots** (heatmaps, PCA): Select color scale (Viridis, RdBu, Spectral, custom)
   - **DE plots**: Customize colors for up-regulated, down-regulated, and non-significant features

4. **Legend**:
   - Position (top-right, bottom-left, outside, etc.)
   - Font size (8–14 pt)
   - Show/hide legend

5. **Grid and background**:
   - Show/hide grid lines
   - Background color (white, light gray, transparent)

6. **Theme presets**:
   - **Default** (Plotly default with sRNAtlas colors)
   - **Minimal** (white background, no grid, thin axes)
   - **Nature style** (white background, grid, serif font)
   - **Dark mode** (dark background, light text)

### 9.2.2 Example: Customizing a volcano plot

**Scenario**: You want a volcano plot for the paper's main figure, following *Nature* style guidelines.

1. Go to **DE Analysis → Results → Volcano Plot**.
2. Click **"Customize plot"**.
3. Set parameters:
   - Width: 1000 px
   - Height: 800 px
   - Title: (leave blank, you'll add it in Illustrator)
   - X-axis label: "log₂ Fold Change (Drought / Control)"
   - Y-axis label: "−log₁₀ Adjusted *P*-value"
   - Up-regulated color: `#d62728` (red)
   - Down-regulated color: `#1f77b4` (blue)
   - Non-significant color: `#7f7f7f` (gray)
   - Alpha (transparency): 0.6
   - Point size: 4
   - Theme: **Nature style**
4. Click **"Apply customization"**.
5. Export:
   - Click **"Export as SVG"** for vector graphics (editable in Illustrator/Inkscape).
   - Or **"Export as PNG"** at 300 DPI for raster (if journal accepts).

**Result**: A clean, publication-ready volcano plot with proper axis labels and colors matching your manuscript style.

### 9.2.3 Example: Customizing a heatmap

**Scenario**: Heatmap of top 50 DE miRNAs, clustered by sample similarity.

1. Go to **DE Analysis → Results → Heatmap**.
2. Select:
   - Number of top features: 50
   - Clustering: Hierarchical (samples and features)
   - Distance metric: Euclidean
   - Linkage method: Average
3. Click **"Customize plot"**:
   - Width: 800 px
   - Height: 1000 px (tall, to accommodate 50 feature names)
   - Color scale: **RdBu (reversed)** (blue = low, red = high)
   - Z-score normalization: **Row** (normalize per feature)
   - Show dendrograms: Yes
   - Cluster row labels font size: 8 pt
4. Export as **PNG (300 DPI)** or **SVG**.

**Tip**: For heatmaps with many features, export as SVG and edit in Illustrator to adjust row label font size and spacing.

***

## 9.3 Plot types and interpretation

### 9.3.1 Quality Control plots

**Length distribution histogram:**

- **What it shows**: Number of reads per read length (nt).
- **Overlays**: Shaded regions for RNA types (miRNA: 18–25 nt, siRNA: 20–24 nt, piRNA: 24–32 nt, etc.).
- **Interpretation**:
  - Strong peak at 21–22 nt → miRNA-enriched library.
  - Peak at 24 nt in plants → hc-siRNAs present.
  - Broad distribution 15–35 nt → degraded RNA or mixed library.
- **Export use**: Figure 1 or Supplementary Figure S1 to demonstrate library quality.

**Quality score per position:**

- **What it shows**: Mean Phred quality score across read positions.
- **Interpretation**:
  - Quality should be ≥30 across the read length.
  - Drop at 3′ end is common (adapter region).
- **Export use**: Supplementary figure to justify trimming parameters.

**GC content distribution:**

- **What it shows**: Histogram of per-read GC%.
- **Expected**: Plant miRNAs typically 40–60% GC.
- **Red flag**: Bimodal distribution may indicate contamination.

**Contamination heatmap:**

- **What it shows**: Presence/absence of adapter signatures and overrepresented sequences.
- **Export use**: QC documentation in methods.

### 9.3.2 Alignment plots

**Mapping rate bar chart:**

- **What it shows**: Percentage of reads aligned per sample.
- **Interpretation**:
  - >70% mapped → good alignment to reference.
  - <50% mapped → wrong reference or poor library quality.
- **Export use**: Main text or supplementary table.

**5′-nucleotide composition (by size class):**

- **What it shows**: Bar chart showing fraction of each nucleotide (A, C, G, U) at the 5′ position, separated by read length.
- **Plant-specific interpretation**:
  - 21-nt reads: 5′-U enrichment (60–80%) → AGO1-loaded miRNAs.
  - 24-nt reads: 5′-A enrichment (60–80%) → AGO4-loaded hc-siRNAs.
- **Export use**: Validation figure showing biological correctness of library (e.g., "5′-nucleotide profiles confirm AGO sorting").

**MAPQ distribution:**

- **What it shows**: Histogram of mapping quality scores.
- **Interpretation**:
  - Peak at high MAPQ (>30) → mostly unique mappers.
  - Uniform distribution → many multi-mappers (expected for miRNAs).

### 9.3.3 Differential expression plots

**Volcano plot:**

- **What it shows**: Scatter plot of log₂ fold change (x-axis) vs. −log₁₀ adjusted *P*-value (y-axis).
- **Points colored by**:
  - Red: significantly up-regulated.
  - Blue: significantly down-regulated.
  - Gray: non-significant.
- **Thresholds**: Vertical lines at ±log₂FC threshold, horizontal line at adjusted *P* threshold.
- **Hover**: miRNA name, log₂FC, *P*-value, adjusted *P*-value.
- **Export use**: Main figure in results section.

**MA plot:**

- **What it shows**: Scatter plot of mean normalized expression (x-axis) vs. log₂ fold change (y-axis).
- **Purpose**: Shows relationship between expression level and fold change; helps identify outliers.
- **Export use**: Supplementary figure or main text if volcano plot is in supplement.

**PCA plot:**

- **What it shows**: Principal component analysis of normalized counts, with samples colored by condition.
- **Interpretation**:
  - Samples cluster by condition → biological signal present.
  - Outlier samples → potential batch effects or mislabeling.
  - PC1 separates conditions → strong differential expression signal.
- **Export use**: Main figure showing experimental design validation, or supplementary figure.

**Heatmap of top DE features:**

- **What it shows**: Hierarchical clustering of top *N* DE features (rows) across samples (columns).
- **Color**: Expression level (z-score normalized per feature).
- **Dendrograms**: Show sample similarity and feature co-expression patterns.
- **Export use**: Main figure showing expression patterns, or supplementary figure.

### 9.3.4 Enrichment plots

**Dot plot:**

- **What it shows**: Enriched GO terms or KEGG pathways (y-axis) vs. adjusted *P*-value (x-axis).
- **Point size**: Number of genes in the term.
- **Point color**: Adjusted *P*-value (gradient from light to dark).
- **Interpretation**: Top terms are most significantly enriched.
- **Export use**: Main figure or supplementary figure showing functional interpretation.

**Bar chart:**

- **What it shows**: Top *N* enriched terms (y-axis) vs. −log₁₀ adjusted *P*-value (x-axis).
- **Simpler than dot plot** but lacks gene count information.
- **Export use**: Main text figure for simplicity.

***

## 9.4 Exporting plots in publication formats

### 9.4.1 Export formats

sRNAtlas supports the following export formats:

| Format | Use case | Advantages | Disadvantages |
|---|---|---|---|
| **PNG** (300 DPI) | Print publication (raster) | Widely accepted, exact reproduction | Not editable, file size |
| **SVG** | Vector graphics | Editable in Illustrator/Inkscape, scalable | Some journals convert to raster |
| **PDF** | Vector graphics | Preferred by some journals, scalable | Requires conversion for web |
| **HTML** | Supplementary materials | Fully interactive, hover tooltips | Not suitable for print |
| **JSON** | Plot data | Re-create plot in R/Python | Requires programming |

### 9.4.2 Recommended export workflow

**For main figures:**

1. Generate plot in sRNAtlas with desired customization.
2. Export as **SVG**.
3. Open in **Adobe Illustrator** or **Inkscape**:
   - Adjust font sizes to match manuscript style (typically 7–9 pt for axis labels, 8–10 pt for axis titles).
   - Fine-tune spacing and alignment.
   - Add panel labels (A, B, C) if creating a multi-panel figure.
4. Export from Illustrator as **PDF** or **TIFF (300 DPI)** for submission.

**For supplementary figures:**

- Export as **PNG (300 DPI)** if raster is acceptable.
- Export as **HTML** if you want reviewers/readers to interact with the plot (e.g., hover over points in a volcano plot to see miRNA names).

**For data sharing:**

- Export plot data as **JSON** or **CSV** (available via "Export data" button below most plots).
- Include in supplementary data tables.

### 9.4.3 Batch export

If you have many plots to export (e.g., QC for 20 samples, volcano plots for 5 comparisons), use the **Batch Export** feature:

1. Go to **Reports → Batch Export**.
2. Select plot types to export:
   - ☐ QC length distributions (all samples)
   - ☐ Volcano plots (all comparisons)
   - ☐ MA plots (all comparisons)
   - ☐ Heatmaps
   - ☐ PCA
3. Select export format (PNG, SVG, or both).
4. Set naming convention (e.g., `volcano_[comparison]_[date].svg`).
5. Click **"Export all"**.
6. Downloads a ZIP file with all plots organized by type.

***

## 9.5 Generating HTML reports

The **Reports** module creates comprehensive, self-contained HTML documents that summarize the entire analysis. These reports are ideal for:

- **Sharing with collaborators** who don't have sRNAtlas installed.
- **Archiving** the analysis for future reference.
- **Supplementary materials** for publications (include the HTML report as supplementary file).

### 9.5.1 Report structure

A standard sRNAtlas HTML report includes:

1. **Title page**: Project name, analysis date, sRNAtlas version.
2. **Summary**: Key statistics (number of samples, reads, DE features).
3. **QC section**:
   - Table of QC metrics per sample (total reads, mapped reads, quality scores).
   - Length distribution plots (embedded as PNG or SVG).
   - Contamination scan results.
4. **Alignment section**:
   - Mapping rate table.
   - 5′-nucleotide composition plots.
   - Post-alignment QC metrics.
5. **Differential expression section**:
   - Comparison summary table (number of up/down-regulated features).
   - Volcano plots (interactive if HTML, static if PDF export).
   - Top DE features table (sortable).
6. **Enrichment section**:
   - Enriched GO terms and KEGG pathways (if enrichment was run).
   - Dot plots or bar charts.
7. **Methods**: Auto-generated methods text with all parameters used.
8. **Appendix**: Software versions, references, provenance log.

### 9.5.2 Generating a report

**Step 1: Go to the Reports module**

Navigate to **Reports** in the left sidebar.

**Step 2: Configure report contents**

Select sections to include:

- ☑ Title page
- ☑ Summary statistics
- ☑ Quality control
- ☑ Alignment and mapping
- ☑ Differential expression
- ☑ Target prediction (if available)
- ☑ Enrichment analysis (if available)
- ☑ Methods (auto-generated)
- ☑ Provenance log

**Step 3: Set report options**

- **Report title**: E.g., "*Medicago truncatula* Drought Stress Small RNA Analysis"
- **Author**: Your name or lab name
- **Date**: Auto-filled (current date) or custom
- **Include interactive plots**: Yes (HTML with embedded Plotly) or No (static PNG)
- **Table of contents**: Yes (clickable navigation) or No
- **Theme**: Light (default) or Dark

**Step 4: Generate report**

Click **"Generate Report"**. The process takes 10–60 seconds depending on the number of samples and plots.

**Step 5: Preview and download**

- **Preview**: Opens the HTML report in a new browser tab. Review all sections.
- **Download**: Click "Download HTML Report" to save the file (e.g., `sRNAtlas_Report_20260213.html`).

### 9.5.3 Customizing the report

**Adding custom sections:**

In the **Reports** module, there is an option to add **custom Markdown sections**:

1. Click **"Add custom section"**.
2. Enter a section title (e.g., "Biological Context").
3. Write Markdown content in the text box:

   ```markdown
   ## Biological Context
   
   This experiment investigates the role of microRNAs in drought stress response in *Medicago truncatula*. Previous studies have shown that miR398 and miR408 are key regulators of oxidative stress responses .
   
   ### Hypotheses
   - Drought stress will up-regulate miR398, leading to down-regulation of copper/zinc superoxide dismutases (CSD1, CSD2).
   - Novel stress-specific miRNAs may be identified from unaligned reads.
   
   ### References
    Sunkar et al. (2006). Novel and stress-regulated microRNAs in Arabidopsis. *Plant Cell*, 18:2001-2017.
   ```

4. Click **"Add"**. The custom section will appear in the report after the Summary and before QC.

**Editing the auto-generated Methods:**

The auto-generated Methods text is based on the actual parameters used in your analysis. You can edit it before finalizing the report:

1. In the Reports module, scroll to **"Auto-generated Methods"**.
2. Click **"Edit"**.
3. Modify the text (e.g., add organism-specific details, clarify experimental design).
4. Click **"Save"**. The edited version will be used in the report.

### 9.5.4 Exporting reports as PDF

HTML reports are great for interactivity, but some journals require PDF supplements.

**Method 1: Print to PDF from browser**

1. Open the HTML report in Chrome or Firefox.
2. Press `Ctrl+P` (Windows/Linux) or `Cmd+P` (Mac).
3. Select "Save as PDF" as the destination.
4. Adjust settings:
   - Layout: Portrait (or Landscape for wide tables)
   - Margins: Default
   - Background graphics: On (to preserve plot backgrounds)
5. Click "Save".

**Limitation**: Interactive plots lose interactivity in PDF (rendered as static images).

**Method 2: Use the built-in PDF export**

If sRNAtlas has a **"Export as PDF"** button in the Reports module:

1. Click **"Export as PDF"** instead of "Generate Report".
2. sRNAtlas will render all plots as static images and compile a PDF using a PDF library (e.g., ReportLab or WeasyPrint).
3. Download the PDF.

***

## 9.6 Creating multi-panel figures

For publications, you often need **composite figures** with multiple panels (e.g., Figure 1: A = QC length distribution, B = volcano plot, C = heatmap, D = enrichment).

**Workflow:**

1. Export individual panels from sRNAtlas as **SVG** or **PNG (300 DPI)**.
2. Use a vector graphics editor (Illustrator, Inkscape, or even PowerPoint) to assemble panels:
   - Import all panels as separate images.
   - Arrange in a grid or custom layout.
   - Add panel labels (A, B, C, D) in bold, consistent font (e.g., Arial Bold 12 pt).
   - Ensure consistent axis label font sizes across panels.
3. Export the composite figure as **TIFF (300 DPI)** or **PDF**.

**sRNAtlas built-in multi-panel export (if available):**

Some versions of sRNAtlas may include a **"Create composite figure"** tool:

1. Go to **Reports → Composite Figures**.
2. Select plots to include (e.g., QC length, volcano, heatmap, enrichment).
3. Choose layout: 2×2 grid, 1×4 horizontal, custom.
4. Set panel labels: A, B, C, D (auto-positioned).
5. Click **"Generate composite"**.
6. Export as **SVG** or **PNG**.

***

## 9.7 Data export for external tools

If you want to create plots in R (ggplot2) or Python (matplotlib/seaborn) using sRNAtlas 

### 9.7.1 Exporting count matrices

1. Go to **Counting → Export**.
2. Select format:
   - **CSV**: Raw counts, suitable for any tool.
   - **DESeq2-compatible**: Includes metadata in the format expected by DESeq2 in R.
   - **Normalized counts**: VST or rlog-transformed counts from pyDESeq2.
3. Click **"Download"**.

### 9.7.2 Exporting DE results

1. Go to **DE Analysis → Results → Export**.
2. Select comparison (or "Export all").
3. Format: **CSV** (includes log₂FC, *P*-value, adjusted *P*-value, baseMean).
4. Click **"Download"**.

### 9.7.3 Exporting enrichment results

1. Go to **Enrichment → Results → Export**.
2. Format: **CSV** or **JSON**.
3. Contains: Term ID, term name, adjusted *P*-value, gene list.

### 9.7.4 Example: Re-creating a volcano plot in R

After exporting DE results as `Drought_vs_Control_DE.csv`:

```r
library(ggplot2)

# Load data
de <- read.csv("Drought_vs_Control_DE.csv")

# Add significance column
de$Significance <- ifelse(de$padj < 0.05 & abs(de$log2FoldChange) > 0.585,
                          ifelse(de$log2FoldChange > 0, "Up", "Down"),
                          "NS")

# Volcano plot
ggplot(de, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray")) +
  theme_classic() +
  labs(title = "Drought vs. Control",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-value") +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")

ggsave("volcano_custom.pdf", width = 6, height = 5)
```

This gives you full control over plot aesthetics using ggplot2.

***

## 9.8 Summary

This tutorial demonstrated how to:

- **Customize plots** in sRNAtlas to match publication style guidelines (colors, fonts, themes).
- **Export plots** in multiple formats (PNG, SVG, PDF, HTML) suitable for print and supplementary materials.
- **Generate HTML reports** that document the entire analysis with embedded plots, tables, and auto-generated methods text.
- **Create multi-panel figures** by combining individual plots in vector graphics editors.
- **Export data** for external visualization tools (R, Python) when more advanced customization is needed.

With these skills, you can produce publication-ready figures and comprehensive reports directly from sRNAtlas, streamlining the path from data analysis to manuscript submission.

***

**End of Tutorial Series**

You have now completed the full sRNAtlas tutorial series:

- **Tutorial 1**: Project setup and sample management
- **Tutorial 2**: Reference database management (miRBase, RNAcentral)
- **Tutorial 3**: Quality control and adapter trimming
- **Tutorial 4**: Alignment and post-alignment QC
- **Tutorial 5**: Read counting and differential expression
- **Tutorial 6**: Target prediction and functional enrichment
- **Tutorial 7**: Batch processing and project save/reload
- **Tutorial 8**: Advanced analysis (hierarchical alignment, isomiRs, novel miRNAs)
- **Tutorial 9**: Visualization and reporting

You are now equipped to perform complete, rigorous small RNA-seq analyses from raw FASTQ files to biological interpretation, all within the sRNAtlas platform.

**Next steps:**

- Apply sRNAtlas to your own datasets.
- Explore advanced parameters and custom references.
- Contribute to the project on GitHub: https://github.com/osman12345/sRNAtlas
- Cite sRNAtlas in your publications (see manuscript for citation details).

Happy analyzing! 🧬📊

Sources
