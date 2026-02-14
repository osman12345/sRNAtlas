# Tutorial 6 – Target Prediction and Functional Enrichment

In this tutorial you will use the **Target Prediction** and **GO/Pathway Enrichment** modules to go from a list of differentially expressed miRNAs to predicted target genes and enriched biological processes.

We assume that you have already:

- Run the pipeline through **DE Analysis** (Tutorial 5) and have at least one comparison with significantly up‑ or down‑regulated miRNAs, and  
- Configured a gene identifier system that is compatible with g:Profiler (e.g. TAIR IDs, Medicago gene IDs, or Ensembl IDs for animal species).

If you are starting from an external list of miRNAs and targets, you can still follow this tutorial, starting from Step 2.

***

## 6.1 Overview of the Target Prediction module

Open the **Target Prediction** page from the left-hand navigation. The module has four logical parts:

1. **miRNA input**  
   - Pull significant miRNAs directly from the DE module, **or**  
   - Upload your own list of miRNAs (with sequences or IDs).

2. **Target prediction method**  
   - **psRNATarget (API)** – recommended for **plant** miRNAs; remote call to the psRNATarget server.  
   - **Local seed matching** – fast, purely local heuristic screen based on seed complementarity.

3. **Target gene set management**  
   - View and filter predicted targets per miRNA or per comparison.  
   - Export target lists for downstream analysis.

4. **Hand‑off to enrichment**  
   - Send a set of gene IDs directly into the **GO/Pathway Enrichment** module.

The module is designed so that you can start either from the DE results inside sRNAtlas or from external miRNA lists.

***

## 6.2 Selecting miRNAs as input

### 6.2.1 Using differentially expressed miRNAs from pyDESeq2

1. Go to the **Input** section of the Target Prediction page.
2. Choose **“Use DE results”** as the miRNA source.
3. Select the DE comparison you are interested in (e.g. `Treatment_vs_Control`).
4. Set the thresholds used to define “significant” miRNAs:
   - **Adjusted P‑value (FDR) threshold** – default 0.05  
   - **log₂ fold‑change threshold** – default |log₂FC| > 0.585 (~1.5‑fold)
5. Click **“Load miRNAs from DE”**.

sRNAtlas will:

- Filter the DE results for that comparison using the thresholds.
- Resolve each significant miRNA ID to its **nucleotide sequence** using the active reference (miRBase or custom FASTA).  
- Show a table with columns such as:
  - `miRNA_id`, `Sequence`, `log2FC`, `padj`, `Direction` (up/down).

If any miRNA cannot be matched to a sequence in the current reference, it will be flagged, and you should either:

- Add a reference containing that miRNA to the **Databases** module, or  
- Provide the sequence manually via a CSV upload (next section).

### 6.2.2 Providing your own miRNA list

If you want to predict targets for miRNAs that did not come from the DE module (e.g. published miRNAs or candidates from another tool):

1. In the same **Input** area, choose **“Upload miRNA list”**.
2. Upload a CSV file with at least one of the following schemas:

   **Option A – IDs + sequences (recommended)**  
   ```text
   miRNA_id,sequence
   mtr-miR160a,UGCCUGGCUCCCUGUAUGCCA
   mtr-miR167a,UAGAAGCUGCCAGCAUGAUCU
   ```

   **Option B – IDs only**  
   ```text
   miRNA_id
   mtr-miR160a
   mtr-miR167a
   ```
   In this case sRNAtlas will attempt to resolve sequences from the currently loaded miRBase or custom reference.

3. After upload, check the preview table to confirm that all miRNAs have valid sequences.

***

## 6.3 Choosing a prediction method

In the **Method** section, select one of the available target prediction engines.

### 6.3.1 psRNATarget (API) – recommended for plants

Use this mode when working with **plant** miRNAs and you want established psRNATarget scoring:

1. Set **Prediction method** to **“psRNATarget (API)”**.
2. Configure the psRNATarget parameters:
   - **Expectation threshold** (default: 3.0–4.0)  
     Lower values are stricter. Typical plant settings are 3.0–4.0.
   - **Maximum unpaired energy (UPE)**, **flanking length**, and **seed region** are exposed in the UI if you need fine control. Defaults mirror the psRNATarget 2017 recommended settings.
3. Choose the **target transcript dataset**:
   - For standard species, psRNATarget’s own reference transcript sets are used.  
   - For non‑model organisms, you may need to upload or map your gene IDs carefully in downstream enrichment.
4. Click **“Run psRNATarget”**.

**Notes:**

- This mode sends the miRNA sequences and identifiers to the psRNATarget server and retrieves the predictions as a table.  
- If the server is temporarily unavailable, sRNAtlas shows a clear error message instead of failing silently.

The results table typically includes:

- `miRNA_id`  
- `Target_gene_id`  
- `Target_accession` / transcript ID  
- `Alignment_score` / `Expectation`  
- Target site position and alignment details.

You can filter the table by expectation score and export it as CSV for record‑keeping.

### 6.3.2 Local seed matching – fast, local heuristic

Use this mode when you:

- Want a very fast, purely local screen (no external API), or  
- Are working with organisms where psRNATarget does not provide a curated transcript set.

1. Set **Prediction method** to **“Local seed matching”**.
2. Configure the seed parameters:

   - **Seed region on the miRNA**  
     - **Start position (nt)** – default 2  
     - **End position (nt)** – default 8  
     This defines the canonical “seed” region.
   - **Maximum mismatches in the seed** (default: 0–1)  
     0 is strict (perfect complementarity), 1 allows a single mismatch.
   - Optionally, set whether G:U wobble pairs are allowed (if exposed in the UI).

3. Select or upload the **target transcript FASTA** if required:
   - For many plant species, you may already have a genome or transcript FASTA loaded as a reference in the **Databases** or **Novel miRNA** modules.
   - Otherwise, upload a transcript FASTA containing your gene models.

4. Click **“Run local seed search”**.

The algorithm will:

- Extract the specified seed from each miRNA sequence.
- Scan the target transcript sequences for complementary sites that meet the mismatch criteria.
- Report candidate miRNA–target pairs with basic alignment metrics (seed match position, number of mismatches, strand, etc.).

This mode is intended as a **screening tool**, not a replacement for psRNATarget’s biophysically informed scoring.

***

## 6.4 Inspecting and exporting target lists

After the chosen prediction method finishes, the **Results** section will show one or more tables. The exact layout depends on the implementation, but typically includes:

- A **miRNA‑centric view** – for each miRNA, list all predicted targets.
- A **target‑centric view** – for each gene, list all miRNAs predicted to target it.

Typical columns include:

- `miRNA_id`  
- `miRNA_sequence`  
- `Target_gene_id` (this is what you will pass to enrichment)  
- Method‑specific metrics:
  - psRNATarget: `Expectation`, `UPE`, site position.  
  - Local seed: `Seed_mismatches`, `Site_position`, `Strand`.

You can:

- Use filters (e.g. Expectation ≤ 3.0, seed mismatches ≤ 1).  
- Export the tables as CSV for external use.

When preparing for enrichment, you usually want a **deduplicated list of target gene IDs** for:

- All miRNAs that are **up‑regulated**,  
- All miRNAs that are **down‑regulated**, or  
- A specific subset of miRNAs of interest.

sRNAtlas therefore offers a **“Export gene list”** or **“Send to enrichment”** button that automatically collects unique `Target_gene_id` values passing your thresholds.

***

## 6.5 GO and KEGG enrichment with g:Profiler

Once you have a target gene list, switch to the **GO/Pathway Enrichment** module.

### 6.5.1 Providing the gene list

1. In the **Input** section, either:
   - Click **“Use targets from Target Prediction”** to pull the most recent gene list from the session, **or**  
   - Upload a text/CSV file with one gene identifier per line, matching your organism’s gene ID convention.

2. Check that the gene IDs look reasonable (no unexpected prefixes or truncated IDs).

### 6.5.2 Configuring enrichment parameters

1. Select the **organism code** used by g:Profiler:

   - Example plant codes:
     - *Arabidopsis thaliana* – `athaliana`  
     - *Medicago truncatula* – `mtruncatula`  
     - *Glycine max* – `gmax`  
   - Example animal codes:
     - Human – `hsapiens`  
     - Mouse – `mmusculus`.

2. Choose the **annotation sources**:
   - `GO:BP`, `GO:MF`, `GO:CC`  
   - `KEGG`  
   - Optionally others (Reactome, etc., if exposed).

3. Set the **multiple testing correction**:
   - Default is FDR (Benjamini–Hochberg), with a typical adjusted *P* threshold of 0.05.

4. Optionally specify a **background gene set**:
   - By default, g:Profiler uses the genome‑wide background.
   - You can provide your own background list (e.g. all expressed genes in your experiment) via upload or by selecting “use all genes detected in this dataset” if that option is exposed.

5. Click **“Run enrichment”**.

The results table includes:

- `Term ID` (e.g. GO:0009733)  
- `Term name`  
- `Source` (GO:BP, KEGG, etc.)  
- `Adjusted P‑value`  
- `Term size`, `Intersection size`, and associated genes.

### 6.5.3 Visualizing and exporting enrichment results

sRNAtlas provides:

- A **dot plot** summarizing the top enriched terms (GO and/or KEGG), implemented via Plotly.  
- Interactive filters (by adjusted *P* value, term size, and category).  
- Export options:
  - CSV of enriched terms  
  - PNG/SVG of plots.

This closes the loop from:

> Differentially expressed miRNAs → predicted targets → GO/KEGG terms enriched among those targets.

***

## 6.6 Recommended workflows

### 6.6.1 Plant miRNA workflow (Medicago, Arabidopsis, etc.)

1. Run DE analysis to identify differentially expressed miRNAs.  
2. In Target Prediction:
   - Source miRNAs from DE comparison.  
   - Use **psRNATarget (API)** with expectation 3.0–3.5.  
3. Filter targets by expectation threshold.  
4. In Enrichment:
   - Use appropriate plant organism code (e.g. `mtruncatula`).  
   - Run GO:BP and KEGG enrichment.  
5. Interpret enriched terms in the context of miRNA direction (e.g. up‑regulated miRNAs with enriched stress‑response targets).

### 6.6.2 Non‑model or custom workflow

1. If psRNATarget does not provide a reference for your species, build or obtain a **transcript FASTA** with stable gene identifiers.  
2. In Target Prediction:
   - Upload miRNAs (with sequences).  
   - Use **Local seed matching** with a strict seed (2–8 nt, ≤1 mismatch).  
3. Inspect and manually filter the resulting miRNA–target pairs.  
4. In Enrichment:
   - Ensure g:Profiler supports your gene ID format; if not, convert IDs upstream or use an alternative enrichment tool.  
   - Run at least GO:BP enrichment with FDR correction.

***

## 6.7 Summary

This tutorial showed how to:

- Feed significant miRNAs from pyDESeq2 into the **Target Prediction** module.
- Run target prediction either via **psRNATarget (API)** for plants or **Local seed matching** for fast, local screening.
- Export predicted target gene sets and send them directly to the **GO/Pathway Enrichment** module.
- Obtain enriched GO and KEGG pathways associated with miRNA target genes.

Together with the earlier tutorials (trimming, QC, alignment, counting, DE), this completes the end‑to‑end workflow from **raw FASTQ files to biological interpretation** inside sRNAtlas.
