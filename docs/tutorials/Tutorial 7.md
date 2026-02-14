# Tutorial 7 – Batch Processing and Project Management

In this tutorial you will learn how to automate the sRNAtlas pipeline for multiple samples and how to save, reload, and share entire analysis projects. This is essential for reproducibility and for working with datasets containing dozens or hundreds of samples.

We assume that you have already worked through Tutorials 1–6 and are familiar with the individual modules. This tutorial focuses on *orchestrating* those modules efficiently.

***

## 7.1 Overview of batch mode

The **Batch Processing** module allows you to:

- **Queue** a sequence of pipeline stages (trimming → alignment → counting → DE) across all samples in your project.
- **Monitor** progress with real-time logs and a visual progress tracker.
- **Resume** interrupted runs without re-executing completed steps.
- **Export** a complete provenance record (tool versions, parameters, file checksums) for reproducibility.

Batch mode is particularly useful when:

- You have 10+ samples and don't want to manually click through each stage.
- You want to run the pipeline overnight or on a compute cluster (via the command-line interface).
- You need to guarantee that all samples receive identical parameter settings.

***

## 7.2 Setting up a batch run

### 7.2.1 Preparing your project structure

Before entering batch mode, ensure that:

1. **Sample metadata is complete**  
   - Go to the **Project** module and verify that all samples are listed with their corresponding FASTQ files, conditions, and any other experimental factors.
   - Use the "Import metadata CSV" option if you have many samples (see Tutorial 1 for CSV format).

2. **References are indexed**  
   - Go to the **Databases** module and ensure that your miRBase/RNAcentral references (or custom FASTA files) are indexed.
   - Batch mode will not automatically build indices; if an index is missing, the alignment step will fail.

3. **Parameter presets are selected**  
   - For trimming: choose your adapter preset (e.g. Illumina TruSeq Small RNA).
   - For alignment: decide on the number of mismatches (`-v`), maximum alignments per read (`-k`), and counting mode (unique vs. fractional).
   - For DE: define your experimental design (main factor, comparisons).

   These parameters will be applied uniformly across all samples in the batch run.

### 7.2.2 Navigating to the Batch module

Open the **Batch Processing** page from the left-hand sidebar. The interface is divided into three sections:

1. **Pipeline Configuration** – select which stages to run and set global parameters.
2. **Job Queue** – view the list of pending and completed tasks.
3. **Execution** – start, pause, or cancel the batch run.

***

## 7.3 Configuring the pipeline stages

### 7.3.1 Selecting stages to execute

In the **Pipeline Configuration** section, you will see a checklist of available stages:

- ☐ Quality Control (pre-trim)
- ☐ Adapter Trimming
- ☐ Quality Control (post-trim)
- ☐ Alignment
- ☐ Post-Alignment QC
- ☐ Read Counting
- ☐ Differential Expression

Check the boxes for the stages you want to include. Typical scenarios:

- **Full pipeline from raw FASTQ**: Check all stages.
- **Re-run only DE with different thresholds**: Uncheck everything except "Differential Expression" (assumes count matrix already exists).
- **Add a secondary alignment**: Check "Alignment" and "Counting" only, and select a different reference (e.g. RNAcentral tRNA after the initial miRBase run).

### 7.3.2 Setting global parameters

For each selected stage, a collapsible parameter panel appears below the checklist. Examples:

**Trimming parameters:**
- Adapter preset (or custom sequence)
- Minimum/maximum length (default: 18–35 nt)
- Quality threshold (default: Phred ≥ 20)

**Alignment parameters:**
- Reference database (select from dropdown)
- Number of mismatches (`-v 1` default)
- Maximum alignments per read (`-k 10` default)
- Counting mode: `unique_only`, `fractional`, or `all`

**DE parameters:**
- Main factor (e.g. `Treatment`)
- Comparisons (e.g. `Treatment_vs_Control`, `Stress_vs_Control`)
- FDR threshold (default: 0.05)
- Log₂ fold-change threshold (default: 0.585)
- Global FDR correction (on/off)

**Important**: These parameters are applied to *all samples* in the batch. If you need sample-specific settings (e.g. different adapters for different libraries), you must either:

- Run batch mode separately for each subset of samples, or
- Use the individual modules interactively for those samples.

### 7.3.3 Output directory structure

Specify an **output root directory** where batch results will be saved. sRNAtlas will create subdirectories for each stage:

```
output_root/
├── trimming/
│   ├── sample1_trimmed.fastq.gz
│   ├── sample2_trimmed.fastq.gz
│   └── ...
├── alignment/
│   ├── sample1.bam
│   ├── sample1.bam.bai
│   └── ...
├── counts/
│   └── count_matrix.csv
├── de_results/
│   ├── Treatment_vs_Control.csv
│   └── ...
├── qc_reports/
│   └── ...
└── provenance.yaml
```

This structure is compatible with downstream tools and can be archived for long-term storage.

***

## 7.4 Running the batch pipeline

### 7.4.1 Starting the job queue

Once parameters are set:

1. Click **"Validate Configuration"** to check for missing references, invalid file paths, or parameter conflicts.
   - sRNAtlas will display warnings or errors (e.g. "Reference index not found for sample X").
   - Fix any issues before proceeding.

2. Click **"Start Batch Run"** to begin execution.

### 7.4.2 Monitoring progress

The **Job Queue** section shows:

- A table of all samples and their current stage (e.g. "Trimming in progress", "Alignment queued", "Counting complete").
- A progress bar for each stage.
- Real-time log output in an expandable text box.

Example log entry:

```
[2026-02-13 21:15:32] INFO: Starting trimming for sample1
[2026-02-13 21:15:45] INFO: Cutadapt: 1,234,567 reads processed, 98.2% retained
[2026-02-13 21:15:50] INFO: Trimming complete for sample1
[2026-02-13 21:15:51] INFO: Starting alignment for sample1
```

If a stage fails for a particular sample, the log will display the error message (e.g. "Bowtie index not found"). The batch run will:

- Mark that sample as "Failed" for that stage.
- Continue processing other samples (default behavior), or
- Halt the entire batch if you selected "Stop on first error" mode.

### 7.4.3 Pausing and resuming

- Click **"Pause"** to temporarily suspend the batch run (e.g. to free up CPU for other tasks).
- Click **"Resume"** to continue from the last completed step.
- Completed stages are *not* re-run; sRNAtlas checks for existing output files and skips them unless you explicitly select "Overwrite existing results".

### 7.4.4 Reviewing results

After the batch run completes, navigate to the **Reports** module or inspect the output directory manually. Key files to check:

- `count_matrix.csv` – contains read counts for all samples × all features.
- `de_results/` – CSV files for each comparison.
- `qc_reports/` – HTML summaries of QC metrics per sample.
- `provenance.yaml` – complete record of software versions, parameters, and file checksums.

***

## 7.5 Project persistence: saving and reloading

### 7.5.1 Saving a project

The **Project** module's "Save Project" function writes the entire analysis session to disk as a compressed archive (`.srnatlasproject` file). This includes:

- Sample metadata (names, conditions, file paths)
- All intermediate results (trimmed FASTQs, BAM files, count matrices)
- Session state (selected references, parameter history, DE results, target predictions, enrichment tables)
- Provenance log

To save:

1. Go to **Project → Save/Load → Save Project**.
2. Enter a project name (e.g. `Medicago_drought_2026`).
3. Click **"Save"**.
4. A `.srnatlasproject` file is created in your chosen directory (default: `~/sRNAtlas_projects/`).

The file is a `.tar.gz` archive and can be shared with collaborators or archived for long-term storage.

### 7.5.2 Reloading a project

To restore a saved project:

1. Open sRNAtlas and go to **Project → Save/Load → Load Project**.
2. Click **"Browse"** and select your `.srnatlasproject` file.
3. Click **"Load"**.

sRNAtlas will:

- Extract all files to a temporary working directory.
- Restore the session state (sample list, references, parameters).
- Populate all modules with the saved results (count matrices, DE tables, plots, etc.).

You can now:

- Re-run specific stages with different parameters (e.g. adjust FDR threshold and re-run DE).
- Export additional plots or tables.
- Continue with downstream analysis (e.g. target prediction if it wasn't done initially).

**Important notes:**

- File paths in the saved project are relative to the project root, so the project is portable across systems.
- If you move raw FASTQ files to a different location, you will need to update the file paths in the **Project → Sample Management** interface before re-running upstream stages.

***

## 7.6 Provenance and reproducibility

### 7.6.1 Understanding the provenance log

The `provenance.yaml` file (saved in the output directory and embedded in the project archive) contains:

- **Run metadata**: date, user, sRNAtlas version
- **Tool versions**: Cutadapt, Bowtie, SAMtools, pyDESeq2, etc.
- **Parameters**: complete parameter sets for each stage
- **File checksums**: SHA-256 hashes of input FASTQs, reference FASTAs, and output count matrices
- **Command-line invocations**: exact commands used for external tools

Example snippet:

```yaml
run_id: sRNAtlas_20260213_211530
sRNAtlas_version: 1.0.0
steps:
  - stage: trimming
    tool: cutadapt
    version: 4.6
    parameters:
      adapter: TGGAATTCTCGGGTGCCAAGG
      min_length: 18
      max_length: 35
      quality_cutoff: 20
    input_files:
      - path: data/sample1.fastq.gz
        checksum: a3f2b1c4d5e6f7...
    output_files:
      - path: output/trimming/sample1_trimmed.fastq.gz
        checksum: 9e8d7c6b5a4f3e2...
  - stage: alignment
    tool: bowtie
    version: 1.3.1
    parameters:
      index: mtr_miRBase_v22
      mismatches: 1
      max_alignments: 10
    command: bowtie -v 1 -k 10 --best --strata ...
```

### 7.6.2 Using the provenance log

The provenance log enables:

1. **Exact reproduction**: Re-run the analysis with identical parameters by referencing the YAML file.
2. **Troubleshooting**: If results differ on a different system, compare tool versions and checksums.
3. **Publication**: Include the provenance file as supplementary material to document your methods.
4. **Compliance**: Satisfy journal or funding agency requirements for computational reproducibility.

***

## 7.7 Command-line interface for HPC environments

For users with access to a high-performance computing (HPC) cluster, sRNAtlas provides a command-line interface (CLI) that bypasses the Streamlit web interface and runs the batch pipeline as a headless job.

### 7.7.1 Installing the CLI

The CLI is included in the main sRNAtlas installation. After installing via pip or Conda, the `srnatlas` command should be available:

```bash
srnatlas --version
```

### 7.7.2 Preparing a configuration file

Instead of configuring parameters in the web UI, create a YAML configuration file:

```yaml
# config.yaml
project_name: Medicago_drought_2026
samples:
  - name: Control_rep1
    fastq: /path/to/Control_rep1.fastq.gz
    condition: Control
  - name: Control_rep2
    fastq: /path/to/Control_rep2.fastq.gz
    condition: Control
  - name: Drought_rep1
    fastq: /path/to/Drought_rep1.fastq.gz
    condition: Drought
  - name: Drought_rep2
    fastq: /path/to/Drought_rep2.fastq.gz
    condition: Drought

references:
  - type: miRBase
    organism: mtr
    version: 22
    path: /path/to/mtr_miRBase_v22

pipeline:
  stages:
    - qc_pre
    - trimming
    - qc_post
    - alignment
    - qc_align
    - counting
    - de_analysis
  
  parameters:
    trimming:
      adapter: TGGAATTCTCGGGTGCCAAGG
      min_length: 18
      max_length: 35
      quality_cutoff: 20
    alignment:
      mismatches: 1
      max_alignments: 10
      counting_mode: fractional
    de_analysis:
      design_factor: condition
      comparisons:
        - name: Drought_vs_Control
          test: Drought
          reference: Control
      fdr_threshold: 0.05
      lfc_threshold: 0.585

output_dir: /path/to/output/
```

### 7.7.3 Running the pipeline

Submit the job to your cluster's scheduler (example for SLURM):

```bash
#!/bin/bash
#SBATCH --job-name=sRNAtlas
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00

srnatlas batch --config config.yaml --threads 8
```

The CLI will:

- Parse the configuration file.
- Execute all selected stages sequentially.
- Write logs to `output_dir/logs/`.
- Save the provenance file to `output_dir/provenance.yaml`.

Progress and errors are logged to `stderr`, making it easy to monitor with `tail -f` or your cluster's log viewer.

### 7.7.4 Importing CLI results into the web interface

After the CLI batch run completes, you can load the results into the web interface:

1. Open sRNAtlas (web UI).
2. Go to **Project → Save/Load → Load from Directory**.
3. Browse to the `output_dir` from your CLI run.
4. Click **"Load"**.

The web interface will import:

- Count matrices
- DE results
- QC reports
- Provenance log

You can then visualize plots, run target prediction, or perform enrichment analysis interactively.

***

## 7.8 Best practices for batch processing

### 7.8.1 Test with a subset first

Before running the full pipeline on 100 samples, test on 2–3 representative samples interactively:

1. Run through the modules manually (Trimming → Alignment → Counting → DE).
2. Inspect QC metrics and alignment rates.
3. Adjust parameters if needed (e.g. increase minimum length if many reads are too short).
4. Once satisfied, apply the same parameters to the full batch.

### 7.8.2 Use checkpoints

If your dataset is large, enable checkpointing:

- In batch mode, set **"Checkpoint after each stage"** to `True`.
- This writes intermediate results to disk after each sample completes each stage.
- If the run is interrupted, you can resume from the last checkpoint without re-processing earlier stages.

### 7.8.3 Monitor disk space

Batch runs can generate large amounts of intermediate 

- Trimmed FASTQs: ~same size as raw FASTQs
- BAM files: ~0.5–2× the size of FASTQs, depending on alignment rate
- Temporary files: SAM files, unsorted BAMs

Estimate disk requirements:

- Raw  10 GB per sample (typical for 10M reads)
- Intermediate files: 20–30 GB per sample
- Final outputs (count matrices, DE results): < 1 GB total

For very large datasets (>100 samples), consider:

- Cleaning up intermediate files after each stage (`--cleanup` flag in CLI).
- Using a scratch directory on high-speed storage (e.g. SSD) for temporary files.

### 7.8.4 Archive completed projects

After publication or project completion:

1. Save the project as a `.srnatlasproject` file.
2. Compress the output directory (`tar -czf output.tar.gz output/`).
3. Store both the project file and the raw FASTQs in long-term archival storage (e.g. institutional repository, Zenodo, or SRA).
4. Delete intermediate files (trimmed FASTQs, BAMs) to free up space; these can be regenerated from the raw data + provenance file if needed.

***

## 7.9 Troubleshooting common batch issues

### Issue 1: "Reference index not found"

**Cause**: The reference FASTA was downloaded but not indexed.

**Fix**:

1. Go to **Databases** module.
2. Select the reference.
3. Click **"Build Index"**.
4. Wait for `bowtie-build` to complete.
5. Return to Batch and restart.

### Issue 2: Batch run hangs at alignment

**Cause**: Bowtie is single-threaded by default in the UI; for batch mode, it should use multiple cores.

**Fix**:

- In batch configuration, set **"Threads per sample"** to 4–8 (adjust based on your CPU).
- Or use the CLI with `--threads 8`.

### Issue 3: DE analysis fails with "Not enough replicates"

**Cause**: At least 2 replicates per condition are required for pyDESeq2.

**Fix**:

- Verify that your metadata has ≥2 samples per condition.
- If you only have 1 replicate per condition, DE testing is not statistically valid; consider using fold-change ranking instead (export normalized counts and compute FC manually).

### Issue 4: Out-of-memory error during counting

**Cause**: Large BAM files (e.g. from genome alignment) can exceed available RAM.

**Fix**:

- Increase memory allocation (CLI: `--mem 32G`; HPC: adjust SLURM `#SBATCH --mem`).
- Use sorted BAM files and enable streaming mode (reduces memory footprint).

***

## 7.10 Summary

This tutorial demonstrated how to:

- Configure and run a **batch pipeline** that automates trimming, alignment, counting, and DE analysis across multiple samples.
- **Save and reload** entire projects for reproducibility and collaboration.
- Extract and interpret the **provenance log** for methods documentation.
- Use the **command-line interface** for HPC cluster execution.

Combined with the earlier tutorials, you now have a complete workflow from raw FASTQ files through differential expression, target prediction, and functional enrichment—executed either interactively or in batch mode.

**Next steps:**

- **Tutorial 8** will cover advanced topics: hierarchical alignment, isomiR analysis, and novel miRNA discovery.
- **Tutorial 9** will demonstrate how to generate publication-ready figures and HTML reports.
