# Tutorial 1 – Getting Started: Installation, Project Setup, and Sample Management

Welcome to the sRNAtlas tutorial series! This tutorial will guide you through installing sRNAtlas, creating your first project, and organizing your samples and metadata. By the end of this tutorial, you will have a configured workspace ready for small RNA-seq analysis.

***

## 1.1 Installation

sRNAtlas can be installed via Conda, pip, or Docker. Choose the method that best fits your system and preferences.

### 1.1.1 Installation via Conda (recommended)

Conda manages all dependencies automatically and works across Linux, macOS, and Windows (via WSL2).

**Step 1: Install Miniconda or Anaconda**

If you don't have Conda installed:

- Download Miniconda from https://docs.conda.io/en/latest/miniconda.html
- Follow the installation instructions for your operating system
- Verify installation: `conda --version`

**Step 2: Create a dedicated environment**

```bash
conda create -n srnatlas python=3.10
conda activate srnatlas
```

**Step 3: Install sRNAtlas**

```bash
conda install -c bioconda -c conda-forge srnatlas
```

This installs sRNAtlas and all dependencies including Cutadapt, Bowtie, SAMtools, and Python packages.

**Step 4: Verify installation**

```bash
srnatlas --version
```

You should see output similar to:

```
sRNAtlas version 1.0.0
```

**Step 5: Test the installation**

```bash
srnatlas test
```

This runs a quick self-test to verify all components are working correctly.

### 1.1.2 Installation via pip

If you prefer pip or already have the external tools installed:

**Prerequisites:**

- Python 3.9 or higher
- Cutadapt (≥4.0): `pip install cutadapt`
- Bowtie v1 (≥1.3.0): Install from http://bowtie-bio.sourceforge.net/
- SAMtools (≥1.15): Install from http://www.htslib.org/

**Install sRNAtlas:**

```bash
pip install srnatlas
```

**Verify external tools are in PATH:**

```bash
cutadapt --version
bowtie --version
samtools --version
```

If any tool is not found, add its installation directory to your `PATH` environment variable.

### 1.1.3 Installation via Docker

Docker provides a pre-configured container with all dependencies.

**Step 1: Install Docker**

- Download Docker Desktop from https://www.docker.com/products/docker-desktop
- Follow installation instructions for your OS
- Verify: `docker --version`

**Step 2: Pull the sRNAtlas image**

```bash
docker pull osmanradwan/srnatlas:latest
```

**Step 3: Run the container**

```bash
docker run -p 8501:8501 -v /path/to/your//data osmanradwan/srnatlas:latest
```

Replace `/path/to/your/data` with the directory containing your FASTQ files.

**Step 4: Access the interface**

Open your web browser and navigate to:

```
http://localhost:8501
```

The sRNAtlas interface should appear.

***

## 1.2 Launching sRNAtlas

### 1.2.1 Starting the web interface

**For Conda or pip installations:**

```bash
srnatlas run
```

This starts the Streamlit web server. You should see output like:

```
You can now view sRNAtlas in your browser.

  Local URL: http://localhost:8501
  Network URL: http://192.168.1.100:8501
```

Open the Local URL in your web browser (Chrome, Firefox, or Safari recommended).

**For Docker installations:**

The container starts automatically with the web interface accessible at `http://localhost:8501`.

### 1.2.2 Interface overview

When you first open sRNAtlas, you'll see:

- **Left sidebar**: Navigation menu with all 14 modules
- **Main panel**: Home page with quick start guide
- **Top bar**: sRNAtlas logo and version

The modules are organized into four phases (matching the manuscript Figure 1):

**Phase 1 — Preprocessing:**
- Project Setup
- Quality Control
- Trimming
- Databases

**Phase 2 — Mapping and Quantification:**
- Alignment
- Post-Align QC
- Counting

**Phase 3 — Statistical Analysis:**
- DE Analysis
- Novel miRNA
- isomiR Analysis

**Phase 4 — Biological Interpretation:**
- Target Prediction
- GO/Pathway Enrichment
- Batch Processing
- Reports

***

## 1.3 Creating your first project

### 1.3.1 Project concept

In sRNAtlas, a **project** is a container for:

- Sample metadata (names, conditions, file paths)
- Analysis parameters (adapter sequences, alignment settings, DE design)
- Results (count matrices, DE tables, plots)
- Session state (so you can pick up where you left off)

All projects are saved as `.srnatlasproject` files that can be shared with collaborators or archived.

### 1.3.2 Starting a new project

**Step 1: Navigate to Project Setup**

Click **"Project Setup"** in the left sidebar.

**Step 2: Create a new project**

- Click the **"New Project"** button
- Enter a project name (e.g., `Medicago_drought_2026`)
- Enter a brief description (optional, e.g., "Small RNA profiling of *M. truncatula* under drought stress")
- Select a working directory where all project files will be stored (default: `~/sRNAtlas_projects/`)
- Click **"Create Project"**

The interface now shows the project dashboard with three tabs:

1. **Samples** – Manage sample metadata and file associations
2. **Metadata** – Define experimental design and factors
3. **Settings** – Configure global project settings

***

## 1.4 Adding samples to your project

### 1.4.1 Sample information requirements

For each sample, you need:

- **Sample name**: Unique identifier (e.g., `Control_rep1`, `Drought_rep2`)
- **FASTQ file path**: Location of the raw sequencing file
- **Condition/treatment**: Experimental group (e.g., `Control`, `Drought`)
- **Optional metadata**: Batch, sequencing run, tissue type, time point, etc.

### 1.4.2 Adding samples manually

**Step 1: Go to the Samples tab**

In the Project Setup module, click the **"Samples"** tab.

**Step 2: Add a sample**

Click **"Add Sample"** and fill in the form:

- **Sample name**: `Control_rep1`
- **FASTQ file**: Click "Browse" and navigate to your file, or paste the full path (e.g., `/data/raw_reads/Control_rep1.fastq.gz`)
- **Condition**: `Control`
- **Replicate**: `1` (optional, for numbering)
- **Notes**: Any sample-specific notes (optional)

Click **"Save"**.

Repeat for all samples. For example, a minimal experiment might have:

| Sample name | FASTQ file | Condition | Replicate |
|---|---|---|---|
| Control_rep1 | /data/Control_rep1.fastq.gz | Control | 1 |
| Control_rep2 | /data/Control_rep2.fastq.gz | Control | 2 |
| Control_rep3 | /data/Control_rep3.fastq.gz | Control | 3 |
| Drought_rep1 | /data/Drought_rep1.fastq.gz | Drought | 1 |
| Drought_rep2 | /data/Drought_rep2.fastq.gz | Drought | 2 |
| Drought_rep3 | /data/Drought_rep3.fastq.gz | Drought | 3 |

### 1.4.3 Importing samples from CSV

For projects with many samples (10+), use the CSV import function:

**Step 1: Prepare a CSV file**

Create a CSV file with the following columns (exact header names required):

```csv
sample_name,fastq_path,condition,replicate,batch,notes
Control_rep1,/data/Control_rep1.fastq.gz,Control,1,Batch1,
Control_rep2,/data/Control_rep2.fastq.gz,Control,2,Batch1,
Control_rep3,/data/Control_rep3.fastq.gz,Control,3,Batch2,
Drought_rep1,/data/Drought_rep1.fastq.gz,Drought,1,Batch1,
Drought_rep2,/data/Drought_rep2.fastq.gz,Drought,2,Batch1,
Drought_rep3,/data/Drought_rep3.fastq.gz,Drought,3,Batch2,High_quality_RNA
```

**Required columns:**
- `sample_name`
- `fastq_path`
- `condition`

**Optional columns:**
- `replicate`
- `batch` (important if you have sequencing batches—used for batch correction)
- `tissue` (if comparing multiple tissues)
- `time_point` (for time-course experiments)
- `notes`

**Step 2: Import the CSV**

In the Samples tab, click **"Import from CSV"**, upload your file, and click **"Load"**. sRNAtlas will validate the file and display a preview. If everything looks correct, click **"Confirm Import"**.

### 1.4.4 Verifying file paths

After adding samples, sRNAtlas automatically checks that all FASTQ files exist:

- ✅ Green checkmark: File found
- ❌ Red X: File not found (check the path)

If files are missing:

1. Click the ❌ icon next to the sample
2. Update the file path
3. Click **"Update"**

***

## 1.5 Defining experimental design

### 1.5.1 Metadata tab

The **Metadata** tab lets you define factors (variables) for downstream analysis.

**Step 1: Go to the Metadata tab**

Click **"Metadata"** in the Project Setup module.

**Step 2: Define factors**

You'll see a table with all samples and their associated metadata. Factors are the columns you imported (e.g., `condition`, `replicate`, `batch`).

**For differential expression analysis (Tutorial 5), you need:**

- **At least one treatment factor** (e.g., `condition` with levels `Control`, `Drought`)
- **At least 2 replicates per level** (statistical requirement for DESeq2)

**Step 3: Set the primary factor**

Select the **main factor** for your experiment from the dropdown:

- Primary factor: `condition`

This is the factor you'll use for differential expression comparisons (e.g., Drought vs. Control).

**Step 4: Add additional factors (optional)**

If you have batch effects or multiple variables, add them:

- Click **"Add Factor"**
- Enter factor name (e.g., `sequencing_batch`)
- For each sample, select or enter the level (e.g., `Batch1`, `Batch2`)

**Example experimental design with batch:**

| Sample | condition | replicate | sequencing_batch |
|---|---|---|---|
| Control_rep1 | Control | 1 | Batch1 |
| Control_rep2 | Control | 2 | Batch1 |
| Control_rep3 | Control | 3 | Batch2 |
| Drought_rep1 | Drought | 1 | Batch1 |
| Drought_rep2 | Drought | 2 | Batch1 |
| Drought_rep3 | Drought | 3 | Batch2 |

In this design:

- Primary factor: `condition` (2 levels: Control, Drought)
- Nuisance factor: `sequencing_batch` (2 levels: Batch1, Batch2)

In the DE analysis (Tutorial 5), you can include `sequencing_batch` as a covariate to account for batch effects.

### 1.5.2 Validating your design

Click **"Validate Design"**. sRNAtlas checks:

- ✅ All samples have values for all factors
- ✅ At least 2 replicates per condition (required for DE)
- ⚠️ Confounding detected (e.g., if all Control samples are in Batch1 and all Drought samples are in Batch2—this makes batch and condition inseparable)

If validation passes, you're ready for analysis. If warnings appear, review your design and adjust if needed.

***

## 1.6 Project settings

### 1.6.1 Global settings

Go to the **Settings** tab to configure project-wide parameters:

**Organism:**
- Select your organism from the dropdown (e.g., *Medicago truncatula*, *Arabidopsis thaliana*, *Homo sapiens*)
- This affects:
  - QC interpretation (plant vs. animal RNA size classes)
  - Default references in the Databases module
  - Organism codes for psRNATarget and g:Profiler

**Working directory structure:**

sRNAtlas creates subdirectories automatically:

```
your_project/
├── raw_reads/          # (you put your FASTQ files here)
├── trimmed/            # Cutadapt output
├── alignment/          # BAM files
├── counts/             # Count matrices
├── de_results/         # Differential expression tables
├── figures/            # Exported plots
├── reports/            # HTML reports
└── project_data.json   # Project metadata
```

**File naming convention:**

- Default: `[sample_name]_[stage].[ext]` (e.g., `Control_rep1_trimmed.fastq.gz`)
- Custom: You can define your own template if needed

**Number of threads:**

- Set the number of CPU cores to use for parallel tasks (Bowtie, Cutadapt)
- Default: 4 (adjust based on your system: 1–16)

**Auto-save interval:**

- sRNAtlas auto-saves your session every *N* minutes (default: 5)
- This prevents data loss if the browser closes unexpectedly

***

## 1.7 Saving and loading projects

### 1.7.1 Saving your project

At any time, click the **"Save Project"** button (top-right corner, or in the Project Setup module).

This creates a `.srnatlasproject` file containing:

- All sample metadata
- Analysis parameters
- Results (count matrices, DE tables, etc.)
- Session state (which modules have been run, selected references, etc.)

**Save location:**

By default, the project file is saved in the working directory:

```
~/sRNAtlas_projects/Medicago_drought_2026/Medicago_drought_2026.srnatlasproject
```

You can also click **"Save Project As..."** to choose a different name or location.

### 1.7.2 Loading an existing project

To resume work on a saved project:

1. Open sRNAtlas
2. Click **"Load Project"** (home page or Project Setup module)
3. Browse to your `.srnatlasproject` file
4. Click **"Load"**

All modules will restore their previous state. You can:

- Continue from where you left off (e.g., if you stopped after trimming, go to Alignment)
- Re-run stages with different parameters (e.g., change alignment settings and re-align)
- Add more samples and re-run analyses

### 1.7.3 Sharing projects

To share your analysis with a collaborator:

1. Save the project
2. Copy the `.srnatlasproject` file **and** the working directory (contains all FASTQ files, BAM files, etc.)
3. Send both to your collaborator (e.g., via shared drive or USB)
4. Collaborator opens sRNAtlas and loads the project file

**File paths are relative** to the project root, so the project should work on any system as long as the directory structure is preserved.

***

## 1.8 Example: Setting up a real project

Let's walk through a complete example.

**Scenario:**

You have small RNA-seq data from *Medicago truncatula* leaves:

- 4 samples: 2 control, 2 drought-stressed
- Goal: Identify drought-responsive miRNAs

**Files:**

```
/home/user/data/
├── Control_rep1.fastq.gz
├── Control_rep2.fastq.gz
├── Drought_rep1.fastq.gz
└── Drought_rep2.fastq.gz
```

**Setup steps:**

1. **Launch sRNAtlas:**

   ```bash
   conda activate srnatlas
   srnatlas run
   ```

   Open http://localhost:8501

2. **Create project:**

   - Click "New Project"
   - Name: `Medicago_drought_stress`
   - Description: "Drought-responsive miRNAs in M. truncatula leaves"
   - Working directory: `/home/user/sRNAtlas_projects/`
   - Click "Create Project"

3. **Add samples (CSV method):**

   Create `samples.csv`:

   ```csv
   sample_name,fastq_path,condition,replicate
   Control_rep1,/home/user/data/Control_rep1.fastq.gz,Control,1
   Control_rep2,/home/user/data/Control_rep2.fastq.gz,Control,2
   Drought_rep1,/home/user/data/Drought_rep1.fastq.gz,Drought,1
   Drought_rep2,/home/user/data/Drought_rep2.fastq.gz,Drought,2
   ```

   - In Samples tab, click "Import from CSV"
   - Upload `samples.csv`
   - Confirm import

4. **Verify files:**

   - Check that all 4 samples show ✅ (files found)

5. **Define experimental design:**

   - Go to Metadata tab
   - Primary factor: `condition`
   - Click "Validate Design"
   - Should pass with message: "Design valid: 2 replicates per condition"

6. **Configure settings:**

   - Go to Settings tab
   - Organism: *Medicago truncatula*
   - Threads: 4 (or adjust for your system)
   - Click "Save Settings"

7. **Save the project:**

   - Click "Save Project"
   - Confirmation: "Project saved successfully"

**Result:**

Your project is now set up and ready for analysis. You can proceed to:

- **Tutorial 2**: Download and index references (miRBase *Medicago* miRNAs)
- **Tutorial 3**: Run quality control and adapter trimming
- **Tutorial 4**: Align reads to the reference

***

## 1.9 Troubleshooting common setup issues

**Issue 1: FASTQ files not found after import**

- **Cause**: Incorrect file paths, or files moved after import
- **Fix**: 
  - In Samples tab, click the ❌ icon next to the affected sample
  - Update the `fastq_path` field with the correct absolute path
  - Click "Update"

**Issue 2: "Not enough replicates" warning in design validation**

- **Cause**: Fewer than 2 samples per condition
- **Fix**: Differential expression requires at least 2 replicates per condition. Either:
  - Add more samples, or
  - Skip DE analysis and proceed with descriptive statistics only

**Issue 3: CSV import fails with "Invalid format"**

- **Cause**: Missing required columns or incorrect column names
- **Fix**: Ensure your CSV has these **exact** headers: `sample_name`, `fastq_path`, `condition`. Other columns are optional.

**Issue 4: sRNAtlas web interface won't load**

- **Cause**: Port 8501 already in use, or Streamlit crashed
- **Fix**:
  - Check if another instance is running: `ps aux | grep streamlit`
  - Kill any existing processes: `pkill -f streamlit`
  - Restart: `srnatlas run`
  - If still fails, try a different port: `streamlit run --server.port 8502`

***

## 1.10 Summary

In this tutorial, you learned how to:

- **Install sRNAtlas** via Conda, pip, or Docker
- **Launch the web interface** and navigate the module structure
- **Create a new project** with a descriptive name and working directory
- **Add samples** manually or via CSV import
- **Define experimental design** with factors and replicates
- **Validate your design** to ensure it meets statistical requirements
- **Save and load projects** for reproducibility and sharing

You now have a configured project ready for the next steps in the analysis pipeline.

**Next tutorial:**

- **Tutorial 2** will guide you through setting up reference databases (miRBase and RNAcentral) for alignment and annotation.
