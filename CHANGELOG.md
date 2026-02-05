# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.3.0] - 2026-01-28

### Added
- **Adapter Trimming Module** - Full Cutadapt integration
  - Common adapter presets (Illumina TruSeq, NEBNext, QIAseq, Lexogen)
  - Quality and length filtering
  - Batch processing support
  - Detailed trimming statistics and visualizations
- **miRNA Target Prediction Module**
  - psRNATarget API integration for plants
  - Local seed matching algorithm for all organisms
  - Support for custom transcript databases
  - Results linked to enrichment analysis
- **Reference Database Builder Module**
  - Automatic miRBase download (mature and hairpin sequences)
  - RNAcentral integration
  - Custom database support
  - Automatic Bowtie index building
  - Database management interface
- **Batch Processing Module**
  - Full pipeline automation
  - Job queue management
  - Multi-comparison DE analysis
  - Progress tracking
- **Project Management Module**
  - Save/load analysis sessions (.srna format)
  - Project export as ZIP archive
  - Recent projects list
  - Auto-save option
  - Session data serialization

### Changed
- Reorganized navigation menu with new modules
- Updated home page with new features overview
- Removed R dependency (clusterProfiler) - now fully Python-based
- Improved g:Profiler integration with more organisms

---

## [1.2.0] - 2026-01-26

### Added
- **Hybrid Settings System** with centralized presets and per-module quick settings
  - 5 analysis presets: miRNA Strict, miRNA Relaxed, siRNA, piRNA, General Small RNA
  - Configuration export/import (JSON/YAML formats)
  - Auto-generated methods text for publications
  - Reset to defaults functionality
- **Professional color palettes** for visualizations (colorblind-friendly)

### Changed
- **Improved all visualizations** with better colors, sizing, and interactivity:
  - **PCA Plot**: Larger markers (size 12), proper color mapping, sample labels, equal axis scaling, better legend
  - **Volcano Plot**: Distinct up/down/non-sig colors, layered plotting, improved annotations
  - **MA Plot**: Direction-based coloring (up=red, down=blue), cleaner layout
  - **Heatmap**: Centered color scale, better sizing, clustering improvements
  - **Sample Correlation**: Correlation values displayed, custom color scale for high correlations
  - **Enrichment Dotplot**: Truncated long terms, Viridis-like colorscale, better sizing
- **Switched from Bowtie2 to Bowtie for small RNA alignment**
  - Bowtie is better optimized for short reads (<50 bp)
  - Provides exact mismatch control with `-v` mode
  - No gap alignment (appropriate for miRNA/siRNA)
  - Faster alignment for small RNA sequences
- Updated alignment module with Bowtie-specific settings:
  - Mismatch control (0-3 mismatches)
  - Best mode and strata mode support
  - Multi-mapper suppression option
- Updated all installation scripts (conda, linux, macos, server)
- Updated Dockerfile to use Bowtie instead of Bowtie2
- Updated documentation (INSTALL.md)

---

## [1.1.0] - 2026-01-26

### Added
- Reference Database Management system
  - Centralized reference selection UI component
  - Support for multiple species/organisms
  - Auto-detection of FASTA files in references folder
  - Reference metadata management (species, source, RNA types, etc.)
  - Integration with Project Setup and Alignment modules
- New "References" navigation menu item for managing reference databases
- `references.json` configuration file for reference metadata
- `reference_manager.py` utility class for backend reference handling
- `reference_module.py` with UI components for reference selection

### Changed
- Alignment module now uses centralized reference selector
- Project Setup page includes reference database selection

---

## [1.0.0] - 2026-01-25

### Added
- Initial release of sRNA-seq Analysis WebTool
- Differential Expression analysis using pyDESeq2
- GO/KEGG pathway enrichment via g:Profiler
- Interactive visualizations (Volcano, MA, PCA, Heatmap, Correlation)
- Demo dataset for quick testing
- Progress bars for long-running analyses
- Publication-ready figure export (HTML, PNG, SVG, PDF)
- Multi-platform installation scripts (macOS, Linux, Windows, Docker)
- Deployment guides (Streamlit Cloud, Debian Server)

### Features
- Quality Control module
- Alignment module (Bowtie)
- Read Counting module
- Differential Expression module
- Enrichment Analysis module
- Report Generation module
- Comprehensive Help & Documentation

---

## How to Update This File

When releasing a new version:

1. Update the `VERSION` file with new version number
2. Add a new section at the top of this file:

```markdown
## [X.Y.Z] - YYYY-MM-DD

### Added
- New features

### Changed
- Changes in existing functionality

### Deprecated
- Features that will be removed in future

### Removed
- Removed features

### Fixed
- Bug fixes

### Security
- Security improvements
```

## Version History

| Version | Date | Description |
|---------|------|-------------|
| 1.3.0 | 2026-01-28 | Adapter trimming, target prediction, database builder, batch processing, project management |
| 1.2.0 | 2026-01-26 | Switched to Bowtie for small RNA alignment |
| 1.1.0 | 2026-01-26 | Reference database management system |
| 1.0.0 | 2026-01-25 | Initial release |
