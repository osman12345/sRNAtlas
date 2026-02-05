"""
Settings Module for sRNAtlas
Centralized configuration with presets and export/import functionality

Hybrid Approach:
- Quick settings: Available within each module for frequently-changed parameters
- Central Settings: Default values, presets, global options, export/import
"""
import streamlit as st
import json
import yaml
from pathlib import Path
from datetime import datetime
from typing import Dict, Any

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from config import config

# ============================================================================
# PRESET DEFINITIONS
# ============================================================================

PRESETS = {
    "miRNA_strict": {
        "name": "miRNA Strict",
        "description": "Strict settings for mature miRNA analysis. Zero mismatches, suppress multi-mappers with >10 hits.",
        "icon": "🎯",
        "settings": {
            "alignment": {
                "mismatches": 0,
                "max_alignments": 5,
                "suppress_multi": 10,
                "best_mode": True,
                "strata": True,
                "min_length": 20,
                "max_length": 24
            },
            "qc": {
                "min_quality": 25,
                "min_length": 18,
                "max_length": 26
            },
            "counting": {
                "min_mapq": 0,
                "min_counts_per_sample": 5,
                "min_samples_detected": 2
            },
            "deseq": {
                "min_cpm": 1.0,
                "padj_threshold": 0.05,
                "lfc_threshold": 1.0
            }
        }
    },
    "miRNA_relaxed": {
        "name": "miRNA Relaxed",
        "description": "Relaxed settings for miRNA discovery. Allows 1 mismatch for isomiR detection.",
        "icon": "🔬",
        "settings": {
            "alignment": {
                "mismatches": 1,
                "max_alignments": 10,
                "suppress_multi": 0,
                "best_mode": True,
                "strata": True,
                "min_length": 18,
                "max_length": 26
            },
            "qc": {
                "min_quality": 20,
                "min_length": 17,
                "max_length": 28
            },
            "counting": {
                "min_mapq": 0,
                "min_counts_per_sample": 1,
                "min_samples_detected": 2
            },
            "deseq": {
                "min_cpm": 0.5,
                "padj_threshold": 0.05,
                "lfc_threshold": 0.585
            }
        }
    },
    "siRNA": {
        "name": "siRNA Analysis",
        "description": "Settings optimized for siRNA (21-24nt). Allows multi-mapping for repetitive targets.",
        "icon": "🧬",
        "settings": {
            "alignment": {
                "mismatches": 1,
                "max_alignments": 20,
                "suppress_multi": 0,
                "best_mode": True,
                "strata": False,
                "min_length": 20,
                "max_length": 25
            },
            "qc": {
                "min_quality": 20,
                "min_length": 20,
                "max_length": 25
            },
            "counting": {
                "min_mapq": 0,
                "min_counts_per_sample": 1,
                "min_samples_detected": 2
            },
            "deseq": {
                "min_cpm": 0.5,
                "padj_threshold": 0.05,
                "lfc_threshold": 0.585
            }
        }
    },
    "piRNA": {
        "name": "piRNA Analysis",
        "description": "Settings for piRNA analysis (24-32nt). Allows more mismatches for ping-pong signature detection.",
        "icon": "🔄",
        "settings": {
            "alignment": {
                "mismatches": 2,
                "max_alignments": 50,
                "suppress_multi": 0,
                "best_mode": True,
                "strata": False,
                "min_length": 24,
                "max_length": 32
            },
            "qc": {
                "min_quality": 20,
                "min_length": 23,
                "max_length": 35
            },
            "counting": {
                "min_mapq": 0,
                "min_counts_per_sample": 1,
                "min_samples_detected": 2
            },
            "deseq": {
                "min_cpm": 0.5,
                "padj_threshold": 0.05,
                "lfc_threshold": 0.585
            }
        }
    },
    "general_sRNA": {
        "name": "General Small RNA",
        "description": "Balanced settings for mixed small RNA analysis. Good starting point for exploratory analysis.",
        "icon": "⚖️",
        "settings": {
            "alignment": {
                "mismatches": 1,
                "max_alignments": 10,
                "suppress_multi": 0,
                "best_mode": True,
                "strata": True,
                "min_length": 18,
                "max_length": 35
            },
            "qc": {
                "min_quality": 20,
                "min_length": 18,
                "max_length": 40
            },
            "counting": {
                "min_mapq": 0,
                "min_counts_per_sample": 1,
                "min_samples_detected": 2
            },
            "deseq": {
                "min_cpm": 0.5,
                "padj_threshold": 0.05,
                "lfc_threshold": 0.585
            }
        }
    }
}


def render_settings_page():
    """Render the settings page with hybrid approach"""
    st.header("⚙️ Settings")

    st.info("""
    **Hybrid Settings Approach**

    - **Quick settings** are available within each module (Alignment, DE Analysis, etc.) for parameters you change frequently
    - **This page** manages default values, presets, and configuration export/import
    """)

    # Main tabs
    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
        "🎛️ Presets",
        "🔗 Alignment",
        "📊 QC & Counting",
        "🔬 DE Analysis",
        "🧬 Enrichment",
        "💾 Import/Export"
    ])

    with tab1:
        render_presets_tab()

    with tab2:
        render_alignment_settings()

    with tab3:
        render_qc_counting_settings()

    with tab4:
        render_de_settings()

    with tab5:
        render_enrichment_settings()

    with tab6:
        render_import_export_tab()


def render_presets_tab():
    """Render preset selection and management"""
    st.subheader("🎛️ Analysis Presets")

    st.markdown("""
    Select a preset to quickly configure all settings for your specific analysis type.
    Presets set optimal defaults for alignment, QC, counting, and DE analysis.
    """)

    # Current preset indicator
    current_preset = st.session_state.get('current_preset', None)
    if current_preset:
        st.success(f"Current preset: **{PRESETS[current_preset]['name']}**")

    st.divider()

    # Display presets as cards
    cols = st.columns(2)

    for idx, (preset_id, preset) in enumerate(PRESETS.items()):
        col = cols[idx % 2]

        with col:
            with st.container():
                st.markdown(f"### {preset['icon']} {preset['name']}")
                st.markdown(f"*{preset['description']}*")

                # Show key settings
                with st.expander("View settings"):
                    settings = preset['settings']

                    st.markdown("**Alignment:**")
                    st.markdown(f"- Mismatches: {settings['alignment']['mismatches']}")
                    st.markdown(f"- Max alignments: {settings['alignment']['max_alignments']}")
                    st.markdown(f"- Size range: {settings['alignment']['min_length']}-{settings['alignment']['max_length']} nt")

                    st.markdown("**DE Analysis:**")
                    st.markdown(f"- FDR: {settings['deseq']['padj_threshold']}")
                    st.markdown(f"- Log2FC: {settings['deseq']['lfc_threshold']}")

                if st.button(f"Apply {preset['name']}", key=f"apply_{preset_id}", type="primary"):
                    apply_preset(preset_id)
                    st.success(f"✅ Applied preset: {preset['name']}")
                    st.rerun()

                st.markdown("---")

    # Reset to defaults
    st.divider()
    if st.button("🔄 Reset All to Defaults", type="secondary"):
        reset_to_defaults()
        st.session_state['current_preset'] = None
        st.success("All settings reset to defaults!")
        st.rerun()


def apply_preset(preset_id: str):
    """Apply a preset to all settings"""
    if preset_id not in PRESETS:
        return

    preset = PRESETS[preset_id]
    settings = preset['settings']

    # Apply alignment settings
    if 'alignment' in settings:
        for key, value in settings['alignment'].items():
            if hasattr(config.alignment, key):
                setattr(config.alignment, key, value)

    # Apply QC settings
    if 'qc' in settings:
        for key, value in settings['qc'].items():
            if hasattr(config.qc, key):
                setattr(config.qc, key, value)

    # Apply counting settings
    if 'counting' in settings:
        for key, value in settings['counting'].items():
            if hasattr(config.counting, key):
                setattr(config.counting, key, value)

    # Apply DE settings
    if 'deseq' in settings:
        for key, value in settings['deseq'].items():
            if hasattr(config.deseq, key):
                setattr(config.deseq, key, value)

    # Store current preset
    st.session_state['current_preset'] = preset_id


def reset_to_defaults():
    """Reset all settings to their default values"""
    # Reset alignment
    config.alignment.mismatches = 1
    config.alignment.max_alignments = 10
    config.alignment.suppress_multi = 0
    config.alignment.best_mode = True
    config.alignment.strata = True
    config.alignment.min_length = 18
    config.alignment.max_length = 35
    config.alignment.threads = 8

    # Reset QC
    config.qc.min_quality = 20
    config.qc.min_length = 18
    config.qc.max_length = 50
    config.qc.trim_quality = 20

    # Reset counting
    config.counting.min_mapq = 0
    config.counting.min_counts_per_sample = 1
    config.counting.min_samples_detected = 2

    # Reset DE
    config.deseq.min_cpm = 0.5
    config.deseq.padj_threshold = 0.05
    config.deseq.lfc_threshold = 0.585


def render_alignment_settings():
    """Render Bowtie alignment default settings"""
    st.subheader("🔗 Alignment Defaults (Bowtie)")

    st.markdown("""
    These are the **default values** used when you start a new alignment.
    You can override them in the Alignment module for individual runs.
    """)

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Core Settings**")

        threads = st.number_input(
            "Default Threads",
            min_value=1,
            max_value=64,
            value=config.alignment.threads,
            help="Number of parallel threads for alignment"
        )

        mismatches = st.selectbox(
            "Default Mismatches (-v)",
            options=[0, 1, 2, 3],
            index=config.alignment.mismatches,
            help="Number of mismatches allowed. 0-1 recommended for miRNA."
        )

        max_alignments = st.number_input(
            "Default Max Alignments (-k)",
            min_value=1,
            max_value=100,
            value=config.alignment.max_alignments,
            help="Report up to k valid alignments per read"
        )

        suppress_multi = st.number_input(
            "Suppress Multi-mappers (-m)",
            min_value=0,
            max_value=100,
            value=config.alignment.suppress_multi,
            help="Suppress reads with more than m alignments (0=disabled)"
        )

    with col2:
        st.markdown("**Mode Settings**")

        best_mode = st.checkbox(
            "Best Mode (--best)",
            value=config.alignment.best_mode,
            help="Report alignments in best-to-worst order"
        )

        strata = st.checkbox(
            "Strata Mode (--strata)",
            value=config.alignment.strata,
            help="Only report alignments in the best stratum"
        )

        st.markdown("**Size Selection**")

        min_length = st.number_input(
            "Min Read Length",
            min_value=10,
            max_value=50,
            value=config.alignment.min_length
        )

        max_length = st.number_input(
            "Max Read Length",
            min_value=20,
            max_value=100,
            value=config.alignment.max_length
        )

    if st.button("💾 Save Alignment Defaults", key="save_alignment"):
        config.alignment.threads = threads
        config.alignment.mismatches = mismatches
        config.alignment.max_alignments = max_alignments
        config.alignment.suppress_multi = suppress_multi
        config.alignment.best_mode = best_mode
        config.alignment.strata = strata
        config.alignment.min_length = min_length
        config.alignment.max_length = max_length
        st.session_state['current_preset'] = None  # Clear preset since user customized
        st.success("✅ Alignment defaults saved!")


def render_qc_counting_settings():
    """Render QC and Counting settings"""
    st.subheader("📊 QC & Counting Defaults")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Quality Control**")

        min_quality = st.number_input(
            "Min Quality Score",
            min_value=0,
            max_value=40,
            value=config.qc.min_quality,
            help="Minimum Phred quality score"
        )

        trim_quality = st.number_input(
            "Trim Quality Threshold",
            min_value=0,
            max_value=40,
            value=config.qc.trim_quality
        )

        adapter_sequence = st.text_input(
            "Adapter Sequence",
            value=config.qc.adapter_sequence,
            help="Illumina small RNA adapter"
        )

    with col2:
        st.markdown("**Read Counting**")

        min_mapq = st.number_input(
            "Min Mapping Quality",
            min_value=0,
            max_value=60,
            value=config.counting.min_mapq,
            help="Set to 0 to include multi-mappers (recommended for small RNA)"
        )

        count_mode = st.selectbox(
            "Count Mode",
            options=['union', 'intersection-strict', 'intersection-nonempty'],
            index=0
        )

        strand_specific = st.checkbox(
            "Strand-Specific",
            value=config.counting.strand_specific
        )

        st.markdown("**Filtering**")

        min_counts = st.number_input(
            "Min Counts per Sample",
            min_value=0,
            max_value=100,
            value=config.counting.min_counts_per_sample
        )

        min_samples = st.number_input(
            "Min Samples Detected",
            min_value=1,
            max_value=20,
            value=config.counting.min_samples_detected
        )

    if st.button("💾 Save QC & Counting Defaults", key="save_qc"):
        config.qc.min_quality = min_quality
        config.qc.trim_quality = trim_quality
        config.qc.adapter_sequence = adapter_sequence
        config.counting.min_mapq = min_mapq
        config.counting.count_mode = count_mode
        config.counting.strand_specific = strand_specific
        config.counting.min_counts_per_sample = min_counts
        config.counting.min_samples_detected = min_samples
        st.session_state['current_preset'] = None
        st.success("✅ QC & Counting defaults saved!")


def render_de_settings():
    """Render DE analysis settings"""
    st.subheader("🔬 Differential Expression Defaults")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Filtering**")

        min_cpm = st.number_input(
            "Minimum CPM",
            min_value=0.0,
            max_value=10.0,
            value=float(config.deseq.min_cpm),
            step=0.1,
            help="Minimum counts per million for filtering"
        )

        min_samples = st.number_input(
            "Min Samples",
            min_value=1,
            max_value=20,
            value=config.deseq.min_samples,
            help="Gene must have min CPM in at least this many samples"
        )

        st.markdown("**Significance Thresholds**")

        padj_threshold = st.number_input(
            "FDR Threshold",
            min_value=0.01,
            max_value=0.20,
            value=config.deseq.padj_threshold,
            step=0.01,
            help="Adjusted p-value threshold"
        )

        lfc_threshold = st.number_input(
            "Log2 Fold Change Threshold",
            min_value=0.0,
            max_value=3.0,
            value=config.deseq.lfc_threshold,
            step=0.1,
            help="Minimum absolute log2 fold change"
        )

    with col2:
        st.markdown("**Batch Correction (SVA)**")

        use_sva = st.checkbox(
            "Use SVA",
            value=config.deseq.use_sva,
            help="Surrogate Variable Analysis for batch correction"
        )

        if use_sva:
            max_sv = st.number_input(
                "Max Surrogate Variables",
                min_value=1,
                max_value=10,
                value=config.deseq.max_surrogate_variables
            )

            sv_threshold = st.number_input(
                "SV Correlation Warning Threshold",
                min_value=0.1,
                max_value=0.9,
                value=config.deseq.sv_correlation_threshold,
                step=0.1
            )

        st.markdown("**Common Thresholds Reference**")
        st.markdown("""
        | Analysis | FDR | Log2FC |
        |----------|-----|--------|
        | Strict | 0.01 | 1.0 |
        | Standard | 0.05 | 0.585 |
        | Relaxed | 0.10 | 0.0 |
        """)

    if st.button("💾 Save DE Defaults", key="save_de"):
        config.deseq.min_cpm = min_cpm
        config.deseq.min_samples = min_samples
        config.deseq.padj_threshold = padj_threshold
        config.deseq.lfc_threshold = lfc_threshold
        config.deseq.use_sva = use_sva
        if use_sva:
            config.deseq.max_surrogate_variables = max_sv
            config.deseq.sv_correlation_threshold = sv_threshold
        st.session_state['current_preset'] = None
        st.success("✅ DE defaults saved!")


def render_enrichment_settings():
    """Render enrichment settings"""
    st.subheader("🧬 GO/Pathway Enrichment Defaults")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**GO Ontologies**")
        go_bp = st.checkbox("Biological Process (BP)", value='BP' in config.enrichment.go_ontologies)
        go_mf = st.checkbox("Molecular Function (MF)", value='MF' in config.enrichment.go_ontologies)
        go_cc = st.checkbox("Cellular Component (CC)", value='CC' in config.enrichment.go_ontologies)

        enrichment_padj = st.number_input(
            "Enrichment FDR Threshold",
            min_value=0.01,
            max_value=0.20,
            value=config.enrichment.enrichment_padj,
            step=0.01
        )

    with col2:
        st.markdown("**Pathways**")

        use_kegg = st.checkbox(
            "Include KEGG Pathways",
            value=config.enrichment.use_kegg
        )

        st.markdown("**Gene Set Size Filters**")

        min_gs_size = st.number_input(
            "Min Gene Set Size",
            min_value=1,
            max_value=50,
            value=config.enrichment.min_gene_set_size
        )

        max_gs_size = st.number_input(
            "Max Gene Set Size",
            min_value=100,
            max_value=2000,
            value=config.enrichment.max_gene_set_size
        )

        top_terms = st.number_input(
            "Top Terms to Display",
            min_value=5,
            max_value=100,
            value=config.enrichment.top_terms_to_show
        )

    if st.button("💾 Save Enrichment Defaults", key="save_enrichment"):
        ontologies = []
        if go_bp:
            ontologies.append('BP')
        if go_mf:
            ontologies.append('MF')
        if go_cc:
            ontologies.append('CC')

        config.enrichment.go_ontologies = ontologies
        config.enrichment.use_kegg = use_kegg
        config.enrichment.enrichment_padj = enrichment_padj
        config.enrichment.min_gene_set_size = min_gs_size
        config.enrichment.max_gene_set_size = max_gs_size
        config.enrichment.top_terms_to_show = top_terms
        st.session_state['current_preset'] = None
        st.success("✅ Enrichment defaults saved!")


def render_import_export_tab():
    """Render import/export functionality"""
    st.subheader("💾 Configuration Management")

    st.markdown("""
    Export your settings to share with collaborators or for reproducibility.
    Import settings from a previous analysis or from a colleague.
    """)

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("### 📤 Export Configuration")

        export_format = st.radio(
            "Export Format",
            options=["JSON", "YAML"],
            horizontal=True
        )

        # Get current settings
        current_settings = get_all_settings()

        # Add metadata
        export_data = {
            "metadata": {
                "tool": "sRNAtlas",
                "version": "1.2.0",
                "exported_at": datetime.now().isoformat(),
                "preset": st.session_state.get('current_preset', None)
            },
            "settings": current_settings
        }

        if export_format == "JSON":
            export_str = json.dumps(export_data, indent=2)
            filename = "srna_config.json"
            mime = "application/json"
        else:
            export_str = yaml.dump(export_data, default_flow_style=False, sort_keys=False)
            filename = "srna_config.yaml"
            mime = "text/yaml"

        st.download_button(
            f"📥 Download {export_format}",
            export_str,
            filename,
            mime,
            type="primary"
        )

        with st.expander("Preview export"):
            st.code(export_str, language=export_format.lower())

    with col2:
        st.markdown("### 📥 Import Configuration")

        uploaded_file = st.file_uploader(
            "Upload configuration file",
            type=['json', 'yaml', 'yml']
        )

        if uploaded_file:
            try:
                content = uploaded_file.read().decode('utf-8')

                if uploaded_file.name.endswith('.json'):
                    imported = json.loads(content)
                else:
                    imported = yaml.safe_load(content)

                st.success(f"✅ Loaded: {uploaded_file.name}")

                # Show metadata if present
                if 'metadata' in imported:
                    meta = imported['metadata']
                    st.info(f"**From:** {meta.get('tool', 'Unknown')} v{meta.get('version', '?')}")
                    if meta.get('preset'):
                        st.info(f"**Preset:** {meta['preset']}")

                with st.expander("Preview imported settings"):
                    st.json(imported.get('settings', imported))

                if st.button("✅ Apply Imported Settings", type="primary"):
                    settings = imported.get('settings', imported)
                    apply_imported_settings(settings)
                    st.success("Settings applied successfully!")
                    st.rerun()

            except Exception as e:
                st.error(f"Error loading file: {e}")

    # Quick export for publications
    st.divider()
    st.markdown("### 📋 Copy Settings for Publication")
    st.markdown("Copy this text to include in your methods section:")

    methods_text = generate_methods_text()
    st.code(methods_text, language=None)

    st.button("📋 Copy to Clipboard", disabled=True, help="Use Ctrl+C to copy from the box above")


def get_all_settings() -> Dict[str, Any]:
    """Get all current settings as a dictionary"""
    return {
        "alignment": {
            "mismatches": config.alignment.mismatches,
            "max_alignments": config.alignment.max_alignments,
            "suppress_multi": config.alignment.suppress_multi,
            "best_mode": config.alignment.best_mode,
            "strata": config.alignment.strata,
            "threads": config.alignment.threads,
            "min_length": config.alignment.min_length,
            "max_length": config.alignment.max_length
        },
        "qc": {
            "min_quality": config.qc.min_quality,
            "trim_quality": config.qc.trim_quality,
            "adapter_sequence": config.qc.adapter_sequence,
            "min_length": config.qc.min_length,
            "max_length": config.qc.max_length
        },
        "counting": {
            "min_mapq": config.counting.min_mapq,
            "count_mode": config.counting.count_mode,
            "strand_specific": config.counting.strand_specific,
            "min_counts_per_sample": config.counting.min_counts_per_sample,
            "min_samples_detected": config.counting.min_samples_detected
        },
        "deseq": {
            "min_cpm": config.deseq.min_cpm,
            "min_samples": config.deseq.min_samples,
            "padj_threshold": config.deseq.padj_threshold,
            "lfc_threshold": config.deseq.lfc_threshold,
            "use_sva": config.deseq.use_sva,
            "normalization_method": config.deseq.normalization_method
        },
        "enrichment": {
            "go_ontologies": config.enrichment.go_ontologies,
            "use_kegg": config.enrichment.use_kegg,
            "enrichment_padj": config.enrichment.enrichment_padj,
            "min_gene_set_size": config.enrichment.min_gene_set_size,
            "max_gene_set_size": config.enrichment.max_gene_set_size
        }
    }


def apply_imported_settings(settings: Dict[str, Any]):
    """Apply imported settings to config"""
    if 'alignment' in settings:
        for key, value in settings['alignment'].items():
            if hasattr(config.alignment, key):
                setattr(config.alignment, key, value)

    if 'qc' in settings:
        for key, value in settings['qc'].items():
            if hasattr(config.qc, key):
                setattr(config.qc, key, value)

    if 'counting' in settings:
        for key, value in settings['counting'].items():
            if hasattr(config.counting, key):
                setattr(config.counting, key, value)

    if 'deseq' in settings:
        for key, value in settings['deseq'].items():
            if hasattr(config.deseq, key):
                setattr(config.deseq, key, value)

    if 'enrichment' in settings:
        for key, value in settings['enrichment'].items():
            if hasattr(config.enrichment, key):
                setattr(config.enrichment, key, value)

    st.session_state['current_preset'] = None


def generate_methods_text() -> str:
    """Generate methods section text for publications"""
    preset = st.session_state.get('current_preset', None)
    preset_text = f" using the '{PRESETS[preset]['name']}' preset" if preset else ""

    return f"""Small RNA-seq analysis was performed using sRNAtlas {config.version}{preset_text}.
Reads were aligned to the reference using Bowtie (v1.x) with {config.alignment.mismatches} mismatch(es) allowed (-v {config.alignment.mismatches}),
reporting up to {config.alignment.max_alignments} alignments per read (-k {config.alignment.max_alignments}){', with --best --strata flags' if config.alignment.best_mode and config.alignment.strata else ''}.
Reads between {config.alignment.min_length}-{config.alignment.max_length} nucleotides were retained.
Differential expression analysis was performed using pyDESeq2 with median ratio normalization.
Genes with CPM < {config.deseq.min_cpm} in fewer than {config.deseq.min_samples} samples were filtered.
Significance was defined as adjusted p-value < {config.deseq.padj_threshold} and |log2FC| > {config.deseq.lfc_threshold}.
GO enrichment analysis was performed using g:Profiler with FDR < {config.enrichment.enrichment_padj}."""
