"""
sRNAtlas - Main Application
Comprehensive Small RNA-seq Analysis Platform
"""
import streamlit as st
from streamlit_option_menu import option_menu
import sys
import base64
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from config import config

# Logo path
LOGO_PATH = Path(__file__).parent.parent / "assets" / "logo.svg"

def get_logo_base64():
    """Read and encode logo as base64"""
    if LOGO_PATH.exists():
        with open(LOGO_PATH, "r") as f:
            svg_content = f.read()
        return base64.b64encode(svg_content.encode()).decode()
    return None

# Get version
def get_version():
    """Read version from VERSION file"""
    version_file = Path(__file__).parent.parent / "VERSION"
    if version_file.exists():
        return version_file.read_text().strip()
    return "dev"

__version__ = get_version()

# Page configuration
st.set_page_config(
    page_title=config.app_name,
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        color: #1E88E5;
        text-align: center;
        padding: 1rem 0;
    }
    .sub-header {
        font-size: 1.2rem;
        color: #666;
        text-align: center;
        margin-bottom: 2rem;
    }
    .metric-card {
        background-color: #f8f9fa;
        border-radius: 10px;
        padding: 1rem;
        margin: 0.5rem 0;
    }
    .stProgress > div > div > div > div {
        background-color: #1E88E5;
    }
    .sidebar .sidebar-content {
        background-color: #f0f2f6;
    }
</style>
""", unsafe_allow_html=True)


def init_session_state():
    """Initialize session state variables"""
    from utils.provenance import ProvenanceTracker

    defaults = {
        # Project info
        'project_name': None,
        'project_dir': None,
        'organism': None,
        'taxid': None,

        # Data
        'uploaded_files': [],
        'sample_metadata': None,
        'reference_fasta': None,
        'reference_index': None,

        # Results
        'qc_results': None,
        'alignment_results': None,
        'count_matrix': None,
        'annotations': None,
        'de_results': None,
        'enrichment_results': None,

        # Job tracking
        'current_job': None,
        'job_history': [],

        # Settings
        'alignment_settings': config.alignment,
        'qc_settings': config.qc,
        'counting_settings': config.counting,
        'deseq_settings': config.deseq,
        'enrichment_settings': config.enrichment,

        # Provenance tracking
        'provenance_tracker': ProvenanceTracker("sRNAtlas", __version__),
    }

    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value


def render_sidebar():
    """Render the sidebar navigation"""
    with st.sidebar:
        # Display logo
        logo_b64 = get_logo_base64()
        if logo_b64:
            st.markdown(f'''
                <div style="display: flex; align-items: center; gap: 10px; margin-bottom: 10px;">
                    <img src="data:image/svg+xml;base64,{logo_b64}" width="40" height="40">
                    <span style="font-size: 1.3rem; font-weight: bold; color: #1E88E5;">sRNAtlas</span>
                </div>
            ''', unsafe_allow_html=True)
        else:
            st.markdown("### 🧬 sRNAtlas")

        selected = option_menu(
            menu_title="Navigation",
            options=[
                "Home",
                "Project",
                "Quality Control",
                "Trimming",
                "Databases",
                "Alignment",
                "Post-Align QC",
                "Counting",
                "DE Analysis",
                "Novel miRNA",
                "isomiR",
                "Targets",
                "GO/Pathway",
                "Batch",
                "Reports",
                "Settings",
                "Help"
            ],
            icons=[
                "house", "folder", "bar-chart", "scissors",
                "database", "link", "clipboard-check", "graph-up", "search",
                "stars", "layers", "bullseye", "diagram-3", "lightning",
                "file-text", "gear", "question-circle"
            ],
            default_index=0,
            styles={
                "container": {"padding": "5px"},
                "icon": {"color": "#1E88E5", "font-size": "18px"},
                "nav-link": {"font-size": "14px", "text-align": "left", "margin": "0px"},
                "nav-link-selected": {"background-color": "#1E88E5"},
            }
        )

        # Show current project info
        st.divider()
        if st.session_state.project_name:
            st.success(f"**Project:** {st.session_state.project_name}")
            if st.session_state.organism:
                st.info(f"**Organism:** {st.session_state.organism}")
        else:
            st.warning("No project loaded")

        # Show job status
        if st.session_state.current_job:
            st.divider()
            st.info(f"**Active Job:** {st.session_state.current_job}")

        # Show version at bottom
        st.divider()
        st.caption(f"Version {__version__}")

        return selected


def render_home():
    """Render the home page"""
    # Display logo and title
    logo_b64 = get_logo_base64()
    if logo_b64:
        st.markdown(f'''
            <div style="text-align: center; padding: 1rem 0;">
                <img src="data:image/svg+xml;base64,{logo_b64}" width="80" height="80" style="margin-bottom: 10px;">
                <p class="main-header" style="margin: 0;">sRNAtlas</p>
            </div>
        ''', unsafe_allow_html=True)
    else:
        st.markdown('<p class="main-header">🧬 sRNAtlas</p>', unsafe_allow_html=True)
    st.markdown('<p class="sub-header">Comprehensive Small RNA-seq Analysis Platform</p>', unsafe_allow_html=True)

    # Pipeline overview
    st.subheader("📋 Analysis Pipeline")

    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown("""
        ### ✂️ Adapter Trimming
        - Cutadapt integration
        - Common adapter presets
        - Quality filtering
        - Length filtering
        """)

        st.markdown("""
        ### 📊 Quality Control
        - Read quality assessment
        - Size distribution analysis
        - Adapter content detection
        - rRNA contamination check
        """)

    with col2:
        st.markdown("""
        ### 🗄️ Reference Databases
        - miRBase download
        - RNAcentral support
        - Custom database builder
        - Automatic indexing
        """)

        st.markdown("""
        ### 🔗 Alignment & Counting
        - Bowtie alignment (optimized for small RNA)
        - Multi-mapper handling
        - RNA type classification
        - Count matrix generation
        """)

    with col3:
        st.markdown("""
        ### 🔬 DE & Target Analysis
        - pyDESeq2 analysis
        - miRNA target prediction
        - psRNATarget integration
        - GO/KEGG enrichment
        """)

        st.markdown("""
        ### ⚡ Batch & Reports
        - Full pipeline automation
        - Project save/load
        - HTML/PDF reports
        - Result export
        """)

    st.divider()

    # Quick start
    st.subheader("🚀 Quick Start")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("""
        **New Analysis:**
        1. **Trimming**: Remove adapters with Cutadapt
        2. **QC**: Check read quality
        3. **Database**: Download/build reference
        4. **Alignment**: Map reads with Bowtie
        5. **Counting**: Generate count matrix
        6. **DE Analysis**: Find differential expression
        7. **Targets**: Predict miRNA targets
        8. **Enrichment**: GO/KEGG analysis
        """)

    with col2:
        st.markdown("""
        **Project Management:**
        - 💾 **Save** your analysis progress
        - 📂 **Load** previous projects
        - 📤 **Export** results as ZIP
        - ⚡ **Batch** process multiple samples

        **Supported Input:**
        - Raw FASTQ files
        - Pre-aligned BAM files
        - Count matrices (CSV)
        """)

    # Supported features
    st.divider()
    st.subheader("✨ Key Features")

    features = [
        ("🧬 Multi-organism", "Plants & animals with miRBase/RNAcentral"),
        ("✂️ Adapter Trimming", "Cutadapt with common presets"),
        ("🎯 Target Prediction", "psRNATarget & local seed matching"),
        ("📊 Interactive Plots", "Plotly-powered visualizations"),
        ("⚡ Batch Processing", "Automate full pipeline"),
        ("💾 Project Management", "Save/load analysis sessions"),
        ("🔬 Robust Statistics", "pyDESeq2 with global FDR"),
        ("📥 Flexible Input", "FASTQ, BAM, or count matrices"),
        ("📋 Reports", "HTML reports & figure export"),
    ]

    cols = st.columns(3)
    for i, (title, desc) in enumerate(features):
        with cols[i % 3]:
            st.markdown(f"**{title}**")
            st.caption(desc)


def render_placeholder(title: str, description: str):
    """Render a placeholder page for modules not yet implemented"""
    st.header(title)
    st.info(description)
    st.markdown("This module will be implemented in the full version.")


def main():
    """Main application entry point"""
    # Initialize session state
    init_session_state()

    # Render sidebar and get selection
    selected = render_sidebar()

    # Route to appropriate page
    if "Home" in selected:
        render_home()
    elif "Project" in selected:
        from modules.project_module import render_project_page
        render_project_page()
    elif "Quality Control" in selected:
        from modules.qc_module import render_qc_page
        render_qc_page()
    elif "Trimming" in selected:
        from modules.trimming_module import render_trimming_page
        render_trimming_page()
    elif "Databases" in selected:
        from modules.database_module import render_database_page
        render_database_page()
    elif "Alignment" in selected and "Post" not in selected:
        from modules.alignment_module import render_alignment_page
        render_alignment_page()
    elif "Post-Align" in selected:
        from modules.post_alignment_qc_module import render_post_alignment_qc_page
        render_post_alignment_qc_page()
    elif "Counting" in selected:
        from modules.counting_module import render_counting_page
        render_counting_page()
    elif "DE Analysis" in selected:
        from modules.de_module import render_de_page
        render_de_page()
    elif "Novel miRNA" in selected:
        from modules.novel_mirna_module import render_novel_mirna_page
        render_novel_mirna_page()
    elif "isomiR" in selected:
        from modules.isomir_module import render_isomir_page
        render_isomir_page()
    elif "Targets" in selected:
        from modules.target_prediction_module import render_target_prediction_page
        render_target_prediction_page()
    elif "GO/Pathway" in selected:
        from modules.enrichment_module import render_enrichment_page
        render_enrichment_page()
    elif "Batch" in selected:
        from modules.batch_module import render_batch_page
        render_batch_page()
    elif "Reports" in selected:
        from modules.reports_module import render_reports_page
        render_reports_page()
    elif "Settings" in selected:
        from modules.settings_module import render_settings_page
        render_settings_page()
    elif "Help" in selected:
        from modules.help_module import render_help_page
        render_help_page()


if __name__ == "__main__":
    main()
