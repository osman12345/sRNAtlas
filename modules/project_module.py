"""
Project Management Module for sRNAtlas
Save, load, and manage analysis projects
"""
import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional, Any
from datetime import datetime
import json
import pickle
import zipfile
import io
import hashlib
import shutil

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from config import config


# Project file format version
PROJECT_VERSION = "1.0"

# Session state keys to save
SAVEABLE_KEYS = [
    # Data
    'count_matrix',
    'annotated_counts',
    'sample_stats',
    'metadata',
    'annotations',

    # DE results
    'de_results',
    'de_comparisons',

    # Enrichment
    'enrichment_genes',
    'enrichment_results',
    'enrichment_background',

    # Target prediction
    'target_results',
    'target_mirna_seqs',

    # Trimming
    'trim_stats',

    # Settings
    'analysis_settings',
    'active_reference',

    # Project info
    'project_name',
    'project_description',
]


def serialize_dataframe(df: pd.DataFrame) -> Dict:
    """Serialize DataFrame to JSON-compatible format"""
    return {
        'type': 'dataframe',
        'data': df.to_json(orient='split', date_format='iso'),
        'index_name': df.index.name
    }


def deserialize_dataframe(data: Dict) -> pd.DataFrame:
    """Deserialize DataFrame from JSON format"""
    df = pd.read_json(io.StringIO(data['data']), orient='split')
    if data.get('index_name'):
        df.index.name = data['index_name']
    return df


def serialize_value(value: Any) -> Any:
    """Serialize a value for JSON storage"""
    if isinstance(value, pd.DataFrame):
        return serialize_dataframe(value)
    elif isinstance(value, np.ndarray):
        return {
            'type': 'ndarray',
            'data': value.tolist(),
            'dtype': str(value.dtype)
        }
    elif isinstance(value, (np.int64, np.int32)):
        return int(value)
    elif isinstance(value, (np.float64, np.float32)):
        return float(value)
    elif isinstance(value, dict):
        return {
            'type': 'dict',
            'data': {k: serialize_value(v) for k, v in value.items()}
        }
    elif isinstance(value, list):
        return {
            'type': 'list',
            'data': [serialize_value(v) for v in value]
        }
    elif isinstance(value, Path):
        return {
            'type': 'path',
            'data': str(value)
        }
    else:
        return value


def deserialize_value(value: Any) -> Any:
    """Deserialize a value from JSON storage"""
    if isinstance(value, dict) and 'type' in value:
        if value['type'] == 'dataframe':
            return deserialize_dataframe(value)
        elif value['type'] == 'ndarray':
            return np.array(value['data'], dtype=value['dtype'])
        elif value['type'] == 'dict':
            return {k: deserialize_value(v) for k, v in value['data'].items()}
        elif value['type'] == 'list':
            return [deserialize_value(v) for v in value['data']]
        elif value['type'] == 'path':
            return Path(value['data'])
    return value


def create_project_snapshot() -> Dict:
    """Create a snapshot of current session state"""
    snapshot = {
        'version': PROJECT_VERSION,
        'created_at': datetime.now().isoformat(),
        'data': {}
    }

    for key in SAVEABLE_KEYS:
        if key in st.session_state:
            value = st.session_state[key]
            if value is not None:
                try:
                    snapshot['data'][key] = serialize_value(value)
                except Exception as e:
                    st.warning(f"Could not serialize {key}: {e}")

    return snapshot


def restore_project_snapshot(snapshot: Dict):
    """Restore session state from a snapshot"""
    if snapshot.get('version') != PROJECT_VERSION:
        st.warning(f"Project version mismatch. Some data may not load correctly.")

    for key, value in snapshot.get('data', {}).items():
        try:
            st.session_state[key] = deserialize_value(value)
        except Exception as e:
            st.warning(f"Could not restore {key}: {e}")


def save_project_to_file(
    filepath: Path,
    include_raw_data: bool = False
) -> Dict:
    """
    Save project to a file

    Args:
        filepath: Output file path (.srna extension)
        include_raw_data: Whether to include raw FASTQ/BAM files

    Returns:
        Dictionary with save status
    """
    try:
        snapshot = create_project_snapshot()

        # Add metadata
        snapshot['metadata'] = {
            'project_name': st.session_state.get('project_name', 'Untitled'),
            'project_description': st.session_state.get('project_description', ''),
            'saved_at': datetime.now().isoformat(),
            'tool_version': config.version
        }

        # Save as JSON
        filepath = Path(filepath)
        if not filepath.suffix:
            filepath = filepath.with_suffix('.srna')

        with open(filepath, 'w') as f:
            json.dump(snapshot, f, indent=2, default=str)

        return {
            'status': 'success',
            'filepath': str(filepath),
            'size_kb': filepath.stat().st_size / 1024
        }

    except Exception as e:
        return {
            'status': 'error',
            'error': str(e)
        }


def load_project_from_file(filepath: Path) -> Dict:
    """
    Load project from a file

    Args:
        filepath: Input file path

    Returns:
        Dictionary with load status
    """
    try:
        filepath = Path(filepath)

        with open(filepath, 'r') as f:
            snapshot = json.load(f)

        # Restore session state
        restore_project_snapshot(snapshot)

        return {
            'status': 'success',
            'project_name': snapshot.get('metadata', {}).get('project_name', 'Loaded Project'),
            'created_at': snapshot.get('created_at'),
            'keys_loaded': len(snapshot.get('data', {}))
        }

    except json.JSONDecodeError:
        return {
            'status': 'error',
            'error': 'Invalid project file format'
        }
    except Exception as e:
        return {
            'status': 'error',
            'error': str(e)
        }


def export_project_zip(
    include_results: bool = True,
    include_figures: bool = True,
    include_config: bool = True
) -> io.BytesIO:
    """
    Export project as a ZIP archive

    Args:
        include_results: Include result CSV files
        include_figures: Include generated figures
        include_config: Include configuration files

    Returns:
        BytesIO buffer containing ZIP file
    """
    buffer = io.BytesIO()

    with zipfile.ZipFile(buffer, 'w', zipfile.ZIP_DEFLATED) as zf:
        # Project file
        snapshot = create_project_snapshot()
        zf.writestr('project.srna', json.dumps(snapshot, indent=2, default=str))

        # README
        readme = f"""# sRNA-seq Analysis Project Export

Project: {st.session_state.get('project_name', 'Untitled')}
Exported: {datetime.now().isoformat()}
Tool Version: {config.version}

## Contents

- project.srna: Project state file (can be loaded back into the tool)
"""

        # Results
        if include_results:
            readme += "- results/: Analysis result files\n"

            # Count matrix
            if 'count_matrix' in st.session_state:
                counts = st.session_state.count_matrix
                zf.writestr('results/count_matrix.csv', counts.to_csv())

            # Annotated counts
            if 'annotated_counts' in st.session_state:
                ann_counts = st.session_state.annotated_counts
                zf.writestr('results/annotated_counts.csv', ann_counts.to_csv())

            # DE results
            if 'de_results' in st.session_state:
                de_results = st.session_state.de_results
                for comp, df in de_results.get('results', {}).items():
                    if isinstance(df, pd.DataFrame):
                        zf.writestr(f'results/de_{comp}.csv', df.to_csv())

            # Enrichment results
            if 'enrichment_results' in st.session_state:
                enrich = st.session_state.enrichment_results
                for source, df in enrich.get('results', {}).items():
                    if isinstance(df, pd.DataFrame):
                        zf.writestr(f'results/enrichment_{source}.csv', df.to_csv(index=False))

            # Target predictions
            if 'target_results' in st.session_state:
                targets = st.session_state.target_results
                if isinstance(targets, pd.DataFrame):
                    zf.writestr('results/target_predictions.csv', targets.to_csv(index=False))

        # Configuration
        if include_config:
            readme += "- config/: Configuration files\n"

            # Analysis settings
            if 'analysis_settings' in st.session_state:
                settings = st.session_state.analysis_settings
                zf.writestr('config/analysis_settings.json', json.dumps(settings, indent=2))

            # Metadata
            if 'metadata' in st.session_state:
                meta = st.session_state.metadata
                if isinstance(meta, pd.DataFrame):
                    zf.writestr('config/sample_metadata.csv', meta.to_csv(index=False))

        # Write README
        zf.writestr('README.md', readme)

    buffer.seek(0)
    return buffer


def list_recent_projects(projects_dir: Path = None) -> List[Dict]:
    """List recently saved projects"""
    if projects_dir is None:
        projects_dir = Path(config.paths.output_dir) / "projects"

    if not projects_dir.exists():
        return []

    projects = []

    for proj_file in projects_dir.glob("*.srna"):
        try:
            with open(proj_file, 'r') as f:
                data = json.load(f)

            projects.append({
                'filepath': proj_file,
                'name': data.get('metadata', {}).get('project_name', proj_file.stem),
                'saved_at': data.get('metadata', {}).get('saved_at', 'Unknown'),
                'size_kb': proj_file.stat().st_size / 1024
            })
        except Exception:
            pass

    # Sort by modification time
    projects.sort(key=lambda p: p['saved_at'], reverse=True)

    return projects


def render_project_page():
    """Render the project management page"""
    st.header("📁 Project Management")

    st.markdown("""
    Upload your data and manage analysis projects.
    Start by uploading raw FASTQ files in the **Data Upload** tab.
    """)

    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "📤 Data Upload",
        "💾 Save",
        "📂 Load",
        "📥 Export",
        "⚙️ Settings"
    ])

    with tab1:
        render_data_upload()

    with tab2:
        render_save_project()

    with tab3:
        render_load_project()

    with tab4:
        render_export_project()

    with tab5:
        render_project_settings()


def render_data_upload():
    """Render data upload section for raw FASTQ files"""
    st.subheader("Upload Raw Data")

    st.markdown("""
    Upload your raw FASTQ files here. These files will be available for all subsequent
    analysis steps (Quality Control, Trimming, Alignment, etc.).

    **Supported formats:** `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`
    """)

    # Project name
    project_name = st.text_input(
        "Project Name",
        value=st.session_state.get('project_name', 'My Analysis'),
        help="Name for this analysis project",
        key="project_name_upload"
    )
    st.session_state.project_name = project_name

    st.divider()

    # File uploader - accepts all types, filters manually
    st.markdown("### Upload FASTQ Files")

    uploaded_files = st.file_uploader(
        "Select FASTQ files",
        type=None,  # Accept all - filter manually for .fastq.gz support
        accept_multiple_files=True,
        help="Upload raw FASTQ files (.fastq, .fq, .fastq.gz, .fq.gz) - up to 2GB per file",
        key="project_fastq_uploader"
    )

    # Filter to valid FASTQ files
    if uploaded_files:
        valid_extensions = ('.fastq', '.fq', '.fastq.gz', '.fq.gz')
        valid_files = []
        invalid_files = []

        for f in uploaded_files:
            fname = f.name.lower()
            if fname.endswith(valid_extensions):
                valid_files.append(f)
            else:
                invalid_files.append(f.name)

        if invalid_files:
            st.warning(f"⚠️ Skipped {len(invalid_files)} non-FASTQ files: {', '.join(invalid_files[:3])}...")

        uploaded_files = valid_files if valid_files else None

    if uploaded_files:
        st.success(f"✅ Selected {len(uploaded_files)} FASTQ files")

        # Show file details
        file_data = []
        total_size = 0
        for f in uploaded_files:
            size_mb = f.size / (1024 * 1024)
            total_size += size_mb
            file_data.append({
                'Filename': f.name,
                'Size (MB)': f"{size_mb:.2f}",
                'Compressed': '✓' if f.name.endswith('.gz') else ''
            })

        st.dataframe(pd.DataFrame(file_data), width="stretch", key="upload_file_table")
        st.caption(f"**Total size:** {total_size:.2f} MB")

        # Save files to session
        st.divider()

        if st.button("📁 Add Files to Project", type="primary", width="stretch", key="add_files_btn"):
            # Store file info in session state
            project_files = st.session_state.get('project_fastq_files', [])

            # Save files temporarily
            import tempfile
            temp_dir = Path(tempfile.mkdtemp())
            saved_paths = []

            progress = st.progress(0)
            for i, uf in enumerate(uploaded_files):
                progress.progress((i + 1) / len(uploaded_files))

                file_path = temp_dir / uf.name
                with open(file_path, 'wb') as f:
                    f.write(uf.getbuffer())
                saved_paths.append(file_path)

            # Store in session state
            st.session_state.project_fastq_files = saved_paths
            st.session_state.project_fastq_dir = temp_dir
            st.session_state.project_fastq_names = [f.name for f in uploaded_files]

            st.success(f"✅ Added {len(saved_paths)} files to project!")
            st.info("Files are now available in **Quality Control**, **Trimming**, and other modules.")
            st.rerun()

    # Show currently loaded files
    st.divider()
    st.markdown("### Current Project Files")

    project_files = st.session_state.get('project_fastq_files', [])
    project_names = st.session_state.get('project_fastq_names', [])

    if project_files:
        st.success(f"📁 **{len(project_files)} files loaded**")

        # Show file list
        for i, (path, name) in enumerate(zip(project_files, project_names)):
            col1, col2 = st.columns([4, 1])
            with col1:
                if isinstance(path, Path) and path.exists():
                    size_mb = path.stat().st_size / (1024 * 1024)
                    st.caption(f"📄 {name} ({size_mb:.2f} MB)")
                else:
                    st.caption(f"📄 {name}")

        # Option to clear files
        if st.button("🗑️ Clear Project Files", key="clear_files_btn"):
            # Clean up temp directory
            temp_dir = st.session_state.get('project_fastq_dir')
            if temp_dir and Path(temp_dir).exists():
                import shutil
                shutil.rmtree(temp_dir, ignore_errors=True)

            st.session_state.project_fastq_files = []
            st.session_state.project_fastq_names = []
            st.session_state.project_fastq_dir = None
            st.success("Project files cleared")
            st.rerun()

        # Quick actions
        st.divider()
        st.markdown("### Quick Actions")

        col1, col2, col3 = st.columns(3)

        with col1:
            if st.button("📊 Run QC", width="stretch", key="quick_qc_btn"):
                st.info("Go to **Quality Control** module to analyze these files")

        with col2:
            if st.button("✂️ Trim Adapters", width="stretch", key="quick_trim_btn"):
                st.info("Go to **Trimming** module to remove adapters")

        with col3:
            if st.button("🔗 Align Reads", width="stretch", key="quick_align_btn"):
                st.info("Go to **Alignment** module after trimming")

    else:
        st.info("No files loaded. Upload FASTQ files above to start your analysis.")


def render_save_project():
    """Render save project section"""
    st.subheader("Save Project")

    # Project info
    project_name = st.text_input(
        "Project Name",
        value=st.session_state.get('project_name', 'My Analysis'),
        help="Name for this analysis project",
        key="project_name_save"
    )
    st.session_state.project_name = project_name

    project_desc = st.text_area(
        "Description (optional)",
        value=st.session_state.get('project_description', ''),
        help="Brief description of the analysis"
    )
    st.session_state.project_description = project_desc

    # Current state summary
    st.markdown("### Current Session Data")

    data_summary = []

    if 'count_matrix' in st.session_state and st.session_state.count_matrix is not None:
        cm = st.session_state.count_matrix
        if hasattr(cm, 'shape'):
            data_summary.append(f"✅ Count Matrix: {cm.shape[0]} features × {cm.shape[1]} samples")

    if 'de_results' in st.session_state and st.session_state.de_results is not None:
        de_res = st.session_state.de_results
        n_comp = len(de_res.get('results', {})) if isinstance(de_res, dict) else 0
        if n_comp > 0:
            data_summary.append(f"✅ DE Results: {n_comp} comparisons")

    if 'enrichment_results' in st.session_state and st.session_state.enrichment_results is not None:
        enr_res = st.session_state.enrichment_results
        n_sources = len(enr_res.get('results', {})) if isinstance(enr_res, dict) else 0
        if n_sources > 0:
            data_summary.append(f"✅ Enrichment: {n_sources} sources")

    if 'target_results' in st.session_state and st.session_state.target_results is not None:
        tr = st.session_state.target_results
        n_targets = len(tr) if isinstance(tr, pd.DataFrame) else 0
        if n_targets > 0:
            data_summary.append(f"✅ Target Predictions: {n_targets}")

    if 'trim_stats' in st.session_state:
        data_summary.append("✅ Trimming Statistics")

    if not data_summary:
        st.info("No analysis data in current session")
    else:
        for item in data_summary:
            st.write(item)

    # Save options
    st.divider()

    projects_dir = Path(config.paths.output_dir) / "projects"
    projects_dir.mkdir(parents=True, exist_ok=True)

    # Generate filename - handle None or empty project_name
    if not project_name:
        project_name = "My_Analysis"
    safe_name = "".join(c if c.isalnum() or c in '-_' else '_' for c in project_name)
    default_filename = f"{safe_name}_{datetime.now().strftime('%Y%m%d')}.srna"

    col1, col2 = st.columns([3, 1])

    with col1:
        save_path = st.text_input(
            "Save Location",
            value=str(projects_dir / default_filename)
        )

    with col2:
        st.write("")  # Spacer
        st.write("")

    if st.button("💾 Save Project", type="primary", width="stretch"):
        result = save_project_to_file(Path(save_path))

        if result['status'] == 'success':
            st.success(f"✅ Project saved: {result['filepath']}")
            st.info(f"Size: {result['size_kb']:.1f} KB")
        else:
            st.error(f"❌ Save failed: {result['error']}")


def render_load_project():
    """Render load project section"""
    st.subheader("Load Project")

    # Recent projects
    st.markdown("### Recent Projects")

    recent = list_recent_projects()

    if recent:
        for proj in recent[:5]:
            col1, col2, col3 = st.columns([3, 1, 1])

            with col1:
                st.write(f"**{proj['name']}**")
                st.caption(f"Saved: {proj['saved_at'][:10]}")

            with col2:
                st.caption(f"{proj['size_kb']:.1f} KB")

            with col3:
                if st.button("Load", key=f"load_{proj['filepath']}"):
                    result = load_project_from_file(proj['filepath'])

                    if result['status'] == 'success':
                        st.success(f"✅ Loaded: {result['project_name']}")
                        st.info(f"Restored {result['keys_loaded']} data items")
                        st.rerun()
                    else:
                        st.error(f"❌ Load failed: {result['error']}")
    else:
        st.info("No recent projects found")

    # Upload project file
    st.divider()
    st.markdown("### Upload Project File")

    st.caption("💡 To upload **FASTQ files**, use the **Trimming** or **Quality Control** modules.")

    uploaded_file = st.file_uploader(
        "Upload saved project file (.srna or .json)",
        type=['srna', 'json'],
        help="Upload a previously saved project file to restore your analysis session"
    )

    if uploaded_file:
        # Save to temp location and load
        temp_path = Path(config.paths.output_dir) / "temp_project.srna"
        temp_path.parent.mkdir(parents=True, exist_ok=True)

        with open(temp_path, 'wb') as f:
            f.write(uploaded_file.getbuffer())

        if st.button("📂 Load Uploaded Project"):
            result = load_project_from_file(temp_path)

            if result['status'] == 'success':
                st.success(f"✅ Loaded: {result['project_name']}")
                st.rerun()
            else:
                st.error(f"❌ Load failed: {result['error']}")

            # Clean up
            temp_path.unlink()


def render_export_project():
    """Render export project section"""
    st.subheader("Export Project")

    st.markdown("""
    Export your project as a ZIP archive containing all results,
    figures, and configuration files.
    """)

    # Export options
    st.markdown("### Export Options")

    include_results = st.checkbox("Include result files (CSV)", value=True)
    include_figures = st.checkbox("Include figures", value=True)
    include_config = st.checkbox("Include configuration", value=True)

    st.divider()

    if st.button("📤 Generate Export", type="primary", width="stretch"):
        with st.spinner("Generating export..."):
            buffer = export_project_zip(
                include_results=include_results,
                include_figures=include_figures,
                include_config=include_config
            )

        project_name = st.session_state.get('project_name', 'analysis')
        safe_name = "".join(c if c.isalnum() or c in '-_' else '_' for c in project_name)

        st.download_button(
            "📥 Download ZIP Archive",
            buffer.getvalue(),
            f"{safe_name}_export.zip",
            "application/zip"
        )


def render_project_settings():
    """Render project settings section"""
    st.subheader("Project Settings")

    # Auto-save
    st.markdown("### Auto-Save")

    auto_save = st.checkbox(
        "Enable auto-save",
        value=st.session_state.get('auto_save_enabled', False),
        help="Automatically save project every 5 minutes"
    )
    st.session_state.auto_save_enabled = auto_save

    if auto_save:
        st.info("Auto-save will save to the projects directory every 5 minutes")

    # Clear session
    st.divider()
    st.markdown("### Session Management")

    if st.button("🗑️ Clear Current Session", type="secondary"):
        confirm = st.checkbox("I understand this will clear all unsaved data")

        if confirm:
            for key in SAVEABLE_KEYS:
                if key in st.session_state:
                    del st.session_state[key]
            st.success("Session cleared")
            st.rerun()

    # Storage info
    st.divider()
    st.markdown("### Storage Information")

    projects_dir = Path(config.paths.output_dir) / "projects"

    if projects_dir.exists():
        total_size = sum(f.stat().st_size for f in projects_dir.glob("*.srna"))
        n_projects = len(list(projects_dir.glob("*.srna")))

        col1, col2 = st.columns(2)

        with col1:
            st.metric("Saved Projects", n_projects)

        with col2:
            st.metric("Total Size", f"{total_size / 1024:.1f} KB")

        # Cleanup old projects
        if n_projects > 10:
            st.warning(f"You have {n_projects} saved projects. Consider cleaning up old ones.")

            if st.button("🧹 Clean Up Old Projects"):
                # Keep only 10 most recent
                projects = list_recent_projects()
                for proj in projects[10:]:
                    Path(proj['filepath']).unlink()
                st.success(f"Removed {len(projects) - 10} old projects")
                st.rerun()
