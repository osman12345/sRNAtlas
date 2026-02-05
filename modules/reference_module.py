"""
Reference Selection Module for sRNAtlas
UI for selecting and managing reference databases
"""
import streamlit as st
import pandas as pd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.reference_manager import ReferenceManager, get_reference_manager


def render_reference_selector(key_prefix: str = "ref") -> dict:
    """
    Render a reference selector widget

    Args:
        key_prefix: Prefix for session state keys

    Returns:
        Selected reference dictionary or None
    """
    ref_manager = get_reference_manager()
    references = ref_manager.list_references()

    if not references:
        st.warning("No reference databases found. Please add references to `data/references/` folder.")
        return None

    # Create selection options
    options = {ref['id']: f"{ref['name']} ({ref.get('sequences', '?')} sequences)" for ref in references}

    # Add "Custom" option
    options['custom'] = "📁 Upload Custom Reference..."

    selected_id = st.selectbox(
        "Select Reference Database",
        options=list(options.keys()),
        format_func=lambda x: options[x],
        key=f"{key_prefix}_selector",
        help="Choose a pre-loaded reference or upload your own"
    )

    if selected_id == 'custom':
        # Custom upload
        uploaded_file = st.file_uploader(
            "Upload Reference FASTA",
            type=['fasta', 'fa', 'fna'],
            key=f"{key_prefix}_upload",
            help="Upload a FASTA file with reference sequences"
        )

        if uploaded_file:
            st.session_state[f'{key_prefix}_custom_file'] = uploaded_file
            st.success(f"✅ Uploaded: {uploaded_file.name}")
            return {
                'id': 'custom',
                'name': uploaded_file.name,
                'file': uploaded_file,
                'is_custom': True
            }
        return None

    else:
        # Pre-loaded reference
        ref = ref_manager.get_reference(selected_id)
        ref_path = ref_manager.get_reference_path(selected_id)

        if ref and ref_path:
            # Show reference info
            with st.expander("📋 Reference Details", expanded=False):
                col1, col2 = st.columns(2)
                with col1:
                    st.markdown(f"**Species:** {ref.get('species', 'Unknown')}")
                    st.markdown(f"**Sequences:** {ref.get('sequences', 'Unknown'):,}")
                    st.markdown(f"**Source:** {ref.get('source', 'Unknown')}")
                with col2:
                    st.markdown(f"**RNA Types:** {', '.join(ref.get('rna_types', ['Unknown']))}")
                    st.markdown(f"**Version:** {ref.get('version', 'Unknown')}")
                    st.markdown(f"**Taxonomy ID:** {ref.get('taxid', 'Unknown')}")

                if ref.get('description'):
                    st.caption(ref['description'])

            ref['path'] = ref_path
            ref['is_custom'] = False
            return ref

        else:
            st.error(f"Reference file not found: {ref.get('file', 'Unknown')}")
            return None


def render_reference_management_page():
    """Render the full reference management page"""
    st.header("📚 Reference Database Management")

    tab1, tab2, tab3 = st.tabs(["📋 Available References", "➕ Add Reference", "ℹ️ Help"])

    with tab1:
        render_reference_list()

    with tab2:
        render_add_reference()

    with tab3:
        render_reference_help()


def render_reference_list():
    """Render list of available references"""
    st.subheader("Available Reference Databases")

    ref_manager = get_reference_manager()
    references = ref_manager.list_references()

    if not references:
        st.info("No references found. Add references using the 'Add Reference' tab.")
        return

    # Create DataFrame for display
    df_data = []
    for ref in references:
        df_data.append({
            'Name': ref.get('name', 'Unknown'),
            'Species': ref.get('species', 'Unknown'),
            'Sequences': ref.get('sequences', 0),
            'RNA Types': ', '.join(ref.get('rna_types', [])),
            'Source': ref.get('source', 'Unknown'),
            'File': ref.get('file', 'Unknown')
        })

    df = pd.DataFrame(df_data)
    st.dataframe(df, width="stretch", hide_index=True)

    # Validate references
    st.divider()
    validation = ref_manager.validate_references()

    all_valid = all(validation.values())
    if all_valid:
        st.success("✅ All reference files are valid and accessible.")
    else:
        st.warning("⚠️ Some reference files are missing:")
        for ref_id, valid in validation.items():
            if not valid:
                st.error(f"  • {ref_id}: File not found")


def render_add_reference():
    """Render form to add a new reference"""
    st.subheader("Add New Reference")

    st.markdown("""
    **To add a new reference:**
    1. Place your FASTA file in the `data/references/` folder
    2. Fill out the form below
    3. Click "Add Reference"
    """)

    ref_manager = get_reference_manager()

    # List available FASTA files not yet registered
    existing_files = {ref['file'] for ref in ref_manager.list_references()}
    available_files = []

    for ext in ['*.fasta', '*.fa', '*.fna']:
        for f in ref_manager.reference_dir.glob(ext):
            if f.name not in existing_files and not f.name.startswith('.'):
                available_files.append(f.name)

    if not available_files:
        st.info("No unregistered FASTA files found in `data/references/`. Upload a file first.")

        uploaded = st.file_uploader(
            "Upload FASTA file",
            type=['fasta', 'fa', 'fna'],
            key="new_ref_upload"
        )

        if uploaded:
            # Save uploaded file
            save_path = ref_manager.reference_dir / uploaded.name
            with open(save_path, 'wb') as f:
                f.write(uploaded.getvalue())
            st.success(f"✅ Saved {uploaded.name} to references folder")
            st.rerun()

        return

    # Form for adding reference
    with st.form("add_reference_form"):
        col1, col2 = st.columns(2)

        with col1:
            file = st.selectbox("FASTA File", options=available_files)
            name = st.text_input("Display Name", placeholder="e.g., Human miRNAs")
            species = st.text_input("Species", placeholder="e.g., Homo sapiens")

        with col2:
            ref_id = st.text_input(
                "Reference ID",
                value=available_files[0].replace('.fasta', '').replace('.fa', '') if available_files else "",
                help="Unique identifier (no spaces)"
            )
            source = st.selectbox("Source", options=["miRBase", "RNAcentral", "Ensembl", "Custom"])
            taxid = st.number_input("Taxonomy ID (optional)", min_value=0, value=0)

        rna_types = st.multiselect(
            "RNA Types",
            options=["miRNA", "siRNA", "piRNA", "tRNA", "rRNA", "snRNA", "snoRNA", "lncRNA", "other"],
            default=["miRNA"]
        )

        description = st.text_area("Description (optional)", placeholder="Brief description of the reference...")

        submitted = st.form_submit_button("➕ Add Reference", type="primary")

        if submitted:
            if not ref_id or not name or not species:
                st.error("Please fill in all required fields (ID, Name, Species)")
            elif ' ' in ref_id:
                st.error("Reference ID cannot contain spaces")
            else:
                success = ref_manager.add_reference(
                    ref_id=ref_id,
                    name=name,
                    species=species,
                    file=file,
                    description=description,
                    source=source,
                    taxid=taxid if taxid > 0 else None,
                    rna_types=rna_types
                )

                if success:
                    st.success(f"✅ Reference '{name}' added successfully!")
                    st.rerun()
                else:
                    st.error("Failed to add reference. Check that the file exists.")


def render_reference_help():
    """Render help for reference management"""
    st.subheader("Reference Database Help")

    st.markdown("""
    ### About Reference Databases

    Reference databases contain known small RNA sequences for alignment and annotation.
    The tool supports multiple species and custom references.

    ### Supported Formats

    - **FASTA** (.fasta, .fa, .fna)

    ### FASTA Format Example

    ```
    >mtr-miR162_MIMAT0001640_Medicago_truncatula_miR162
    TCGATAAACCTCTGCATCCAG
    >mtr-miR160a_MIMAT0001641_Medicago_truncatula_miR160a
    TGCCTGGCTCCCTGTATGCCA
    ```

    ### Adding New References

    **Method 1: Through the UI**
    1. Go to "Add Reference" tab
    2. Upload your FASTA file
    3. Fill in the metadata
    4. Click "Add Reference"

    **Method 2: Manual**
    1. Place your FASTA file in `data/references/`
    2. Edit `data/references/references.json`
    3. Add an entry for your reference

    ### Example references.json Entry

    ```json
    {
        "id": "hsa_mirbase",
        "name": "Human miRNAs (miRBase v22)",
        "species": "Homo sapiens",
        "common_name": "Human",
        "taxid": 9606,
        "file": "hsa_mirbase_v22.fasta",
        "description": "Human miRNA sequences from miRBase v22",
        "source": "miRBase",
        "version": "22",
        "rna_types": ["miRNA"],
        "sequences": 2654
    }
    ```

    ### Where to Get References

    - **miRBase**: https://www.mirbase.org/ (miRNAs)
    - **RNAcentral**: https://rnacentral.org/ (all ncRNAs)
    - **Ensembl**: https://www.ensembl.org/ (genome annotations)
    - **NCBI**: https://www.ncbi.nlm.nih.gov/ (various sequences)
    """)


# Convenience function for other modules
def get_selected_reference() -> dict:
    """Get the currently selected reference from session state"""
    return st.session_state.get('selected_reference', None)


def set_selected_reference(ref: dict):
    """Set the selected reference in session state"""
    st.session_state['selected_reference'] = ref
