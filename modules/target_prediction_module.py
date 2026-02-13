"""
miRNA Target Prediction Module for sRNAtlas
Integrates psRNATarget (plants) and TargetScan (animals) APIs
"""
import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import requests
import time
import json
import io

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from config import config
from utils.miranda import (
    check_miranda_installed,
    run_miranda,
    filter_miranda_results,
    get_miranda_summary_stats,
    get_miranda_install_instructions
)


# psRNATarget API endpoint
PSRNATARGET_API = "https://www.zhaolab.org/psRNATarget/api"

# Common plant species for psRNATarget
PLANT_SPECIES = {
    "Arabidopsis thaliana": {
        "transcript_db": "Arabidopsis thaliana, TAIR10, cDNA",
        "code": "ath"
    },
    "Medicago truncatula": {
        "transcript_db": "Medicago truncatula, Mt4.0v1, cDNA",
        "code": "mtr"
    },
    "Oryza sativa": {
        "transcript_db": "Oryza sativa, MSU7, cDNA",
        "code": "osa"
    },
    "Zea mays": {
        "transcript_db": "Zea mays, AGPv4, cDNA",
        "code": "zma"
    },
    "Glycine max": {
        "transcript_db": "Glycine max, Wm82.a2.v1, cDNA",
        "code": "gma"
    },
    "Solanum lycopersicum": {
        "transcript_db": "Solanum lycopersicum, ITAG2.4, cDNA",
        "code": "sly"
    },
    "Vitis vinifera": {
        "transcript_db": "Vitis vinifera, Genoscope.12X, cDNA",
        "code": "vvi"
    },
    "Custom": {
        "transcript_db": "Upload your own",
        "code": "custom"
    }
}

# Animal species for TargetScan-like prediction
ANIMAL_SPECIES = {
    "Homo sapiens": "hsa",
    "Mus musculus": "mmu",
    "Rattus norvegicus": "rno",
    "Drosophila melanogaster": "dme",
    "Caenorhabditis elegans": "cel",
    "Danio rerio": "dre"
}


def submit_psrnatarget_job(
    mirna_sequences: Dict[str, str],
    species: str,
    expectation: float = 5.0,
    max_mismatches: int = 4,
    hsp_size: int = 19,
    seed_region: Tuple[int, int] = (2, 13),
    custom_transcripts: Optional[str] = None
) -> Dict:
    """
    Submit a job to psRNATarget server

    Args:
        mirna_sequences: Dictionary of {mirna_id: sequence}
        species: Target species
        expectation: Maximum expectation score (lower = stricter)
        max_mismatches: Maximum allowed mismatches
        hsp_size: Minimum HSP size
        seed_region: Seed region for mismatch penalty (start, end)
        custom_transcripts: Optional custom transcript FASTA

    Returns:
        Dictionary with job ID or error
    """
    # Format miRNA sequences as FASTA
    mirna_fasta = ""
    for mirna_id, seq in mirna_sequences.items():
        mirna_fasta += f">{mirna_id}\n{seq}\n"

    # Prepare submission data
    data = {
        'smallRNA': mirna_fasta,
        'expectation': expectation,
        'penalty_mismatch_seed': 1.0,
        'penalty_mismatch_nonseed': 0.5,
        'penalty_gap_open': 2.0,
        'penalty_gap_ext': 0.5,
        'hsp_size': hsp_size,
        'seed_start': seed_region[0],
        'seed_end': seed_region[1],
    }

    if custom_transcripts:
        data['targetSeq'] = custom_transcripts
    else:
        species_info = PLANT_SPECIES.get(species, PLANT_SPECIES["Arabidopsis thaliana"])
        data['targetDB'] = species_info['transcript_db']

    try:
        response = requests.post(
            f"{PSRNATARGET_API}/submit",
            data=data,
            timeout=60
        )

        if response.status_code == 200:
            result = response.json()
            return {
                'status': 'submitted',
                'job_id': result.get('jobId'),
                'message': 'Job submitted successfully'
            }
        else:
            return {
                'status': 'error',
                'error': f"Server returned status {response.status_code}"
            }

    except requests.exceptions.Timeout:
        return {
            'status': 'error',
            'error': 'Connection timed out'
        }
    except Exception as e:
        return {
            'status': 'error',
            'error': str(e)
        }


def check_psrnatarget_job(job_id: str) -> Dict:
    """Check status of a psRNATarget job"""
    try:
        response = requests.get(
            f"{PSRNATARGET_API}/status/{job_id}",
            timeout=30
        )

        if response.status_code == 200:
            result = response.json()
            return {
                'status': result.get('status', 'unknown'),
                'progress': result.get('progress', 0)
            }
        else:
            return {
                'status': 'error',
                'error': f"Status check failed: {response.status_code}"
            }

    except Exception as e:
        return {
            'status': 'error',
            'error': str(e)
        }


def get_psrnatarget_results(job_id: str) -> Dict:
    """Retrieve results from a completed psRNATarget job"""
    try:
        response = requests.get(
            f"{PSRNATARGET_API}/result/{job_id}",
            timeout=60
        )

        if response.status_code == 200:
            # Parse TSV results
            content = response.text
            df = pd.read_csv(io.StringIO(content), sep='\t')

            return {
                'status': 'success',
                'results': df
            }
        else:
            return {
                'status': 'error',
                'error': f"Failed to retrieve results: {response.status_code}"
            }

    except Exception as e:
        return {
            'status': 'error',
            'error': str(e)
        }


def predict_targets_local(
    mirna_sequences: Dict[str, str],
    transcript_sequences: Dict[str, str],
    seed_start: int = 2,
    seed_end: int = 8,
    max_mismatches: int = 3,
    min_complementarity: float = 0.7
) -> pd.DataFrame:
    """
    Local target prediction using seed matching algorithm

    This is a simplified local implementation for when API is unavailable.
    Uses seed region complementarity scoring.

    Args:
        mirna_sequences: Dictionary of {mirna_id: sequence}
        transcript_sequences: Dictionary of {transcript_id: sequence}
        seed_start: Seed region start (1-indexed)
        seed_end: Seed region end (1-indexed)
        max_mismatches: Maximum mismatches in seed
        min_complementarity: Minimum overall complementarity score

    Returns:
        DataFrame with predicted targets
    """
    def reverse_complement(seq):
        """Get reverse complement of a sequence"""
        complement = {'A': 'U', 'T': 'A', 'U': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(b, 'N') for b in seq.upper()[::-1])

    def count_mismatches(seq1, seq2):
        """Count mismatches between two sequences"""
        return sum(a != b for a, b in zip(seq1, seq2))

    def find_binding_sites(mirna_seq, transcript_seq, seed_start, seed_end):
        """Find potential binding sites in transcript"""
        mirna_rc = reverse_complement(mirna_seq)
        seed = mirna_rc[seed_start-1:seed_end]
        seed_len = len(seed)

        sites = []
        transcript_upper = transcript_seq.upper().replace('T', 'U')

        for i in range(len(transcript_upper) - seed_len + 1):
            site = transcript_upper[i:i+seed_len]
            mismatches = count_mismatches(seed, site)

            if mismatches <= max_mismatches:
                # Calculate extended complementarity
                ext_start = max(0, i - (len(mirna_seq) - seed_end))
                ext_end = min(len(transcript_upper), i + seed_len + seed_start - 1)
                extended_site = transcript_upper[ext_start:ext_end]

                # Score based on complementarity
                if len(extended_site) >= len(mirna_seq) * 0.8:
                    comp_score = 1 - (mismatches / seed_len)
                    sites.append({
                        'position': i + 1,
                        'seed_mismatches': mismatches,
                        'complementarity': comp_score,
                        'target_site': site
                    })

        return sites

    results = []

    for mirna_id, mirna_seq in mirna_sequences.items():
        for trans_id, trans_seq in transcript_sequences.items():
            sites = find_binding_sites(mirna_seq, trans_seq, seed_start, seed_end)

            for site in sites:
                if site['complementarity'] >= min_complementarity:
                    results.append({
                        'miRNA': mirna_id,
                        'Target': trans_id,
                        'Position': site['position'],
                        'Seed_Mismatches': site['seed_mismatches'],
                        'Complementarity': site['complementarity'],
                        'Target_Site': site['target_site']
                    })

    return pd.DataFrame(results)


def load_mirna_from_de_results(
    de_results: Dict,
    padj_threshold: float = 0.05,
    reference_sequences: Optional[Dict[str, str]] = None
) -> Dict[str, str]:
    """
    Extract significant miRNAs from DE results and look up their sequences

    Args:
        de_results: DE results dictionary
        padj_threshold: FDR threshold
        reference_sequences: Dictionary of {mirna_id: sequence} from reference FASTA

    Returns:
        Dictionary of {mirna_id: sequence} for significant miRNAs
    """
    import streamlit as st

    # Get reference sequences from session state if not provided
    if reference_sequences is None:
        reference_sequences = st.session_state.get('reference_sequences', {})

        # Also try mirbase_sequences or loaded FASTA
        if not reference_sequences:
            reference_sequences = st.session_state.get('mirbase_sequences', {})

        # Try to build from reference_fasta if available
        if not reference_sequences and st.session_state.get('reference_fasta'):
            try:
                from Bio import SeqIO
                import io
                fasta_content = st.session_state.get('reference_fasta')
                if isinstance(fasta_content, str):
                    handle = io.StringIO(fasta_content)
                else:
                    handle = io.StringIO(fasta_content.decode('utf-8'))
                reference_sequences = {
                    rec.id: str(rec.seq).upper().replace('T', 'U')
                    for rec in SeqIO.parse(handle, "fasta")
                }
                # Cache for future use
                st.session_state.reference_sequences = reference_sequences
            except Exception as e:
                st.warning(f"Could not parse reference FASTA: {e}")

    significant_mirnas = {}
    missing_sequences = []

    for comparison, df in de_results.get('results', {}).items():
        padj_col = 'padj' if 'padj' in df.columns else 'global_FDR'
        if padj_col in df.columns:
            sig = df[df[padj_col] < padj_threshold].index.tolist()
            for mirna in sig:
                if mirna in significant_mirnas:
                    continue  # Already added

                # Look up sequence in reference
                if mirna in reference_sequences:
                    significant_mirnas[mirna] = reference_sequences[mirna]
                else:
                    # Try fuzzy matching (strip version suffix like .1, -5p, etc.)
                    mirna_base = mirna.split('.')[0].split('-5p')[0].split('-3p')[0]
                    found = False
                    for ref_id, ref_seq in reference_sequences.items():
                        ref_base = ref_id.split('.')[0].split('-5p')[0].split('-3p')[0]
                        if mirna_base.lower() == ref_base.lower():
                            significant_mirnas[mirna] = ref_seq
                            found = True
                            break
                    if not found:
                        missing_sequences.append(mirna)

    # Warn about missing sequences
    if missing_sequences:
        st.warning(
            f"⚠️ Could not find sequences for {len(missing_sequences)} miRNAs. "
            f"Please load reference sequences in the Databases module.\n"
            f"Missing: {', '.join(missing_sequences[:5])}{'...' if len(missing_sequences) > 5 else ''}"
        )

    return significant_mirnas


def render_target_prediction_page():
    """Render the miRNA target prediction page"""
    st.header("🎯 miRNA Target Prediction")

    st.markdown("""
    Predict mRNA targets for your miRNAs using established algorithms.

    **Available methods:**
    - **psRNATarget** (plants): Web API for plant miRNA target prediction
    - **miRanda** (animals): Thermodynamics-based algorithm for animal 3' UTR targeting
    - **Local Seed Matching**: Fast local algorithm based on seed complementarity
    """)

    tab1, tab2, tab3 = st.tabs([
        "📁 Input",
        "🔮 Predict",
        "📊 Results"
    ])

    with tab1:
        render_target_input()

    with tab2:
        render_target_prediction()

    with tab3:
        render_target_results()


def render_target_input():
    """Render input section for target prediction"""
    st.subheader("miRNA Input")

    input_source = st.radio(
        "miRNA Source",
        options=["From DE Results", "Upload miRNA List", "Manual Entry"]
    )

    if input_source == "From DE Results":
        if st.session_state.get('de_results') is None:
            st.warning("No DE results available. Run DE analysis first.")
            return

        results = st.session_state.de_results['results']
        comparison = st.selectbox(
            "Select Comparison",
            options=list(results.keys())
        )

        padj_threshold = st.slider("FDR Threshold", 0.01, 0.20, 0.05, 0.01)
        lfc_threshold = st.slider("LFC Threshold", 0.0, 2.0, 0.585, 0.1)

        direction = st.selectbox(
            "Direction",
            options=["Both", "Up-regulated", "Down-regulated"]
        )

        res_df = results[comparison]
        padj_col = 'padj' if 'padj' in res_df.columns else 'global_FDR'
        lfc_col = 'log2FoldChange' if 'log2FoldChange' in res_df.columns else 'lfc'

        # Filter significant miRNAs
        sig_mask = res_df[padj_col] < padj_threshold

        if direction == "Up-regulated":
            sig_mask = sig_mask & (res_df[lfc_col] > lfc_threshold)
        elif direction == "Down-regulated":
            sig_mask = sig_mask & (res_df[lfc_col] < -lfc_threshold)
        else:
            sig_mask = sig_mask & (res_df[lfc_col].abs() > lfc_threshold)

        sig_mirnas = res_df[sig_mask].index.tolist()
        st.info(f"Selected {len(sig_mirnas)} significant miRNAs")

        if sig_mirnas:
            st.session_state.target_mirnas = sig_mirnas

            with st.expander("View miRNA List"):
                st.text('\n'.join(sig_mirnas[:50]))
                if len(sig_mirnas) > 50:
                    st.text(f"... and {len(sig_mirnas) - 50} more")

    elif input_source == "Upload miRNA List":
        st.markdown("""
        Upload a file with miRNA IDs and sequences in FASTA format or tab-separated format:
        ```
        >miR-156a
        UGACAGAAGAGAGUGAGCAC
        >miR-156b
        UGACAGAAGAGAGUGAGCAU
        ```
        """)

        mirna_file = st.file_uploader(
            "Upload miRNA sequences",
            type=['fasta', 'fa', 'txt', 'csv', 'tsv']
        )

        if mirna_file:
            content = mirna_file.getvalue().decode('utf-8')

            # Parse FASTA or tabular
            mirnas = {}
            if content.startswith('>'):
                # FASTA format
                current_id = None
                current_seq = ""
                for line in content.split('\n'):
                    line = line.strip()
                    if line.startswith('>'):
                        if current_id:
                            mirnas[current_id] = current_seq
                        current_id = line[1:].split()[0]
                        current_seq = ""
                    else:
                        current_seq += line
                if current_id:
                    mirnas[current_id] = current_seq
            else:
                # Tabular format
                for line in content.split('\n'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        mirnas[parts[0]] = parts[1]

            st.info(f"Loaded {len(mirnas)} miRNA sequences")
            st.session_state.target_mirna_seqs = mirnas

    else:  # Manual entry
        mirna_text = st.text_area(
            "Enter miRNA sequences (FASTA format)",
            height=200,
            placeholder=">miR-156a\nUGACAGAAGAGAGUGAGCAC\n>miR-156b\nUGACAGAAGAGAGUGAGCAU"
        )

        if mirna_text:
            mirnas = {}
            current_id = None
            current_seq = ""
            for line in mirna_text.split('\n'):
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        mirnas[current_id] = current_seq
                    current_id = line[1:].split()[0]
                    current_seq = ""
                else:
                    current_seq += line.upper().replace('T', 'U')
            if current_id:
                mirnas[current_id] = current_seq

            st.info(f"Parsed {len(mirnas)} miRNA sequences")
            st.session_state.target_mirna_seqs = mirnas


def render_target_prediction():
    """Render target prediction section"""
    st.subheader("Run Target Prediction")

    # Initialize all parameter variables with defaults to avoid NameError
    # These will be overwritten by the appropriate UI elements based on method
    score_threshold = 140.0
    energy_threshold = -20.0
    miranda_strict = False
    seed_start = 2
    seed_end = 8
    max_mismatches = 2
    expectation = 5.0

    # Check for input
    mirna_seqs = st.session_state.get('target_mirna_seqs', {})
    mirna_ids = st.session_state.get('target_mirnas', [])

    if not mirna_seqs and not mirna_ids:
        st.warning("Please provide miRNA sequences in the Input tab first.")
        return

    # Prediction method
    kingdom = st.radio(
        "Select Kingdom",
        options=["Plants", "Animals"],
        help="Different algorithms are optimized for different kingdoms"
    )

    if kingdom == "Plants":
        method = st.radio(
            "Prediction Method",
            options=["psRNATarget (API)", "Local Seed Matching"],
            help="psRNATarget requires internet; local method is faster but less accurate"
        )

        species = st.selectbox(
            "Target Species",
            options=list(PLANT_SPECIES.keys())
        )

    else:  # Animals
        # Check if miRanda is available
        miranda_available, miranda_version = check_miranda_installed()

        if miranda_available:
            method = st.radio(
                "Prediction Method",
                options=["miRanda (Recommended)", "Local Seed Matching"],
                help="miRanda is a well-established tool for animal miRNA target prediction"
            )
            st.success(f"✅ {miranda_version} detected")
        else:
            method = "Local Seed Matching"
            st.warning("⚠️ miRanda not installed. Using local seed matching.")
            with st.expander("📦 Install miRanda"):
                st.markdown(get_miranda_install_instructions())

        species = st.selectbox(
            "Target Species",
            options=list(ANIMAL_SPECIES.keys())
        )

    # Custom transcripts option
    use_custom = st.checkbox("Use custom transcript sequences")

    if use_custom:
        custom_file = st.file_uploader(
            "Upload transcript FASTA",
            type=['fasta', 'fa', 'txt']
        )

        if custom_file:
            st.session_state.custom_transcripts = custom_file.getvalue().decode('utf-8')

    # Parameters
    st.markdown("### Parameters")

    col1, col2 = st.columns(2)

    with col1:
        if method == "miRanda (Recommended)":
            score_threshold = st.slider("Min Score", 100.0, 200.0, 140.0, 5.0,
                help="Minimum alignment score (higher = stricter)")
            energy_threshold = st.slider("Max Energy (kcal/mol)", -50.0, 0.0, -20.0, 1.0,
                help="Maximum free energy (more negative = stronger binding)")
        else:
            seed_start = st.slider("Seed Start Position", 1, 5, 2)
            seed_end = st.slider("Seed End Position", 7, 13, 8)

    with col2:
        if method == "miRanda (Recommended)":
            miranda_strict = st.checkbox("Strict Seed Pairing",
                help="Require strict complementarity in seed region")
        else:
            max_mismatches = st.slider("Max Seed Mismatches", 0, 5, 2)
        if method == "psRNATarget (API)":
            expectation = st.slider("Expectation Cutoff", 1.0, 10.0, 5.0, 0.5)

    # Run prediction
    st.divider()

    if st.button("🔮 Predict Targets", type="primary", width="stretch"):
        if not mirna_seqs:
            st.error("No miRNA sequences available. Please provide sequences in FASTA format.")
            return

        with st.spinner("Running target prediction..."):
            if method == "psRNATarget (API)" and kingdom == "Plants":
                # Submit to psRNATarget
                st.info("Submitting job to psRNATarget server...")

                result = submit_psrnatarget_job(
                    mirna_sequences=mirna_seqs,
                    species=species,
                    expectation=expectation,
                    seed_region=(seed_start, seed_end),
                    custom_transcripts=st.session_state.get('custom_transcripts')
                )

                if result['status'] == 'submitted':
                    job_id = result['job_id']
                    st.session_state.psrnatarget_job_id = job_id
                    st.info(f"Job submitted! ID: {job_id}")

                    # Poll for results
                    progress_bar = st.progress(0)
                    status_text = st.empty()

                    max_wait = 300  # 5 minutes
                    waited = 0

                    while waited < max_wait:
                        status = check_psrnatarget_job(job_id)

                        if status['status'] == 'completed':
                            results = get_psrnatarget_results(job_id)

                            if results['status'] == 'success':
                                st.session_state.target_results = results['results']
                                st.success(f"Found {len(results['results'])} target predictions!")
                                break
                            else:
                                st.error(f"Failed to get results: {results.get('error')}")
                                break

                        elif status['status'] == 'error':
                            st.error(f"Job failed: {status.get('error')}")
                            break

                        progress = status.get('progress', 0) / 100
                        progress_bar.progress(progress)
                        status_text.text(f"Job status: {status['status']} ({status.get('progress', 0)}%)")

                        time.sleep(5)
                        waited += 5

                    if waited >= max_wait:
                        st.warning("Job is taking too long. Check back later with job ID.")

                else:
                    st.error(f"Submission failed: {result.get('error')}")

            elif method == "miRanda (Recommended)":
                # miRanda prediction for animals
                custom_transcripts = st.session_state.get('custom_transcripts')

                if not custom_transcripts:
                    st.error("miRanda requires 3' UTR sequences. Please upload a FASTA file with target UTR sequences.")
                    return

                # Format miRNA sequences as FASTA
                mirna_fasta = ""
                for mirna_id, seq in mirna_seqs.items():
                    mirna_fasta += f">{mirna_id}\n{seq}\n"

                st.info(f"Running miRanda: {len(mirna_seqs)} miRNAs against target UTRs...")

                # Run miRanda
                result = run_miranda(
                    mirna_fasta=mirna_fasta,
                    target_fasta=custom_transcripts,
                    score_threshold=score_threshold,
                    energy_threshold=energy_threshold,
                    strict=miranda_strict
                )

                if result['status'] == 'success':
                    results_df = result['results']
                    stats = result['stats']

                    st.session_state.target_results = results_df
                    st.session_state.miranda_stats = stats

                    st.success(f"✅ Found {stats['total_predictions']} target predictions!")
                    st.info(f"📊 {stats['unique_mirnas']} miRNAs targeting {stats['unique_targets']} genes")
                else:
                    st.error(f"miRanda failed: {result.get('error')}")
                    if 'install_instructions' in result:
                        with st.expander("📦 Installation Help"):
                            st.markdown(result['install_instructions'])

            else:
                # Local seed matching prediction
                custom_transcripts = st.session_state.get('custom_transcripts')

                if not custom_transcripts:
                    st.error("Local prediction requires custom transcript sequences. Please upload a FASTA file.")
                    return

                # Parse transcripts
                transcripts = {}
                current_id = None
                current_seq = ""
                for line in custom_transcripts.split('\n'):
                    line = line.strip()
                    if line.startswith('>'):
                        if current_id:
                            transcripts[current_id] = current_seq
                        current_id = line[1:].split()[0]
                        current_seq = ""
                    else:
                        current_seq += line
                if current_id:
                    transcripts[current_id] = current_seq

                st.info(f"Searching {len(mirna_seqs)} miRNAs against {len(transcripts)} transcripts...")

                results_df = predict_targets_local(
                    mirna_sequences=mirna_seqs,
                    transcript_sequences=transcripts,
                    seed_start=seed_start,
                    seed_end=seed_end,
                    max_mismatches=max_mismatches
                )

                st.session_state.target_results = results_df
                st.success(f"Found {len(results_df)} target predictions!")


def render_target_results():
    """Render target prediction results"""
    st.subheader("Target Prediction Results")

    results = st.session_state.get('target_results')

    if results is None:
        st.info("No results yet. Run target prediction first.")
        return

    if isinstance(results, pd.DataFrame) and len(results) == 0:
        st.warning("No targets found with current parameters.")
        return

    # Summary
    col1, col2, col3 = st.columns(3)

    with col1:
        st.metric("Total Predictions", len(results))

    with col2:
        n_mirnas = results['miRNA'].nunique() if 'miRNA' in results.columns else results.iloc[:, 0].nunique()
        st.metric("miRNAs with Targets", n_mirnas)

    with col3:
        target_col = 'Target' if 'Target' in results.columns else results.columns[1]
        n_targets = results[target_col].nunique()
        st.metric("Unique Targets", n_targets)

    # Filter options
    st.subheader("Filter Results")

    col1, col2 = st.columns(2)

    with col1:
        if 'miRNA' in results.columns:
            mirnas = ['All'] + sorted(results['miRNA'].unique().tolist())
            selected_mirna = st.selectbox("Filter by miRNA", mirnas)
        else:
            selected_mirna = 'All'

    with col2:
        if 'Seed_Mismatches' in results.columns:
            max_mm = st.slider("Max Seed Mismatches", 0, 5, 3)
        elif 'Expectation' in results.columns:
            max_exp = st.slider("Max Expectation", 0.0, 10.0, 5.0, 0.5)

    # Apply filters
    filtered = results.copy()

    if selected_mirna != 'All' and 'miRNA' in filtered.columns:
        filtered = filtered[filtered['miRNA'] == selected_mirna]

    if 'Seed_Mismatches' in filtered.columns:
        filtered = filtered[filtered['Seed_Mismatches'] <= max_mm]
    elif 'Expectation' in filtered.columns:
        filtered = filtered[filtered['Expectation'] <= max_exp]

    st.info(f"Showing {len(filtered)} predictions after filtering")

    # Results table
    st.dataframe(filtered, width="stretch", height=400)

    # Visualization
    st.subheader("Visualization")

    import plotly.express as px

    col1, col2 = st.columns(2)

    with col1:
        # Targets per miRNA
        if 'miRNA' in filtered.columns:
            targets_per_mirna = filtered['miRNA'].value_counts().head(20)

            fig = px.bar(
                x=targets_per_mirna.values,
                y=targets_per_mirna.index,
                orientation='h',
                title="Top 20 miRNAs by Number of Targets",
                labels={'x': 'Number of Targets', 'y': 'miRNA'}
            )
            fig.update_layout(yaxis={'categoryorder': 'total ascending'})
            st.plotly_chart(fig, width="stretch")

    with col2:
        # Score distribution
        if 'Complementarity' in filtered.columns:
            fig = px.histogram(
                filtered,
                x='Complementarity',
                nbins=30,
                title="Complementarity Score Distribution"
            )
            st.plotly_chart(fig, width="stretch")
        elif 'Expectation' in filtered.columns:
            fig = px.histogram(
                filtered,
                x='Expectation',
                nbins=30,
                title="Expectation Score Distribution"
            )
            st.plotly_chart(fig, width="stretch")

    # Download
    st.divider()

    col1, col2 = st.columns(2)

    with col1:
        csv = filtered.to_csv(index=False)
        st.download_button(
            "📥 Download Results (CSV)",
            csv,
            "target_predictions.csv",
            "text/csv"
        )

    with col2:
        # Get unique targets for enrichment
        target_col = 'Target' if 'Target' in filtered.columns else filtered.columns[1]
        unique_targets = filtered[target_col].unique().tolist()

        targets_text = '\n'.join(unique_targets)
        st.download_button(
            "📥 Download Target List",
            targets_text,
            "target_genes.txt",
            "text/plain"
        )

    # Link to enrichment
    st.divider()
    if st.button("➡️ Send Targets to Enrichment Analysis"):
        target_col = 'Target' if 'Target' in filtered.columns else filtered.columns[1]
        st.session_state.enrichment_genes = filtered[target_col].unique().tolist()
        st.success(f"Sent {len(st.session_state.enrichment_genes)} unique targets to enrichment analysis!")
