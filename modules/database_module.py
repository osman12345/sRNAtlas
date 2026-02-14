"""
Reference Database Builder Module for sRNAtlas
Download and index reference databases from miRBase, RNAcentral, etc.
"""
import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import subprocess
import tempfile
import shutil
import gzip
import requests
import io
import time

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from config import config
from utils.file_handlers import validate_safe_path, sanitize_filename


# Alternative miRBase URLs (GitHub mirror and FTP)
MIRBASE_URLS = [
    "https://mirbase.org/download/mature.fa",
    "https://mirbase.org/download/hairpin.fa",
    "ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz",
    "ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz",
]

# RNAcentral API
RNACENTRAL_API = "https://rnacentral.org/api/v1"

# Available organisms in miRBase (common ones)
MIRBASE_ORGANISMS = {
    "Arabidopsis thaliana": "ath",
    "Medicago truncatula": "mtr",
    "Oryza sativa": "osa",
    "Zea mays": "zma",
    "Glycine max": "gma",
    "Solanum lycopersicum": "sly",
    "Vitis vinifera": "vvi",
    "Homo sapiens": "hsa",
    "Mus musculus": "mmu",
    "Rattus norvegicus": "rno",
    "Drosophila melanogaster": "dme",
    "Caenorhabditis elegans": "cel",
    "Danio rerio": "dre",
}

# RNAcentral species with taxon IDs
RNACENTRAL_SPECIES = {
    "Arabidopsis thaliana": 3702,
    "Medicago truncatula": 3880,
    "Oryza sativa": 4530,
    "Zea mays": 4577,
    "Glycine max": 3847,
    "Homo sapiens": 9606,
    "Mus musculus": 10090,
    "Drosophila melanogaster": 7227,
    "Caenorhabditis elegans": 6239,
    "Danio rerio": 7955,
}

# RNA type mapping for RNAcentral
RNACENTRAL_RNA_TYPES = {
    "miRNA": "miRNA",
    "tRNA": "tRNA",
    "rRNA": "rRNA",
    "snoRNA": "snoRNA",
    "snRNA": "snRNA",
    "lncRNA": "lncRNA",
    "piRNA": "piRNA",
}

# Direct FTP download URLs for RNAcentral (FAST)
RNACENTRAL_FTP_BASE = "https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/by-database"

# Pre-built species FASTA files from RNAcentral FTP
RNACENTRAL_SPECIES_FILES = {
    "Arabidopsis thaliana": "https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_species_specific_ids.fasta.gz",
    "Homo sapiens": "https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_species_specific_ids.fasta.gz",
}


def check_bowtie() -> Tuple[bool, str]:
    """Check if bowtie is installed"""
    try:
        result = subprocess.run(
            ['bowtie', '--version'],
            capture_output=True,
            text=True
        )
        if result.returncode == 0:
            version = result.stdout.split('\n')[0]
            return True, version
        return False, "Not found"
    except FileNotFoundError:
        return False, "Not installed"


def download_with_retry(url: str, timeout: int = 120, retries: int = 3) -> Optional[requests.Response]:
    """Download URL with retry logic"""
    headers = {
        'User-Agent': 'Mozilla/5.0 (compatible; sRNA-WebTool/1.3)'
    }

    for attempt in range(retries):
        try:
            response = requests.get(url, timeout=timeout, headers=headers, stream=True)
            response.raise_for_status()
            return response
        except requests.exceptions.RequestException as e:
            if attempt < retries - 1:
                time.sleep(2 ** attempt)  # Exponential backoff
                continue
            raise e
    return None


def fetch_mirbase_sequences(species_code: str, seq_type: str = "mature") -> Dict:
    """
    Fetch miRNA sequences from miRBase using their API/direct download

    Args:
        species_code: 3-letter species code (e.g., 'ath', 'hsa')
        seq_type: 'mature' or 'hairpin'

    Returns:
        Dictionary with sequences
    """
    # Try multiple URL formats
    urls_to_try = [
        f"https://mirbase.org/download/{seq_type}.fa",
        f"https://www.mirbase.org/download/{seq_type}.fa",
        f"https://mirbase.org/ftp/CURRENT/{seq_type}.fa.gz",
    ]

    sequences = {}
    failed_urls = []

    for url in urls_to_try:
        try:
            response = download_with_retry(url, timeout=120)
            if response is None:
                continue

            # Handle gzipped content
            if url.endswith('.gz'):
                content = gzip.decompress(response.content).decode('utf-8')
            else:
                content = response.text

            # Parse FASTA and filter by species
            current_id = None
            current_seq = ""

            for line in content.split('\n'):
                line = line.strip()
                if line.startswith('>'):
                    if current_id and current_seq:
                        # Check if this is our species
                        if current_id.lower().startswith(f"{species_code.lower()}-"):
                            sequences[current_id] = current_seq.replace('U', 'T').replace('u', 't')

                    # Parse header: >cel-let-7-5p MIMAT0000001
                    parts = line[1:].split()
                    current_id = parts[0]
                    current_seq = ""
                else:
                    current_seq += line

            # Don't forget last sequence
            if current_id and current_seq:
                if current_id.lower().startswith(f"{species_code.lower()}-"):
                    sequences[current_id] = current_seq.replace('U', 'T').replace('u', 't')

            if sequences:
                return {
                    'status': 'success',
                    'sequences': sequences,
                    'count': len(sequences),
                    'source': url
                }

        except requests.exceptions.Timeout:
            failed_urls.append((url, "Connection timed out"))
            continue
        except requests.exceptions.ConnectionError:
            failed_urls.append((url, "Connection refused or network error"))
            continue
        except Exception as e:
            failed_urls.append((url, str(e)))
            continue

    error_details = '; '.join([f'{url}: {err}' for url, err in failed_urls]) if failed_urls else 'Unknown error'
    return {
        'status': 'error',
        'error': f'Could not download from any miRBase source. Tried {len(urls_to_try)} URLs.',
        'details': error_details,
        'urls_tried': [u for u, _ in failed_urls],
        'recommendation': 'Check your internet connection or try again later. Visit https://mirbase.org to download manually.'
    }


def fetch_rnacentral_ftp(
    taxon_id: int,
    rna_types: List[str] = None,
    max_sequences: int = 10000,
    progress_callback=None
) -> Dict:
    """
    FAST: Download sequences from RNAcentral FTP/direct files.
    This bypasses the slow API and downloads pre-built files.

    Args:
        taxon_id: NCBI taxonomy ID
        rna_types: List of RNA types to filter (applied after download)
        max_sequences: Maximum sequences to return
        progress_callback: Optional callback

    Returns:
        Dictionary with sequences
    """
    if progress_callback:
        progress_callback(0.1, "Downloading from RNAcentral (direct download)...")

    sequences = {}

    # Try ENA/EBI direct sequence fetch - much faster
    # Use the EBI proteins API which is more reliable
    ena_url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/textsearch"

    # Build query for ncRNA from specific taxon
    rna_filter = ""
    if rna_types:
        rna_terms = " OR ".join([f'"{rt}"' for rt in rna_types])
        rna_filter = f" AND ({rna_terms})"

    query = f'tax_tree({taxon_id}) AND molecule_type="ncRNA"{rna_filter}'

    try:
        if progress_callback:
            progress_callback(0.2, "Querying ENA database...")

        # ENA text search for ncRNA
        search_url = "https://www.ebi.ac.uk/ena/browser/api/embl/textsearch"
        params = {
            'query': f'tax_tree({taxon_id})',
            'result': 'sequence',
            'limit': min(max_sequences, 5000),
            'format': 'fasta',
            'dataclass': 'STD'
        }

        # Alternative: use Ensembl BioMart for plant species
        if taxon_id in [3702, 4530, 4577, 3847]:  # Plants
            return fetch_ensembl_ncrna(taxon_id, rna_types, max_sequences, progress_callback)

        # For other species, try RFam database
        return fetch_rfam_sequences(taxon_id, rna_types, max_sequences, progress_callback)

    except Exception as e:
        if progress_callback:
            progress_callback(None, f"Direct download failed: {str(e)[:50]}")

    return {
        'status': 'error',
        'error': 'Direct download failed. Try the Standard API method.'
    }


def fetch_ensembl_ncrna(
    taxon_id: int,
    rna_types: List[str] = None,
    max_sequences: int = 5000,
    progress_callback=None
) -> Dict:
    """
    Fetch ncRNA from Ensembl - Note: This has limitations.
    For miRNA, miRBase is recommended instead.
    """
    # Ensembl REST API is limited for bulk ncRNA downloads
    # Return error to trigger fallback to RNAcentral API
    return {'status': 'error', 'error': 'Ensembl method limited. Using RNAcentral API instead.'}


def fetch_rfam_sequences(
    taxon_id: int,
    rna_types: List[str] = None,
    max_sequences: int = 2000,
    progress_callback=None
) -> Dict:
    """
    Fetch sequences from Rfam database - good for ncRNA families
    """
    if progress_callback:
        progress_callback(0.2, "Fetching from Rfam database...")

    sequences = {}

    # Rfam family IDs for common ncRNA types
    rfam_families = {
        "miRNA": ["RF00001", "RF00002", "RF00003"],  # Example families
        "tRNA": ["RF00005", "RF01852"],
        "snoRNA": ["RF00012", "RF00013"],
        "rRNA": ["RF00001", "RF00002"],
    }

    try:
        # Use Rfam API to get sequences for taxon
        rfam_url = f"https://rfam.org/taxonomy/{taxon_id}/sequences"

        response = requests.get(rfam_url, timeout=60, headers={'Accept': 'text/plain'})

        if response.status_code == 200:
            # Parse FASTA response
            content = response.text
            current_id = None
            current_seq = ""

            for line in content.split('\n'):
                line = line.strip()
                if line.startswith('>'):
                    if current_id and current_seq:
                        sequences[current_id] = current_seq.replace('U', 'T').replace('u', 't')
                        if len(sequences) >= max_sequences:
                            break
                    current_id = line[1:].split()[0][:50]
                    current_seq = ""
                else:
                    current_seq += line

            if current_id and current_seq and len(sequences) < max_sequences:
                sequences[current_id] = current_seq.replace('U', 'T').replace('u', 't')

            if progress_callback:
                progress_callback(1.0, f"Downloaded {len(sequences)} sequences from Rfam")

    except Exception as e:
        if progress_callback:
            progress_callback(None, f"Rfam fetch failed: {str(e)[:40]}")

    if sequences:
        return {'status': 'success', 'sequences': sequences, 'count': len(sequences)}

    # Fallback to API method
    return {'status': 'error', 'error': 'Rfam download failed'}


def fetch_rnacentral_sequences(
    taxon_id: int,
    rna_types: List[str] = None,
    max_sequences: int = 5000,
    progress_callback=None
) -> Dict:
    """
    Fetch ncRNA sequences from RNAcentral API - FAST batch method.
    Uses the /rna endpoint which returns sequences directly in the response.

    Args:
        taxon_id: NCBI taxonomy ID
        rna_types: List of RNA types to fetch
        max_sequences: Maximum number of sequences to fetch
        progress_callback: Optional progress callback

    Returns:
        Dictionary with sequences
    """
    if rna_types is None:
        rna_types = ["miRNA"]

    sequences = {}
    total_fetched = 0

    # Direct RNAcentral API - sequences included in response
    base_url = "https://rnacentral.org/api/v1/rna"

    for rna_type in rna_types:
        if total_fetched >= max_sequences:
            break

        if progress_callback:
            progress_callback(
                0.05 + (total_fetched / max_sequences) * 0.9,
                f"Fetching {rna_type} sequences..."
            )

        page = 1
        consecutive_failures = 0
        page_size = 100  # Max allowed by API

        while total_fetched < max_sequences:
            params = {
                'taxid': taxon_id,
                'rna_type': rna_type.lower(),
                'format': 'json',
                'page_size': page_size,
                'page': page
            }

            try:
                response = requests.get(base_url, params=params, timeout=60)

                if response.status_code != 200:
                    consecutive_failures += 1
                    if consecutive_failures >= 3:
                        break
                    time.sleep(2)
                    continue

                consecutive_failures = 0
                data = response.json()
                results = data.get('results', [])

                if not results:
                    break

                for entry in results:
                    if total_fetched >= max_sequences:
                        break

                    urs_id = entry.get('upi', '')
                    sequence = entry.get('sequence', '')
                    description = entry.get('description', '')

                    if urs_id and sequence:
                        clean_desc = description.replace(' ', '_').replace('/', '_')[:40] if description else rna_type
                        seq_id = f"{clean_desc}_{urs_id}_{taxon_id}"
                        sequences[seq_id] = sequence.replace('U', 'T').replace('u', 't')
                        total_fetched += 1

                if progress_callback:
                    progress_callback(
                        0.05 + (total_fetched / max_sequences) * 0.9,
                        f"Downloaded {total_fetched} {rna_type} sequences..."
                    )

                # Check if more pages
                if not data.get('next'):
                    break

                page += 1
                time.sleep(0.3)  # Rate limiting

            except requests.exceptions.Timeout:
                consecutive_failures += 1
                if consecutive_failures >= 3:
                    if progress_callback:
                        progress_callback(None, f"Timeout after {total_fetched} sequences")
                    break
                time.sleep(3)
                continue
            except requests.exceptions.ConnectionError as e:
                consecutive_failures += 1
                if consecutive_failures >= 3:
                    if progress_callback:
                        progress_callback(None, f"Network error after {total_fetched} sequences: {str(e)[:50]}")
                    break
                time.sleep(3)
                continue
            except json.JSONDecodeError as e:
                consecutive_failures += 1
                if consecutive_failures >= 3:
                    if progress_callback:
                        progress_callback(None, f"Invalid API response after {total_fetched} sequences")
                    break
                continue
            except Exception as e:
                consecutive_failures += 1
                if consecutive_failures >= 3:
                    if progress_callback:
                        progress_callback(None, f"Error after {total_fetched} sequences: {str(e)[:50]}")
                    break
                continue

    if progress_callback:
        progress_callback(1.0, f"Complete! Downloaded {len(sequences)} sequences")

    if sequences:
        return {
            'status': 'success',
            'sequences': sequences,
            'count': len(sequences)
        }
    else:
        return {
            'status': 'error',
            'error': 'No sequences found. Check species/RNA type or try again.'
        }


def fetch_rnacentral_via_text_search(
    taxon_id: int,
    rna_types: List[str] = None,
    max_sequences: int = 2000,
    progress_callback=None
) -> Dict:
    """
    Fetch sequences using RNAcentral text search API.
    This method uses a different endpoint that may work when the standard fails.

    Args:
        taxon_id: NCBI taxonomy ID
        rna_types: List of RNA types
        max_sequences: Max sequences to fetch
        progress_callback: Optional callback

    Returns:
        Dictionary with sequences
    """
    sequences = {}

    if rna_types is None:
        rna_types = ["miRNA"]

    if progress_callback:
        progress_callback(0.05, "Connecting to RNAcentral text search...")

    total_fetched = 0

    for rna_type in rna_types:
        if total_fetched >= max_sequences:
            break

        # Build text search query
        query = f'taxid:{taxon_id} AND rna_type:"{rna_type}"'
        search_url = "https://rnacentral.org/api/v1/rna"

        if progress_callback:
            progress_callback(
                0.1 + (total_fetched / max_sequences) * 0.85,
                f"Searching for {rna_type}..."
            )

        page = 1
        consecutive_failures = 0

        while total_fetched < max_sequences:
            params = {
                'query': query,
                'format': 'json',
                'page_size': 50,
                'page': page
            }

            try:
                response = requests.get(search_url, params=params, timeout=45)

                if response.status_code != 200:
                    consecutive_failures += 1
                    if consecutive_failures >= 3:
                        break
                    time.sleep(2)
                    continue

                consecutive_failures = 0
                data = response.json()
                results = data.get('results', [])

                if not results:
                    break

                for entry in results:
                    if total_fetched >= max_sequences:
                        break

                    urs_id = entry.get('upi', '')
                    sequence = entry.get('sequence', '')
                    description = entry.get('description', '')

                    if urs_id and sequence:
                        clean_desc = description.replace(' ', '_').replace('/', '_')[:40] if description else rna_type
                        seq_id = f"{clean_desc}_{urs_id}_{taxon_id}"
                        sequences[seq_id] = sequence.replace('U', 'T').replace('u', 't')
                        total_fetched += 1

                if progress_callback:
                    progress_callback(
                        0.1 + (total_fetched / max_sequences) * 0.85,
                        f"Downloaded {total_fetched} {rna_type} sequences..."
                    )

                if not data.get('next'):
                    break

                page += 1
                time.sleep(0.5)

            except requests.exceptions.Timeout:
                consecutive_failures += 1
                if consecutive_failures >= 3:
                    break
                time.sleep(3)
                continue
            except requests.exceptions.ConnectionError as e:
                consecutive_failures += 1
                if consecutive_failures >= 3:
                    if progress_callback:
                        progress_callback(None, f"Network error after {total_fetched} sequences: {str(e)[:50]}")
                    break
                time.sleep(3)
                continue
            except json.JSONDecodeError as e:
                consecutive_failures += 1
                if consecutive_failures >= 3:
                    if progress_callback:
                        progress_callback(None, f"Invalid API response after {total_fetched} sequences")
                    break
                continue
            except Exception as e:
                consecutive_failures += 1
                if consecutive_failures >= 3:
                    if progress_callback:
                        progress_callback(None, f"Error after {total_fetched} sequences: {str(e)[:50]}")
                    break
                continue

    if progress_callback:
        progress_callback(1.0, f"Complete! Downloaded {len(sequences)} sequences")

    if sequences:
        return {
            'status': 'success',
            'sequences': sequences,
            'count': len(sequences)
        }
    else:
        return {
            'status': 'error',
            'error': 'No sequences found. The API may be slow - try again later.'
        }


def save_sequences_to_fasta(sequences: Dict[str, str], output_file: Path) -> Dict:
    """Save sequences dictionary to FASTA file"""
    try:
        output_file = validate_safe_path(output_file)
    except ValueError as e:
        return {'status': 'error', 'error': f'Invalid output path: {e}'}

    output_file.parent.mkdir(parents=True, exist_ok=True)

    try:
        # Validate sequences contain only valid nucleotide characters
        valid_chars = set('ATCGNatcgnURYSWKMBDHVuryswkmbdhv')
        for seq_id, sequence in sequences.items():
            invalid = set(sequence) - valid_chars
            if invalid:
                return {
                    'status': 'error',
                    'error': f'Sequence {seq_id} contains invalid characters: {invalid}'
                }

        with open(output_file, 'w') as f:
            for seq_id, sequence in sequences.items():
                f.write(f">{seq_id}\n{sequence}\n")

        return {
            'status': 'success',
            'file': output_file,
            'count': len(sequences)
        }
    except Exception as e:
        return {
            'status': 'error',
            'error': str(e)
        }


def build_bowtie_index(
    fasta_file: Path,
    index_prefix: Path,
    threads: int = 4,
    progress_callback=None
) -> Dict:
    """
    Build Bowtie index from FASTA file
    """
    try:
        fasta_file = validate_safe_path(fasta_file, must_exist=True)
        index_prefix = validate_safe_path(index_prefix.parent) / index_prefix.name
    except ValueError as e:
        return {'status': 'error', 'error': f'Invalid path: {e}'}

    try:
        if progress_callback:
            progress_callback(0.1, "Building Bowtie index...")

        cmd = [
            'bowtie-build',
            '-f',
            str(fasta_file),
            str(index_prefix)
        ]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=1800
        )

        if result.returncode != 0:
            return {
                'status': 'error',
                'error': result.stderr
            }

        if progress_callback:
            progress_callback(1.0, "Index built successfully!")

        index_files = list(index_prefix.parent.glob(f"{index_prefix.name}*.ebwt"))

        return {
            'status': 'success',
            'index_prefix': str(index_prefix),
            'index_files': [str(f) for f in index_files],
            'message': f"Built index with {len(index_files)} files"
        }

    except subprocess.TimeoutExpired:
        return {
            'status': 'error',
            'error': 'Index building timed out'
        }
    except FileNotFoundError:
        return {
            'status': 'error',
            'error': 'bowtie-build not found. Please install Bowtie.'
        }
    except Exception as e:
        return {
            'status': 'error',
            'error': str(e)
        }


def merge_fasta_files(
    input_files: List[Path],
    output_file: Path,
    deduplicate: bool = True
) -> Dict:
    """Merge multiple FASTA files into one"""
    sequences = {}
    total_input = 0

    for fasta_file in input_files:
        with open(fasta_file, 'r') as f:
            current_id = None
            current_seq = ""

            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id and current_seq:
                        total_input += 1
                        if deduplicate:
                            sequences[current_seq] = current_id
                        else:
                            sequences[f"{current_id}_{total_input}"] = current_seq
                    current_id = line[1:]
                    current_seq = ""
                else:
                    current_seq += line

            if current_id and current_seq:
                total_input += 1
                if deduplicate:
                    sequences[current_seq] = current_id
                else:
                    sequences[f"{current_id}_{total_input}"] = current_seq

    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, 'w') as out:
        if deduplicate:
            for seq, seq_id in sequences.items():
                out.write(f">{seq_id}\n{seq}\n")
        else:
            for seq_id, seq in sequences.items():
                out.write(f">{seq_id}\n{seq}\n")

    return {
        'status': 'success',
        'input_sequences': total_input,
        'output_sequences': len(sequences),
        'file': output_file
    }


def render_database_page():
    """Render the database builder page"""
    st.header("🗄️ Reference Databases")

    st.markdown("""
    Select or download reference databases for small RNA alignment.
    Pre-loaded databases are ready to use - just select and build the index.
    """)

    # Check bowtie
    bowtie_ok, version = check_bowtie()

    if bowtie_ok:
        st.success(f"✅ Bowtie installed: {version}")
    else:
        st.warning("""
        ⚠️ Bowtie not found. You can still view references, but need Bowtie to build indices.

        Install with: `conda install -c bioconda bowtie` or `brew install bowtie`
        """)

    tab1, tab2, tab3, tab4 = st.tabs([
        "📚 Pre-loaded References",
        "📁 Custom Upload",
        "🌐 Download Online",
        "🔧 Manage"
    ])

    with tab1:
        render_preloaded_references(bowtie_ok)

    with tab2:
        render_custom_database(bowtie_ok)

    with tab3:
        render_online_download(bowtie_ok)

    with tab4:
        render_database_management(bowtie_ok)


def render_preloaded_references(bowtie_ok: bool):
    """Render pre-loaded reference database selection"""
    st.subheader("Pre-loaded Reference Databases")

    st.markdown("""
    Select from bundled reference databases. These include sequences from miRBase
    and other curated sources for common model organisms.

    **Supported RNA types:** miRNA, tRNA, rRNA, snoRNA, snRNA, piRNA, siRNA
    """)

    # Load references.json
    ref_dir = Path(config.paths.reference_dir)
    json_file = ref_dir / "references.json"

    if not json_file.exists():
        st.warning("No reference database index found. Use Custom Upload or Download Online tabs.")
        return

    import json
    try:
        with open(json_file, 'r') as f:
            ref_data = json.load(f)
    except Exception as e:
        st.error(f"Error loading references: {e}")
        return

    references = ref_data.get('references', [])

    if not references:
        st.info("No pre-loaded references available. Download from the Online tab or upload custom files.")
        return

    # Group by species
    species_refs = {}
    for ref in references:
        species = ref.get('species', 'Unknown')
        if species not in species_refs:
            species_refs[species] = []
        species_refs[species].append(ref)

    # Species filter
    all_species = sorted(species_refs.keys())
    selected_species = st.selectbox(
        "Filter by Species",
        options=["All"] + all_species
    )

    # Display references
    st.markdown("### Available References")

    display_refs = references if selected_species == "All" else species_refs.get(selected_species, [])

    for ref in display_refs:
        ref_id = ref.get('id', 'unknown')
        ref_name = ref.get('name', ref_id)
        ref_file = ref.get('file', '')
        ref_source = ref.get('source', 'Unknown')
        ref_rna_types = ref.get('rna_types', [])
        ref_seqs = ref.get('sequences', 0)
        ref_recommended = ref.get('recommended', False)

        # Check if file exists
        file_path = ref_dir / ref_file
        file_exists = file_path.exists()

        # Check if index exists
        if file_exists:
            index_prefix = file_path.with_suffix('')
            has_index = any(index_prefix.parent.glob(f"{index_prefix.name}*.ebwt"))
        else:
            has_index = False

        # Display card
        badge = "⭐ " if ref_recommended else ""
        status_icon = "✅" if has_index else ("📄" if file_exists else "❌")

        with st.expander(f"{status_icon} {badge}{ref_name}", expanded=ref_recommended):
            col1, col2, col3 = st.columns([2, 2, 1])

            with col1:
                st.caption(f"**Source:** {ref_source}")
                st.caption(f"**RNA Types:** {', '.join(ref_rna_types)}")

            with col2:
                st.caption(f"**Sequences:** {ref_seqs:,}")
                if file_exists:
                    size_mb = file_path.stat().st_size / (1024 * 1024)
                    st.caption(f"**Size:** {size_mb:.2f} MB")

            with col3:
                if has_index:
                    st.success("Ready")
                elif file_exists:
                    st.warning("Needs index")
                else:
                    st.error("Missing")

            # Actions
            if file_exists and not has_index and bowtie_ok:
                if st.button(f"🔨 Build Index", key=f"build_{ref_id}"):
                    with st.spinner(f"Building Bowtie index for {ref_name}..."):
                        index_prefix = file_path.with_suffix('')
                        result = build_bowtie_index(file_path, index_prefix)
                        if result['status'] == 'success':
                            st.success("✅ Index built successfully!")
                            st.rerun()
                        else:
                            st.error(f"❌ {result.get('error')}")

            if has_index:
                if st.button(f"✅ Select for Alignment", key=f"select_{ref_id}"):
                    st.session_state.active_reference = str(file_path.with_suffix(''))
                    st.session_state.selected_ref_name = ref_name
                    st.success(f"Selected **{ref_name}** for alignment")

    # Show current selection
    st.divider()
    active_ref = st.session_state.get('active_reference')
    ref_name = st.session_state.get('selected_ref_name', '')

    if active_ref:
        st.info(f"**Current Selection:** {ref_name or Path(active_ref).name}")
    else:
        st.warning("No reference selected. Build an index and click 'Select for Alignment'.")


def render_online_download(bowtie_ok: bool):
    """Render online download options (miRBase and RNAcentral)"""
    st.subheader("Download from Online Databases")

    st.markdown("""
    Download additional reference sequences from online databases.

    **Recommended approach:**
    - Use **miRBase** for miRNA sequences (fast and reliable)
    - Use **RNAcentral** for other ncRNA types (tRNA, snoRNA, etc.)
    """)

    download_source = st.radio(
        "Select Source",
        options=["miRBase (miRNA)", "RNAcentral (all ncRNA types)"],
        horizontal=True
    )

    if download_source == "miRBase (miRNA)":
        render_mirbase_download(bowtie_ok)
    else:
        render_rnacentral_download(bowtie_ok)


def render_mirbase_download(bowtie_ok: bool):
    """Render miRBase download section"""
    st.subheader("miRBase Download")

    st.markdown("""
    Download miRNA sequences from [miRBase](https://www.mirbase.org/).
    Sequences are fetched via API and filtered for your species.
    """)

    col1, col2 = st.columns(2)

    with col1:
        species = st.selectbox(
            "Species",
            options=list(MIRBASE_ORGANISMS.keys()),
            key="mirbase_species"
        )
        species_code = MIRBASE_ORGANISMS[species]
        st.caption(f"Species code: {species_code}")

    with col2:
        seq_type = st.radio(
            "Sequence Type",
            options=["Mature", "Hairpin", "Both"],
            help="Mature sequences are ~22nt, hairpins are ~70-100nt"
        )

    # Output directory
    db_dir = Path(config.paths.reference_dir) / "mirbase"
    st.text_input("Output Directory", value=str(db_dir), disabled=True)

    # Download button
    if st.button("📥 Download from miRBase", type="primary", key="dl_mirbase"):
        db_dir.mkdir(parents=True, exist_ok=True)

        progress_bar = st.progress(0)
        status_text = st.empty()

        downloaded_files = []

        if seq_type in ["Mature", "Both"]:
            status_text.text("Fetching mature miRNA sequences...")
            progress_bar.progress(0.1)

            result = fetch_mirbase_sequences(species_code, "mature")

            if result['status'] == 'success':
                output_file = db_dir / f"mirbase_{species_code}_mature.fa"
                save_result = save_sequences_to_fasta(result['sequences'], output_file)

                if save_result['status'] == 'success':
                    st.success(f"✅ Downloaded {result['count']} mature miRNA sequences")
                    downloaded_files.append(output_file)
                    progress_bar.progress(0.5)
                else:
                    st.error(f"❌ Failed to save: {save_result.get('error')}")
            else:
                st.error(f"❌ {result.get('error', 'Download failed')}")

        if seq_type in ["Hairpin", "Both"]:
            status_text.text("Fetching hairpin sequences...")
            progress_bar.progress(0.6)

            result = fetch_mirbase_sequences(species_code, "hairpin")

            if result['status'] == 'success':
                output_file = db_dir / f"mirbase_{species_code}_hairpin.fa"
                save_result = save_sequences_to_fasta(result['sequences'], output_file)

                if save_result['status'] == 'success':
                    st.success(f"✅ Downloaded {result['count']} hairpin sequences")
                    downloaded_files.append(output_file)
                    progress_bar.progress(0.9)
                else:
                    st.error(f"❌ Failed to save: {save_result.get('error')}")
            else:
                st.error(f"❌ {result.get('error', 'Download failed')}")

        progress_bar.progress(1.0)
        status_text.text("Download complete!")

        # Build index option
        if downloaded_files and bowtie_ok:
            st.divider()
            st.markdown("### Build Bowtie Index")

            for fasta_file in downloaded_files:
                col1, col2 = st.columns([3, 1])
                with col1:
                    st.text(f"📄 {fasta_file.name}")
                with col2:
                    if st.button("Build Index", key=f"idx_{fasta_file.name}"):
                        index_prefix = fasta_file.with_suffix('')
                        with st.spinner(f"Building index..."):
                            idx_result = build_bowtie_index(fasta_file, index_prefix)
                            if idx_result['status'] == 'success':
                                st.success(f"✅ Index built!")
                            else:
                                st.error(f"❌ {idx_result.get('error')}")


def render_rnacentral_download(bowtie_ok: bool):
    """Render RNAcentral download section"""
    st.subheader("RNAcentral / ncRNA Download")

    st.markdown("""
    Download ncRNA sequences from RNAcentral database.

    💡 **For miRNA:** Use the **miRBase** tab instead - it's faster and more complete.

    **RNAcentral** is useful for other ncRNA types (tRNA, snoRNA, rRNA, piRNA, etc.)
    """)

    col1, col2 = st.columns(2)

    with col1:
        species = st.selectbox(
            "Species",
            options=list(RNACENTRAL_SPECIES.keys()),
            key="rnacentral_species"
        )
        taxon_id = RNACENTRAL_SPECIES[species]
        st.caption(f"Taxon ID: {taxon_id}")

    with col2:
        max_seqs = st.slider(
            "Max Sequences",
            min_value=100,
            max_value=5000,
            value=1000,
            step=100,
            help="Maximum sequences to download"
        )

    rna_types = st.multiselect(
        "RNA Types",
        options=list(RNACENTRAL_RNA_TYPES.keys()),
        default=["miRNA"],
        help="Select RNA types to download"
    )

    # Note about speed
    st.info("⏱️ RNAcentral API can be slow (1-5 minutes for 500+ sequences). For miRNA, use miRBase tab instead.")

    # Output directory
    db_dir = Path(config.paths.reference_dir) / "rnacentral"
    st.text_input("Output Directory", value=str(db_dir), disabled=True, key="rnacentral_dir")

    if st.button("📥 Download ncRNA Sequences", type="primary", key="dl_rnacentral"):
        if not rna_types:
            st.error("Please select at least one RNA type")
            return

        db_dir.mkdir(parents=True, exist_ok=True)

        progress_bar = st.progress(0)
        status_text = st.empty()

        def update_progress(progress, message):
            if progress is not None:
                progress_bar.progress(min(progress, 1.0))
            status_text.text(message)

        update_progress(0.05, "Connecting to RNAcentral API...")

        # Use RNAcentral API directly
        result = fetch_rnacentral_sequences(
            taxon_id=taxon_id,
            rna_types=rna_types,
            max_sequences=max_seqs,
            progress_callback=update_progress
        )

        if result['status'] == 'success' and result['count'] > 0:
            # Generate filename
            species_short = species.replace(' ', '_')
            rna_str = '_'.join(rna_types[:3])
            output_file = db_dir / f"rnacentral_{species_short}_{rna_str}.fa"

            save_result = save_sequences_to_fasta(result['sequences'], output_file)

            if save_result['status'] == 'success':
                st.success(f"✅ Downloaded {result['count']} sequences to {output_file.name}")

                # Offer to build index
                if bowtie_ok:
                    st.divider()
                    if st.button("🔨 Build Bowtie Index", key="build_rnacentral_idx"):
                        index_prefix = output_file.with_suffix('')
                        with st.spinner("Building index..."):
                            idx_result = build_bowtie_index(output_file, index_prefix)
                            if idx_result['status'] == 'success':
                                st.success("✅ Index built successfully!")
                            else:
                                st.error(f"❌ {idx_result.get('error')}")
            else:
                st.error(f"❌ Failed to save: {save_result.get('error')}")
        else:
            st.error(f"❌ {result.get('error', 'No sequences found.')}")
            st.warning("💡 **Tips:** Try reducing max sequences, selecting fewer RNA types, or use the Text Search method.")

            # Provide manual download link as fallback
            st.info(f"""
            **Alternative:** Download manually from RNAcentral:

            1. Visit [RNAcentral Search](https://rnacentral.org/search?q=taxonomy:{taxon_id})
            2. Filter by RNA types
            3. Click "Download" → "FASTA"
            4. Upload the file in the "Custom" tab
            """)


def render_custom_database(bowtie_ok: bool):
    """Render custom database upload section"""
    st.subheader("Custom Reference Database")

    st.markdown("""
    Upload your own reference sequences in FASTA format.
    You can combine multiple files and build a unified index.
    """)

    uploaded_files = st.file_uploader(
        "Upload FASTA files",
        type=['fasta', 'fa', 'fna', 'txt'],
        accept_multiple_files=True,
        key="custom_fasta_upload"
    )

    if uploaded_files:
        st.info(f"Uploaded {len(uploaded_files)} files")

        # Show file info
        total_seqs = 0
        for uf in uploaded_files:
            content = uf.getvalue().decode('utf-8')
            seq_count = content.count('>')
            total_seqs += seq_count
            st.caption(f"📄 {uf.name}: {seq_count} sequences")

        st.write(f"**Total: {total_seqs} sequences**")

        # Options
        col1, col2 = st.columns(2)

        with col1:
            db_name = st.text_input(
                "Database Name",
                value="custom_reference",
                help="Name for the output database"
            )

        with col2:
            deduplicate = st.checkbox(
                "Remove duplicate sequences",
                value=True
            )

        if st.button("🔨 Build Custom Database", type="primary", key="build_custom"):
            db_dir = Path(config.paths.reference_dir) / "custom"
            db_dir.mkdir(parents=True, exist_ok=True)

            # Save uploaded files
            temp_files = []
            for uf in uploaded_files:
                temp_path = db_dir / f"temp_{uf.name}"
                with open(temp_path, 'wb') as f:
                    f.write(uf.getbuffer())
                temp_files.append(temp_path)

            # Merge files
            merged_file = db_dir / f"{db_name}.fa"

            with st.spinner("Merging FASTA files..."):
                result = merge_fasta_files(
                    input_files=temp_files,
                    output_file=merged_file,
                    deduplicate=deduplicate
                )

            # Clean up temp files
            for tf in temp_files:
                tf.unlink()

            if result['status'] == 'success':
                st.success(f"✅ Merged {result['input_sequences']} → {result['output_sequences']} sequences")

                # Build index
                if bowtie_ok:
                    index_prefix = merged_file.with_suffix('')

                    with st.spinner("Building Bowtie index..."):
                        idx_result = build_bowtie_index(
                            fasta_file=merged_file,
                            index_prefix=index_prefix
                        )

                        if idx_result['status'] == 'success':
                            st.success(f"✅ Index built: {index_prefix.name}")
                        else:
                            st.error(f"❌ Index failed: {idx_result.get('error')}")
            else:
                st.error(f"❌ Merge failed: {result.get('error')}")


def render_database_management(bowtie_ok: bool):
    """Render database management section"""
    st.subheader("Manage Databases")

    ref_dir = Path(config.paths.reference_dir)

    if not ref_dir.exists():
        st.info("No reference databases found. Download or upload sequences first.")
        return

    # Find all FASTA files and indices
    fasta_files = list(ref_dir.rglob("*.fa")) + list(ref_dir.rglob("*.fasta"))

    if not fasta_files:
        st.info("No databases found. Download or upload reference sequences first.")
        return

    # Group by database
    databases = {}

    for fa in fasta_files:
        db_name = fa.stem
        db_path = fa.parent

        # Check if index exists
        has_index = any(db_path.glob(f"{db_name}*.ebwt"))

        # Count sequences
        try:
            with open(fa, 'r') as f:
                seq_count = sum(1 for line in f if line.startswith('>'))
        except:
            seq_count = 0

        databases[db_name] = {
            'fasta': fa,
            'indexed': has_index,
            'sequences': seq_count,
            'size_mb': fa.stat().st_size / (1024 * 1024),
            'path': db_path
        }

    # Display databases
    st.markdown("### Available Databases")

    for db_name, info in databases.items():
        with st.expander(f"📁 {db_name}", expanded=False):
            col1, col2, col3 = st.columns(3)

            with col1:
                st.metric("Sequences", f"{info['sequences']:,}")

            with col2:
                st.metric("Size", f"{info['size_mb']:.2f} MB")

            with col3:
                if info['indexed']:
                    st.success("✅ Indexed")
                else:
                    st.warning("⚠️ Not indexed")

            st.text(f"Path: {info['fasta']}")

            col1, col2, col3 = st.columns(3)

            with col1:
                if not info['indexed'] and bowtie_ok:
                    if st.button("Build Index", key=f"idx_{db_name}"):
                        idx_prefix = info['fasta'].with_suffix('')
                        with st.spinner("Building index..."):
                            result = build_bowtie_index(info['fasta'], idx_prefix)
                            if result['status'] == 'success':
                                st.success("Index built!")
                                st.rerun()
                            else:
                                st.error(result.get('error'))

            with col2:
                if st.button("Set Active", key=f"active_{db_name}"):
                    st.session_state.active_reference = str(info['fasta'].with_suffix(''))
                    st.success(f"✅ Set {db_name} as active database")

            with col3:
                if st.button("🗑️ Delete", key=f"del_{db_name}"):
                    # Delete FASTA and index files
                    info['fasta'].unlink()
                    for idx in info['path'].glob(f"{db_name}*.ebwt"):
                        idx.unlink()
                    st.success(f"Deleted {db_name}")
                    st.rerun()

    # Show active database
    st.divider()
    active = st.session_state.get('active_reference')
    if active:
        st.info(f"**Active Database:** {Path(active).name}")
    else:
        st.warning("No active database set. Click 'Set Active' on a database above.")
