"""
Reference Manager for sRNAtlas
Handles loading, listing, and managing reference databases
"""
import json
from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd


class ReferenceManager:
    """Manages reference databases for sRNA analysis"""

    def __init__(self, reference_dir: Optional[Path] = None):
        """
        Initialize the reference manager

        Args:
            reference_dir: Path to references directory. If None, uses default location.
        """
        if reference_dir is None:
            # Default: data/references relative to project root
            self.reference_dir = Path(__file__).parent.parent / "data" / "references"
        else:
            self.reference_dir = Path(reference_dir)

        self.config_file = self.reference_dir / "references.json"
        self.references = {}
        self._load_config()

    def _load_config(self):
        """Load reference configuration from JSON file"""
        if self.config_file.exists():
            with open(self.config_file, 'r') as f:
                data = json.load(f)
                self.references = {ref['id']: ref for ref in data.get('references', [])}
                self.metadata = data.get('metadata', {})
        else:
            self.references = {}
            self.metadata = {}
            # Try to auto-detect references
            self._auto_detect_references()

    def _auto_detect_references(self):
        """Auto-detect FASTA files in reference directory and subdirectories"""
        if not self.reference_dir.exists():
            return

        # Scan main directory and subdirectories (mirbase, rnacentral, custom)
        search_dirs = [self.reference_dir]
        for subdir in ['mirbase', 'rnacentral', 'custom']:
            subdir_path = self.reference_dir / subdir
            if subdir_path.exists():
                search_dirs.append(subdir_path)

        for search_dir in search_dirs:
            # Determine source based on directory
            if 'mirbase' in str(search_dir):
                source = 'miRBase'
            elif 'rnacentral' in str(search_dir):
                source = 'RNAcentral'
            else:
                source = 'Custom'

            # Search for FASTA files
            for ext in ['*.fasta', '*.fa', '*.fna']:
                for fasta_file in search_dir.glob(ext):
                    ref_id = fasta_file.stem
                    # Use relative path from reference_dir
                    try:
                        rel_path = fasta_file.relative_to(self.reference_dir)
                    except ValueError:
                        rel_path = fasta_file.name

                    if ref_id not in self.references:
                        # Try to parse species and RNA type from filename
                        species, rna_types = self._parse_filename(fasta_file.stem)

                        # Check if bowtie index exists
                        has_index = (fasta_file.parent / f"{ref_id}.1.ebwt").exists()

                        self.references[ref_id] = {
                            'id': ref_id,
                            'name': ref_id.replace('_', ' ').title(),
                            'species': species,
                            'file': str(rel_path),
                            'full_path': str(fasta_file),
                            'description': f'Auto-detected from {source}: {fasta_file.name}',
                            'source': source,
                            'rna_types': rna_types,
                            'has_index': has_index,
                            'sequences': self._count_sequences(fasta_file)
                        }

    def _parse_filename(self, filename: str) -> tuple:
        """Parse species and RNA type from filename"""
        species = 'Unknown'
        rna_types = ['unknown']

        # Common species codes
        species_codes = {
            'ath': 'Arabidopsis thaliana',
            'mtr': 'Medicago truncatula',
            'osa': 'Oryza sativa',
            'zma': 'Zea mays',
            'gma': 'Glycine max',
            'sly': 'Solanum lycopersicum',
            'vvi': 'Vitis vinifera',
            'hsa': 'Homo sapiens',
            'mmu': 'Mus musculus',
            'rno': 'Rattus norvegicus',
            'dme': 'Drosophila melanogaster',
            'cel': 'Caenorhabditis elegans',
            'dre': 'Danio rerio',
        }

        # Full species names
        species_names = {
            'arabidopsis_thaliana': 'Arabidopsis thaliana',
            'medicago_truncatula': 'Medicago truncatula',
            'oryza_sativa': 'Oryza sativa',
            'zea_mays': 'Zea mays',
            'glycine_max': 'Glycine max',
            'homo_sapiens': 'Homo sapiens',
            'mus_musculus': 'Mus musculus',
        }

        filename_lower = filename.lower()

        # Check for species code (e.g., mirbase_ath_mature)
        for code, name in species_codes.items():
            if f'_{code}_' in filename_lower or filename_lower.startswith(f'{code}_') or filename_lower.startswith(f'mirbase_{code}'):
                species = name
                break

        # Check for full species name
        for name_key, name in species_names.items():
            if name_key in filename_lower:
                species = name
                break

        # Detect RNA types
        if 'mirna' in filename_lower or 'mature' in filename_lower or 'hairpin' in filename_lower:
            rna_types = ['miRNA']
        elif 'trna' in filename_lower:
            rna_types = ['tRNA']
        elif 'snorna' in filename_lower:
            rna_types = ['snoRNA']
        elif 'rrna' in filename_lower:
            rna_types = ['rRNA']
        elif 'pirna' in filename_lower:
            rna_types = ['piRNA']

        return species, rna_types

    def save_config(self):
        """Save reference configuration to JSON file"""
        data = {
            'references': list(self.references.values()),
            'metadata': self.metadata
        }

        self.reference_dir.mkdir(parents=True, exist_ok=True)

        with open(self.config_file, 'w') as f:
            json.dump(data, f, indent=4)

    def list_references(self) -> List[Dict]:
        """
        List all available references

        Returns:
            List of reference dictionaries
        """
        return list(self.references.values())

    def get_reference(self, ref_id: str) -> Optional[Dict]:
        """
        Get a specific reference by ID

        Args:
            ref_id: Reference identifier

        Returns:
            Reference dictionary or None if not found
        """
        return self.references.get(ref_id)

    def get_reference_path(self, ref_id: str) -> Optional[Path]:
        """
        Get the file path for a reference

        Args:
            ref_id: Reference identifier

        Returns:
            Path to reference file or None if not found
        """
        ref = self.get_reference(ref_id)
        if ref:
            # Check full_path first (for auto-detected refs in subdirs)
            if 'full_path' in ref:
                full_path = Path(ref['full_path'])
                if full_path.exists():
                    return full_path

            # Fall back to relative path
            path = self.reference_dir / ref['file']
            if path.exists():
                return path
        return None

    def get_index_prefix(self, ref_id: str) -> Optional[Path]:
        """
        Get the Bowtie index prefix for a reference

        Args:
            ref_id: Reference identifier

        Returns:
            Path to index prefix or None if not found
        """
        ref_path = self.get_reference_path(ref_id)
        if ref_path:
            index_prefix = ref_path.with_suffix('')
            # Check if index exists
            if (index_prefix.parent / f"{index_prefix.name}.1.ebwt").exists():
                return index_prefix
        return None

    def refresh(self):
        """Refresh the reference list by re-scanning directories"""
        self.references = {}
        self._load_config()
        self._auto_detect_references()

    def get_reference_names(self) -> Dict[str, str]:
        """
        Get a mapping of reference IDs to display names

        Returns:
            Dictionary of id -> name
        """
        return {ref_id: ref['name'] for ref_id, ref in self.references.items()}

    def get_species_list(self) -> List[str]:
        """
        Get list of available species

        Returns:
            List of species names
        """
        return list(set(ref.get('species', 'Unknown') for ref in self.references.values()))

    def add_reference(
        self,
        ref_id: str,
        name: str,
        species: str,
        file: str,
        description: str = "",
        source: str = "Custom",
        taxid: Optional[int] = None,
        rna_types: Optional[List[str]] = None,
        **kwargs
    ) -> bool:
        """
        Add a new reference to the configuration

        Args:
            ref_id: Unique identifier for the reference
            name: Display name
            species: Species name
            file: Filename (must exist in reference_dir)
            description: Optional description
            source: Data source (e.g., miRBase, RNAcentral, Custom)
            taxid: NCBI Taxonomy ID
            rna_types: List of RNA types included
            **kwargs: Additional metadata

        Returns:
            True if successful, False otherwise
        """
        # Check if file exists
        file_path = self.reference_dir / file
        if not file_path.exists():
            return False

        # Count sequences
        seq_count = self._count_sequences(file_path)

        self.references[ref_id] = {
            'id': ref_id,
            'name': name,
            'species': species,
            'file': file,
            'description': description,
            'source': source,
            'taxid': taxid,
            'rna_types': rna_types or ['unknown'],
            'sequences': seq_count,
            **kwargs
        }

        self.save_config()
        return True

    def remove_reference(self, ref_id: str) -> bool:
        """
        Remove a reference from the configuration (does not delete file)

        Args:
            ref_id: Reference identifier

        Returns:
            True if removed, False if not found
        """
        if ref_id in self.references:
            del self.references[ref_id]
            self.save_config()
            return True
        return False

    def _count_sequences(self, fasta_path: Path) -> int:
        """Count sequences in a FASTA file"""
        count = 0
        try:
            with open(fasta_path, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        count += 1
        except Exception:
            pass
        return count

    def get_reference_info(self, ref_id: str) -> Optional[pd.DataFrame]:
        """
        Get detailed information about a reference as a DataFrame

        Args:
            ref_id: Reference identifier

        Returns:
            DataFrame with reference info
        """
        ref = self.get_reference(ref_id)
        if ref:
            # Update sequence count if not set
            if ref.get('sequences', 0) == 0:
                path = self.get_reference_path(ref_id)
                if path:
                    ref['sequences'] = self._count_sequences(path)
                    self.save_config()

            return pd.DataFrame([ref]).T.rename(columns={0: 'Value'})
        return None

    def validate_references(self) -> Dict[str, bool]:
        """
        Validate that all reference files exist

        Returns:
            Dictionary of ref_id -> exists (bool)
        """
        results = {}
        for ref_id, ref in self.references.items():
            path = self.reference_dir / ref['file']
            results[ref_id] = path.exists()
        return results


# Convenience function for quick access
def get_reference_manager() -> ReferenceManager:
    """Get the default reference manager instance"""
    return ReferenceManager()
