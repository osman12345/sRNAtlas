#!/usr/bin/env python3
"""
Reference Data Download Script for sRNAtlas

This script downloads reference sequences from miRBase and RNAcentral
for bundling with the application. Run this script to update reference data.

Usage:
    python download_references.py --species ath mtr hsa --output ../data/references/

Author: sRNA WebTool Team
"""

import argparse
import requests
import time
import sys
from pathlib import Path
from typing import List, Optional
import json

# Species configurations
SPECIES_CONFIG = {
    # Plants
    'ath': {
        'name': 'Arabidopsis thaliana',
        'common': 'Thale Cress',
        'taxid': 3702,
        'mirbase_code': 'ath',
        'kingdom': 'plant'
    },
    'mtr': {
        'name': 'Medicago truncatula',
        'common': 'Barrel Medic',
        'taxid': 3880,
        'mirbase_code': 'mtr',
        'kingdom': 'plant'
    },
    'osa': {
        'name': 'Oryza sativa',
        'common': 'Rice',
        'taxid': 39947,
        'mirbase_code': 'osa',
        'kingdom': 'plant'
    },
    'zma': {
        'name': 'Zea mays',
        'common': 'Maize',
        'taxid': 4577,
        'mirbase_code': 'zma',
        'kingdom': 'plant'
    },
    'sly': {
        'name': 'Solanum lycopersicum',
        'common': 'Tomato',
        'taxid': 4081,
        'mirbase_code': 'sly',
        'kingdom': 'plant'
    },
    'gma': {
        'name': 'Glycine max',
        'common': 'Soybean',
        'taxid': 3847,
        'mirbase_code': 'gma',
        'kingdom': 'plant'
    },
    'vvi': {
        'name': 'Vitis vinifera',
        'common': 'Grape',
        'taxid': 29760,
        'mirbase_code': 'vvi',
        'kingdom': 'plant'
    },
    # Animals
    'hsa': {
        'name': 'Homo sapiens',
        'common': 'Human',
        'taxid': 9606,
        'mirbase_code': 'hsa',
        'kingdom': 'animal'
    },
    'mmu': {
        'name': 'Mus musculus',
        'common': 'Mouse',
        'taxid': 10090,
        'mirbase_code': 'mmu',
        'kingdom': 'animal'
    },
    'dre': {
        'name': 'Danio rerio',
        'common': 'Zebrafish',
        'taxid': 7955,
        'mirbase_code': 'dre',
        'kingdom': 'animal'
    },
    'dme': {
        'name': 'Drosophila melanogaster',
        'common': 'Fruit Fly',
        'taxid': 7227,
        'mirbase_code': 'dme',
        'kingdom': 'animal'
    },
    'cel': {
        'name': 'Caenorhabditis elegans',
        'common': 'Roundworm',
        'taxid': 6239,
        'mirbase_code': 'cel',
        'kingdom': 'animal'
    },
    'rno': {
        'name': 'Rattus norvegicus',
        'common': 'Rat',
        'taxid': 10116,
        'mirbase_code': 'rno',
        'kingdom': 'animal'
    },
}

# RNA types available from RNAcentral
RNA_TYPES = ['miRNA', 'tRNA', 'rRNA', 'snRNA', 'snoRNA', 'piRNA', 'lncRNA']


def download_mirbase(species_code: str, output_dir: Path, seq_type: str = 'mature') -> Optional[Path]:
    """
    Download miRNA sequences from miRBase.

    Args:
        species_code: 3-letter species code (e.g., 'ath', 'hsa')
        output_dir: Directory to save files
        seq_type: 'mature' or 'hairpin'

    Returns:
        Path to downloaded file or None if failed
    """
    if species_code not in SPECIES_CONFIG:
        print(f"  ❌ Unknown species code: {species_code}")
        return None

    mirbase_code = SPECIES_CONFIG[species_code]['mirbase_code']
    url = f"https://www.mirbase.org/download/{seq_type}.fa"

    print(f"  📥 Downloading miRBase {seq_type} for {species_code}...")

    try:
        response = requests.get(url, timeout=60)
        response.raise_for_status()

        # Filter for species
        lines = response.text.split('\n')
        filtered_lines = []
        include_seq = False

        for line in lines:
            if line.startswith('>'):
                # Check if this sequence belongs to our species
                include_seq = line.lower().startswith(f'>{mirbase_code}-')
                if include_seq:
                    filtered_lines.append(line)
            elif include_seq and line.strip():
                filtered_lines.append(line)

        if not filtered_lines:
            print(f"  ⚠️ No sequences found for {species_code}")
            return None

        # Save to file
        output_file = output_dir / f"mirbase_{species_code}_{seq_type}.fa"
        with open(output_file, 'w') as f:
            f.write('\n'.join(filtered_lines) + '\n')

        seq_count = len([l for l in filtered_lines if l.startswith('>')])
        print(f"  ✅ Saved {seq_count} {seq_type} miRNAs to {output_file.name}")
        return output_file

    except Exception as e:
        print(f"  ❌ Failed to download miRBase {seq_type}: {e}")
        return None


def download_rnacentral(species_code: str, output_dir: Path, rna_type: str, max_seqs: int = 5000) -> Optional[Path]:
    """
    Download RNA sequences from RNAcentral.

    Args:
        species_code: 3-letter species code
        output_dir: Directory to save files
        rna_type: RNA type (e.g., 'tRNA', 'rRNA')
        max_seqs: Maximum number of sequences to download

    Returns:
        Path to downloaded file or None if failed
    """
    if species_code not in SPECIES_CONFIG:
        print(f"  ❌ Unknown species code: {species_code}")
        return None

    config = SPECIES_CONFIG[species_code]
    taxid = config['taxid']
    species_name = config['name'].replace(' ', '%20')

    print(f"  📥 Downloading RNAcentral {rna_type} for {species_code}...")

    # RNAcentral text search API
    base_url = "https://rnacentral.org/api/v1/rna"

    sequences = []
    page = 1
    page_size = 100

    try:
        while len(sequences) < max_seqs:
            params = {
                'taxid': taxid,
                'rna_type': rna_type,
                'page': page,
                'page_size': page_size,
                'format': 'json'
            }

            response = requests.get(base_url, params=params, timeout=60)

            if response.status_code != 200:
                break

            data = response.json()
            results = data.get('results', [])

            if not results:
                break

            for item in results:
                rna_id = item.get('rnacentral_id', '')
                seq = item.get('sequence', '')
                desc = item.get('description', '')

                if rna_id and seq and len(seq) <= 200:  # Small RNA length filter
                    header = f">{rna_id}|{rna_type}|{config['name']}|{desc[:50]}"
                    sequences.append((header, seq))

            if not data.get('next'):
                break

            page += 1
            time.sleep(0.5)  # Rate limiting

            if page > 50:  # Safety limit
                break

        if not sequences:
            print(f"  ⚠️ No {rna_type} sequences found for {species_code}")
            return None

        # Save to file
        output_file = output_dir / f"rnacentral_{species_code}_{rna_type.lower()}.fa"
        with open(output_file, 'w') as f:
            for header, seq in sequences:
                f.write(f"{header}\n{seq}\n")

        print(f"  ✅ Saved {len(sequences)} {rna_type} sequences to {output_file.name}")
        return output_file

    except Exception as e:
        print(f"  ❌ Failed to download RNAcentral {rna_type}: {e}")
        return None


def create_combined_reference(species_code: str, output_dir: Path, files: List[Path]) -> Optional[Path]:
    """
    Combine multiple FASTA files into a single reference file.

    Args:
        species_code: Species code
        output_dir: Output directory
        files: List of FASTA files to combine

    Returns:
        Path to combined file or None
    """
    if not files:
        return None

    output_file = output_dir / f"{species_code}_smallRNA_combined.fa"

    print(f"  📦 Creating combined reference for {species_code}...")

    total_seqs = 0
    with open(output_file, 'w') as out:
        for fasta_file in files:
            if fasta_file and fasta_file.exists():
                with open(fasta_file, 'r') as f:
                    for line in f:
                        out.write(line)
                        if line.startswith('>'):
                            total_seqs += 1

    print(f"  ✅ Combined reference: {total_seqs} total sequences")
    return output_file


def update_references_json(output_dir: Path, downloaded_refs: dict):
    """Update references.json with downloaded references."""
    json_file = output_dir / 'references.json'

    # Load existing
    if json_file.exists():
        with open(json_file, 'r') as f:
            data = json.load(f)
    else:
        data = {'references': [], 'metadata': {}}

    # Create reference entries
    existing_ids = {ref['id'] for ref in data['references']}

    for species_code, files_info in downloaded_refs.items():
        config = SPECIES_CONFIG[species_code]

        for file_info in files_info:
            ref_id = file_info['id']
            if ref_id in existing_ids:
                # Update existing
                for ref in data['references']:
                    if ref['id'] == ref_id:
                        ref.update(file_info)
                        break
            else:
                # Add new
                data['references'].append(file_info)
                existing_ids.add(ref_id)

    # Update metadata
    data['metadata'] = {
        'last_updated': time.strftime('%Y-%m-%d'),
        'total_references': len(data['references']),
        'instructions': 'Reference database for sRNAtlas'
    }

    with open(json_file, 'w') as f:
        json.dump(data, f, indent=4)

    print(f"\n✅ Updated {json_file}")


def main():
    parser = argparse.ArgumentParser(description='Download sRNA reference sequences')
    parser.add_argument('--species', nargs='+', default=['ath', 'mtr', 'hsa'],
                       help='Species codes to download (default: ath mtr hsa)')
    parser.add_argument('--output', type=str, default='../data/references',
                       help='Output directory')
    parser.add_argument('--rna-types', nargs='+', default=['miRNA'],
                       help='RNA types to download (default: miRNA)')
    parser.add_argument('--include-all-rna', action='store_true',
                       help='Download all RNA types (miRNA, tRNA, rRNA, snRNA, snoRNA)')
    parser.add_argument('--mirbase-only', action='store_true',
                       help='Only download from miRBase (faster, miRNA only)')

    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create subdirectories
    mirbase_dir = output_dir / 'mirbase'
    rnacentral_dir = output_dir / 'rnacentral'
    mirbase_dir.mkdir(exist_ok=True)
    rnacentral_dir.mkdir(exist_ok=True)

    rna_types = RNA_TYPES if args.include_all_rna else args.rna_types

    print("=" * 60)
    print("sRNA Reference Database Downloader")
    print("=" * 60)
    print(f"Species: {', '.join(args.species)}")
    print(f"RNA Types: {', '.join(rna_types)}")
    print(f"Output: {output_dir.absolute()}")
    print("=" * 60)

    downloaded_refs = {}

    for species_code in args.species:
        if species_code not in SPECIES_CONFIG:
            print(f"\n⚠️ Unknown species: {species_code}, skipping...")
            continue

        config = SPECIES_CONFIG[species_code]
        print(f"\n📌 Processing {config['name']} ({species_code})...")

        downloaded_refs[species_code] = []
        downloaded_files = []

        # Download miRBase (always for miRNA)
        if 'miRNA' in rna_types:
            mature_file = download_mirbase(species_code, mirbase_dir, 'mature')
            hairpin_file = download_mirbase(species_code, mirbase_dir, 'hairpin')

            if mature_file:
                downloaded_files.append(mature_file)
                seq_count = sum(1 for line in open(mature_file) if line.startswith('>'))
                downloaded_refs[species_code].append({
                    'id': f'mirbase_{species_code}_mature',
                    'name': f"{config['name']} - miRNA (mature)",
                    'species': config['name'],
                    'common_name': config['common'],
                    'taxid': config['taxid'],
                    'file': f"mirbase/{mature_file.name}",
                    'description': f"Mature miRNA sequences from miRBase",
                    'source': 'miRBase',
                    'rna_types': ['miRNA'],
                    'sequences': seq_count
                })

            if hairpin_file:
                downloaded_files.append(hairpin_file)
                seq_count = sum(1 for line in open(hairpin_file) if line.startswith('>'))
                downloaded_refs[species_code].append({
                    'id': f'mirbase_{species_code}_hairpin',
                    'name': f"{config['name']} - miRNA (hairpin)",
                    'species': config['name'],
                    'common_name': config['common'],
                    'taxid': config['taxid'],
                    'file': f"mirbase/{hairpin_file.name}",
                    'description': f"Hairpin/precursor miRNA sequences from miRBase",
                    'source': 'miRBase',
                    'rna_types': ['miRNA'],
                    'sequences': seq_count
                })

        # Download other RNA types from RNAcentral
        if not args.mirbase_only:
            for rna_type in rna_types:
                if rna_type == 'miRNA':
                    continue  # Already downloaded from miRBase

                rna_file = download_rnacentral(species_code, rnacentral_dir, rna_type)
                if rna_file:
                    downloaded_files.append(rna_file)
                    seq_count = sum(1 for line in open(rna_file) if line.startswith('>'))
                    downloaded_refs[species_code].append({
                        'id': f'rnacentral_{species_code}_{rna_type.lower()}',
                        'name': f"{config['name']} - {rna_type}",
                        'species': config['name'],
                        'common_name': config['common'],
                        'taxid': config['taxid'],
                        'file': f"rnacentral/{rna_file.name}",
                        'description': f"{rna_type} sequences from RNAcentral",
                        'source': 'RNAcentral',
                        'rna_types': [rna_type],
                        'sequences': seq_count
                    })

                time.sleep(1)  # Rate limiting between RNA types

    # Update references.json
    update_references_json(output_dir, downloaded_refs)

    print("\n" + "=" * 60)
    print("Download complete!")
    print("=" * 60)


if __name__ == '__main__':
    main()
