"""
Novel miRNA Discovery Module for sRNAtlas
Identify unannotated small RNAs with miRNA-like characteristics
"""
import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import subprocess
import tempfile

try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False


def render_novel_mirna_page():
    """Render the novel miRNA discovery page"""
    st.header("🔍 Novel miRNA Discovery")
    st.markdown("Identify unannotated small RNAs with miRNA-like characteristics")
    
    tab1, tab2, tab3, tab4 = st.tabs([
        "📥 Input",
        "⚙️ Settings", 
        "🔬 Analysis",
        "📊 Results"
    ])
    
    with tab1:
        render_input_tab()
    
    with tab2:
        render_settings_tab()
    
    with tab3:
        render_analysis_tab()
    
    with tab4:
        render_results_tab()


def render_input_tab():
    """Input selection for novel miRNA discovery"""
    st.subheader("Input Data")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("#### Unaligned Reads")
        st.info("Use reads that did NOT align to known miRNAs")
        
        # Check for unaligned files from alignment step
        if 'alignment_results' in st.session_state and st.session_state.alignment_results:
            unaligned_files = []
            for result in st.session_state.alignment_results:
                if 'unaligned_file' in result and Path(result['unaligned_file']).exists():
                    unaligned_files.append(result['unaligned_file'])
            
            if unaligned_files:
                st.success(f"Found {len(unaligned_files)} unaligned read files")
                use_unaligned = st.checkbox("Use unaligned reads from alignment", value=True)
                if use_unaligned:
                    st.session_state.novel_input_files = unaligned_files
        
        # Manual upload option
        st.markdown("**Or upload FASTQ files:**")
        uploaded = st.file_uploader(
            "Upload unaligned reads",
            type=['fastq', 'fq', 'fastq.gz', 'fq.gz'],
            accept_multiple_files=True,
            key="novel_fastq_upload"
        )
    
    with col2:
        st.markdown("#### Reference Genome")
        st.info("Required for hairpin structure prediction")
        
        genome_source = st.radio(
            "Genome source",
            ["Use existing", "Upload FASTA"],
            key="novel_genome_source"
        )
        
        if genome_source == "Upload FASTA":
            genome_file = st.file_uploader(
                "Upload genome FASTA",
                type=['fasta', 'fa', 'fna'],
                key="novel_genome_upload"
            )


def render_settings_tab():
    """Settings for novel miRNA discovery"""
    st.subheader("Discovery Parameters")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("#### Size Selection")
        min_len = st.slider("Minimum length (nt)", 18, 25, 20)
        max_len = st.slider("Maximum length (nt)", 22, 30, 24)
        st.session_state.novel_min_len = min_len
        st.session_state.novel_max_len = max_len
        
        st.markdown("#### Read Filtering")
        min_reads = st.number_input("Minimum read count", 10, 10000, 100)
        st.session_state.novel_min_reads = min_reads
    
    with col2:
        st.markdown("#### miRNA Criteria")
        
        st.markdown("**5' Nucleotide Preference:**")
        prefer_u = st.checkbox("Prefer 5' Uracil (U/T)", value=True)
        st.session_state.novel_prefer_u = prefer_u
        
        st.markdown("**Hairpin Structure:**")
        check_hairpin = st.checkbox("Check for hairpin potential", value=True)
        st.session_state.novel_check_hairpin = check_hairpin
        
        if check_hairpin:
            min_mfe = st.slider(
                "Max MFE (kcal/mol)", 
                -50.0, -10.0, -20.0,
                help="More negative = more stable hairpin"
            )
            st.session_state.novel_min_mfe = min_mfe
    
    st.markdown("---")
    st.markdown("#### Filtering Criteria Summary")
    st.markdown(f"""
    - Length: {min_len}-{max_len} nt
    - Minimum reads: {min_reads}
    - 5' U preference: {'Yes' if prefer_u else 'No'}
    - Hairpin check: {'Yes' if check_hairpin else 'No'}
    """)


def render_analysis_tab():
    """Run novel miRNA discovery analysis"""
    st.subheader("Run Analysis")
    
    # Check prerequisites
    ready = True
    
    if not PYSAM_AVAILABLE:
        st.warning("⚠️ pysam not installed - some features limited")
    
    input_files = st.session_state.get('novel_input_files', [])
    if not input_files:
        st.error("❌ No input files selected")
        ready = False
    else:
        st.success(f"✅ {len(input_files)} input files ready")
    
    if ready:
        if st.button("🚀 Start Discovery", type="primary"):
            with st.spinner("Analyzing reads for novel miRNAs..."):
                results = run_novel_discovery(input_files)
                st.session_state.novel_results = results
                st.success("✅ Analysis complete!")
                st.rerun()


def run_novel_discovery(input_files: List[str]) -> Dict:
    """Run novel miRNA discovery pipeline"""
    results = {
        'candidates': [],
        'stats': {},
        'sequences': {}
    }
    
    min_len = st.session_state.get('novel_min_len', 20)
    max_len = st.session_state.get('novel_max_len', 24)
    min_reads = st.session_state.get('novel_min_reads', 100)
    prefer_u = st.session_state.get('novel_prefer_u', True)
    
    # Collect and count unique sequences
    seq_counts = {}
    
    for file_path in input_files:
        try:
            # Read FASTQ and count sequences
            import gzip
            opener = gzip.open if str(file_path).endswith('.gz') else open
            
            with opener(file_path, 'rt') as f:
                line_num = 0
                for line in f:
                    line_num += 1
                    if line_num % 4 == 2:  # Sequence line
                        seq = line.strip()
                        if min_len <= len(seq) <= max_len:
                            seq_counts[seq] = seq_counts.get(seq, 0) + 1
        except Exception as e:
            st.warning(f"Error reading {file_path}: {e}")
    
    # Filter by count and criteria
    candidates = []
    for seq, count in seq_counts.items():
        if count >= min_reads:
            candidate = {
                'sequence': seq,
                'length': len(seq),
                'count': count,
                'first_nt': seq[0],
                'u_bias': seq[0] in ['T', 'U'],
                'gc_content': (seq.count('G') + seq.count('C')) / len(seq) * 100
            }
            
            # Apply U preference filter
            if prefer_u and not candidate['u_bias']:
                candidate['score'] = 0.5
            else:
                candidate['score'] = 1.0
            
            # Score based on characteristics
            candidate['score'] *= min(count / 1000, 1.0)  # Read count contribution
            
            candidates.append(candidate)
    
    # Sort by score
    candidates.sort(key=lambda x: x['score'], reverse=True)
    
    results['candidates'] = candidates[:100]  # Top 100
    results['stats'] = {
        'total_unique': len(seq_counts),
        'passed_filters': len(candidates),
        'with_u_bias': sum(1 for c in candidates if c['u_bias'])
    }
    
    return results


def render_results_tab():
    """Display novel miRNA discovery results"""
    st.subheader("Discovery Results")
    
    if 'novel_results' not in st.session_state or not st.session_state.novel_results:
        st.info("Run analysis to see results")
        return
    
    results = st.session_state.novel_results
    
    # Summary stats
    col1, col2, col3 = st.columns(3)
    col1.metric("Total Unique Sequences", f"{results['stats']['total_unique']:,}")
    col2.metric("Passed Filters", results['stats']['passed_filters'])
    col3.metric("With 5' U Bias", results['stats']['with_u_bias'])
    
    st.markdown("---")
    
    # Candidates table
    if results['candidates']:
        st.markdown("#### Top Candidate Sequences")
        
        df = pd.DataFrame(results['candidates'])
        df = df[['sequence', 'length', 'count', 'first_nt', 'u_bias', 'gc_content', 'score']]
        df.columns = ['Sequence', 'Length', 'Reads', '5\' nt', 'U Bias', 'GC %', 'Score']
        df['GC %'] = df['GC %'].round(1)
        df['Score'] = df['Score'].round(3)
        df['U Bias'] = df['U Bias'].map({True: '✅', False: '❌'})
        
        st.dataframe(df, width="stretch", hide_index=True)
        
        # Download button
        csv = df.to_csv(index=False)
        st.download_button(
            "📥 Download Candidates (CSV)",
            csv,
            file_name="novel_mirna_candidates.csv",
            mime="text/csv"
        )
    else:
        st.warning("No candidates found matching criteria")
