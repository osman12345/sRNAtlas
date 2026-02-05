"""
Batch Processing Module for sRNAtlas
Run multiple analyses and manage job queues
"""
import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Callable
from datetime import datetime
import json
import time
import threading
import queue
from dataclasses import dataclass, asdict
from enum import Enum

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from config import config


class JobStatus(Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


@dataclass
class BatchJob:
    """Represents a batch processing job"""
    job_id: str
    job_type: str  # 'trimming', 'alignment', 'counting', 'de_analysis', 'full_pipeline'
    status: JobStatus
    created_at: str
    started_at: Optional[str] = None
    completed_at: Optional[str] = None
    progress: float = 0.0
    message: str = ""
    parameters: Dict = None
    input_files: List[str] = None
    output_files: List[str] = None
    error: Optional[str] = None
    results: Optional[Dict] = None


class BatchProcessor:
    """Manages batch job processing"""

    def __init__(self):
        self.jobs: Dict[str, BatchJob] = {}
        self.job_queue = queue.Queue()
        self.is_processing = False
        self._load_jobs()

    def _load_jobs(self):
        """Load jobs from storage"""
        jobs_file = Path(config.paths.output_dir) / "batch_jobs.json"
        if jobs_file.exists():
            try:
                with open(jobs_file, 'r') as f:
                    data = json.load(f)
                    for job_id, job_data in data.items():
                        job_data['status'] = JobStatus(job_data['status'])
                        self.jobs[job_id] = BatchJob(**job_data)
            except Exception:
                pass

    def _save_jobs(self):
        """Save jobs to storage"""
        jobs_file = Path(config.paths.output_dir) / "batch_jobs.json"
        jobs_file.parent.mkdir(parents=True, exist_ok=True)

        data = {}
        for job_id, job in self.jobs.items():
            job_dict = asdict(job)
            job_dict['status'] = job.status.value
            data[job_id] = job_dict

        with open(jobs_file, 'w') as f:
            json.dump(data, f, indent=2)

    def create_job(
        self,
        job_type: str,
        parameters: Dict,
        input_files: List[str]
    ) -> BatchJob:
        """Create a new batch job"""
        job_id = f"{job_type}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

        job = BatchJob(
            job_id=job_id,
            job_type=job_type,
            status=JobStatus.PENDING,
            created_at=datetime.now().isoformat(),
            parameters=parameters,
            input_files=input_files,
            output_files=[]
        )

        self.jobs[job_id] = job
        self._save_jobs()

        return job

    def update_job(self, job_id: str, **kwargs):
        """Update job attributes"""
        if job_id in self.jobs:
            job = self.jobs[job_id]
            for key, value in kwargs.items():
                if hasattr(job, key):
                    setattr(job, key, value)
            self._save_jobs()

    def get_job(self, job_id: str) -> Optional[BatchJob]:
        """Get job by ID"""
        return self.jobs.get(job_id)

    def list_jobs(self, status: Optional[JobStatus] = None) -> List[BatchJob]:
        """List all jobs, optionally filtered by status"""
        jobs = list(self.jobs.values())
        if status:
            jobs = [j for j in jobs if j.status == status]
        return sorted(jobs, key=lambda j: j.created_at, reverse=True)

    def cancel_job(self, job_id: str):
        """Cancel a pending job"""
        if job_id in self.jobs:
            job = self.jobs[job_id]
            if job.status == JobStatus.PENDING:
                job.status = JobStatus.CANCELLED
                self._save_jobs()

    def delete_job(self, job_id: str):
        """Delete a job"""
        if job_id in self.jobs:
            del self.jobs[job_id]
            self._save_jobs()


def create_pipeline_config(
    samples: List[Dict],
    reference_db: str,
    output_dir: str,
    run_trimming: bool = True,
    run_qc: bool = True,
    run_alignment: bool = True,
    run_counting: bool = True,
    run_de: bool = True,
    trimming_params: Dict = None,
    alignment_params: Dict = None,
    de_params: Dict = None
) -> Dict:
    """
    Create a pipeline configuration for batch processing

    Args:
        samples: List of sample dictionaries with 'name', 'fastq', 'group'
        reference_db: Path to reference database
        output_dir: Output directory
        run_*: Flags for each pipeline step
        *_params: Parameters for each step

    Returns:
        Pipeline configuration dictionary
    """
    return {
        'samples': samples,
        'reference_db': reference_db,
        'output_dir': output_dir,
        'steps': {
            'trimming': run_trimming,
            'qc': run_qc,
            'alignment': run_alignment,
            'counting': run_counting,
            'de_analysis': run_de
        },
        'parameters': {
            'trimming': trimming_params or {},
            'alignment': alignment_params or {},
            'de': de_params or {}
        }
    }


def render_batch_page():
    """Render the batch processing page"""
    st.header("⚡ Batch Processing")

    st.markdown("""
    Run multiple analyses in batch mode or process entire pipelines automatically.
    """)

    # Initialize batch processor
    if 'batch_processor' not in st.session_state:
        st.session_state.batch_processor = BatchProcessor()

    processor = st.session_state.batch_processor

    tab1, tab2, tab3 = st.tabs([
        "🚀 New Batch",
        "📋 Job Queue",
        "📊 Results"
    ])

    with tab1:
        render_new_batch(processor)

    with tab2:
        render_job_queue(processor)

    with tab3:
        render_batch_results(processor)


def render_new_batch(processor: BatchProcessor):
    """Render new batch job creation"""
    st.subheader("Create New Batch Job")

    job_type = st.selectbox(
        "Job Type",
        options=[
            "Full Pipeline",
            "Trimming Only",
            "Alignment Only",
            "DE Analysis Only",
            "Multi-comparison DE"
        ]
    )

    if job_type == "Full Pipeline":
        render_full_pipeline_config(processor)
    elif job_type == "Trimming Only":
        render_trimming_batch_config(processor)
    elif job_type == "Alignment Only":
        render_alignment_batch_config(processor)
    elif job_type == "DE Analysis Only":
        render_de_batch_config(processor)
    else:
        render_multi_comparison_config(processor)


def render_full_pipeline_config(processor: BatchProcessor):
    """Render full pipeline configuration"""
    st.markdown("### Full Pipeline Configuration")

    st.markdown("""
    Run the complete sRNA-seq pipeline:
    **Trimming → QC → Alignment → Counting → DE Analysis**
    """)

    # Sample sheet upload
    st.markdown("#### 1. Sample Sheet")

    sample_sheet_example = """SampleID,FASTQ,Group,Replicate
Control_1,/path/to/control_1.fastq.gz,Control,1
Control_2,/path/to/control_2.fastq.gz,Control,2
Control_3,/path/to/control_3.fastq.gz,Control,3
Treatment_1,/path/to/treatment_1.fastq.gz,Treatment,1
Treatment_2,/path/to/treatment_2.fastq.gz,Treatment,2
Treatment_3,/path/to/treatment_3.fastq.gz,Treatment,3"""

    with st.expander("Sample Sheet Format"):
        st.code(sample_sheet_example, language="csv")

    sample_file = st.file_uploader(
        "Upload Sample Sheet (CSV)",
        type=['csv', 'tsv', 'txt']
    )

    samples = []
    if sample_file:
        content = sample_file.getvalue().decode('utf-8')
        sep = '\t' if '\t' in content.split('\n')[0] else ','
        sample_file.seek(0)
        sample_df = pd.read_csv(sample_file, sep=sep)

        st.dataframe(sample_df, width="stretch")

        # Validate required columns
        required = ['SampleID', 'Group']
        if all(col in sample_df.columns for col in required):
            samples = sample_df.to_dict('records')
            st.success(f"✅ Loaded {len(samples)} samples")
        else:
            st.error(f"Missing required columns: {required}")

    # Reference database
    st.markdown("#### 2. Reference Database")

    available_refs = st.session_state.get('available_indices', [])

    if available_refs:
        reference = st.selectbox(
            "Select Reference Database",
            options=available_refs
        )
    else:
        reference = st.text_input(
            "Reference Database Path",
            help="Path to Bowtie index prefix"
        )

    # Pipeline steps
    st.markdown("#### 3. Pipeline Steps")

    col1, col2 = st.columns(2)

    with col1:
        run_trimming = st.checkbox("Adapter Trimming", value=True)
        run_qc = st.checkbox("Quality Control", value=True)
        run_alignment = st.checkbox("Alignment", value=True)

    with col2:
        run_counting = st.checkbox("Read Counting", value=True)
        run_de = st.checkbox("DE Analysis", value=True)

    # Parameters
    st.markdown("#### 4. Parameters")

    with st.expander("Trimming Parameters"):
        adapter_3 = st.text_input(
            "3' Adapter",
            value="TGGAATTCTCGGGTGCCAAGG"
        )
        min_length = st.slider("Min Length", 10, 30, 18)
        max_length = st.slider("Max Length", 25, 50, 35)

    with st.expander("Alignment Parameters"):
        mismatches = st.slider("Mismatches", 0, 3, 1)
        multi_alignments = st.slider("Max Multi-alignments", 1, 50, 10)

    with st.expander("DE Analysis Parameters"):
        fdr_threshold = st.slider("FDR Threshold", 0.01, 0.20, 0.05)
        lfc_threshold = st.slider("LFC Threshold", 0.0, 2.0, 0.585)

    # Output directory
    st.markdown("#### 5. Output")

    output_dir = st.text_input(
        "Output Directory",
        value=str(Path(config.paths.output_dir) / "batch_results")
    )

    # Create job
    st.divider()

    if st.button("🚀 Submit Batch Job", type="primary", width="stretch"):
        if not samples:
            st.error("Please upload a sample sheet")
            return

        if not reference:
            st.error("Please specify a reference database")
            return

        # Create pipeline config
        pipeline_config = create_pipeline_config(
            samples=samples,
            reference_db=reference,
            output_dir=output_dir,
            run_trimming=run_trimming,
            run_qc=run_qc,
            run_alignment=run_alignment,
            run_counting=run_counting,
            run_de=run_de,
            trimming_params={
                'adapter_3': adapter_3,
                'min_length': min_length,
                'max_length': max_length
            },
            alignment_params={
                'mismatches': mismatches,
                'max_alignments': multi_alignments
            },
            de_params={
                'fdr_threshold': fdr_threshold,
                'lfc_threshold': lfc_threshold
            }
        )

        # Create job
        job = processor.create_job(
            job_type='full_pipeline',
            parameters=pipeline_config,
            input_files=[s.get('FASTQ', s.get('SampleID')) for s in samples]
        )

        st.success(f"✅ Created batch job: {job.job_id}")
        st.info("Go to Job Queue tab to monitor progress")


def render_trimming_batch_config(processor: BatchProcessor):
    """Render trimming batch configuration"""
    st.markdown("### Batch Trimming Configuration")

    # Upload multiple FASTQ files
    fastq_files = st.file_uploader(
        "Upload FASTQ Files",
        type=['fastq', 'fq', 'gz'],
        accept_multiple_files=True
    )

    if fastq_files:
        st.info(f"Selected {len(fastq_files)} files")

        # Adapter settings
        adapter_3 = st.text_input("3' Adapter", value="TGGAATTCTCGGGTGCCAAGG")

        col1, col2 = st.columns(2)
        with col1:
            min_length = st.slider("Min Length", 10, 30, 18)
        with col2:
            max_length = st.slider("Max Length", 25, 50, 35)

        if st.button("🚀 Submit Trimming Job", type="primary"):
            job = processor.create_job(
                job_type='trimming',
                parameters={
                    'adapter_3': adapter_3,
                    'min_length': min_length,
                    'max_length': max_length
                },
                input_files=[f.name for f in fastq_files]
            )
            st.success(f"✅ Created job: {job.job_id}")


def render_alignment_batch_config(processor: BatchProcessor):
    """Render alignment batch configuration"""
    st.markdown("### Batch Alignment Configuration")

    # Upload FASTQ files
    fastq_files = st.file_uploader(
        "Upload FASTQ Files",
        type=['fastq', 'fq', 'gz'],
        accept_multiple_files=True,
        key="align_fastq"
    )

    # Reference selection
    reference = st.text_input("Reference Database Path")

    # Alignment parameters
    col1, col2 = st.columns(2)
    with col1:
        mismatches = st.slider("Mismatches", 0, 3, 1)
    with col2:
        multi_alignments = st.slider("Max Multi-alignments", 1, 50, 10)

    if fastq_files and reference:
        if st.button("🚀 Submit Alignment Job", type="primary"):
            job = processor.create_job(
                job_type='alignment',
                parameters={
                    'reference': reference,
                    'mismatches': mismatches,
                    'max_alignments': multi_alignments
                },
                input_files=[f.name for f in fastq_files]
            )
            st.success(f"✅ Created job: {job.job_id}")


def render_de_batch_config(processor: BatchProcessor):
    """Render DE analysis batch configuration"""
    st.markdown("### Batch DE Analysis Configuration")

    # Upload count matrix
    count_file = st.file_uploader(
        "Upload Count Matrix (CSV)",
        type=['csv', 'tsv', 'txt']
    )

    # Upload metadata
    meta_file = st.file_uploader(
        "Upload Sample Metadata (CSV)",
        type=['csv', 'tsv', 'txt']
    )

    if count_file and meta_file:
        # Load and preview
        counts = pd.read_csv(count_file, index_col=0)
        metadata = pd.read_csv(meta_file)

        st.write("Count Matrix:", counts.shape)
        st.write("Metadata:", metadata.shape)

        # DE parameters
        design_factor = st.selectbox(
            "Design Factor",
            options=metadata.columns.tolist()
        )

        fdr_threshold = st.slider("FDR Threshold", 0.01, 0.20, 0.05)

        if st.button("🚀 Submit DE Job", type="primary"):
            job = processor.create_job(
                job_type='de_analysis',
                parameters={
                    'design_factor': design_factor,
                    'fdr_threshold': fdr_threshold
                },
                input_files=[count_file.name, meta_file.name]
            )
            st.success(f"✅ Created job: {job.job_id}")


def render_multi_comparison_config(processor: BatchProcessor):
    """Render multi-comparison DE configuration"""
    st.markdown("### Multi-Comparison DE Analysis")

    st.markdown("""
    Run multiple DE comparisons in one batch job.
    Useful for time-series or multi-condition experiments.
    """)

    # Upload files
    count_file = st.file_uploader(
        "Count Matrix",
        type=['csv', 'tsv'],
        key="multi_counts"
    )

    meta_file = st.file_uploader(
        "Metadata",
        type=['csv', 'tsv'],
        key="multi_meta"
    )

    if count_file and meta_file:
        metadata = pd.read_csv(meta_file)

        # Select factor
        factor = st.selectbox(
            "Comparison Factor",
            options=metadata.columns.tolist()
        )

        groups = metadata[factor].unique().tolist()

        # Define comparisons
        st.markdown("#### Define Comparisons")

        comparisons = []
        n_comparisons = st.number_input("Number of Comparisons", 1, 10, 1)

        for i in range(int(n_comparisons)):
            col1, col2 = st.columns(2)
            with col1:
                treatment = st.selectbox(
                    f"Treatment {i+1}",
                    options=groups,
                    key=f"treat_{i}"
                )
            with col2:
                control = st.selectbox(
                    f"Control {i+1}",
                    options=groups,
                    key=f"ctrl_{i}"
                )

            comparisons.append({
                'name': f"{treatment}_vs_{control}",
                'treatment': treatment,
                'control': control
            })

        if st.button("🚀 Submit Multi-Comparison Job", type="primary"):
            job = processor.create_job(
                job_type='multi_comparison',
                parameters={
                    'factor': factor,
                    'comparisons': comparisons,
                    'fdr_threshold': 0.05
                },
                input_files=[count_file.name, meta_file.name]
            )
            st.success(f"✅ Created job: {job.job_id}")


def render_job_queue(processor: BatchProcessor):
    """Render job queue view"""
    st.subheader("Job Queue")

    # Refresh button
    if st.button("🔄 Refresh"):
        st.rerun()

    jobs = processor.list_jobs()

    if not jobs:
        st.info("No jobs in queue")
        return

    # Summary
    col1, col2, col3, col4 = st.columns(4)

    pending = sum(1 for j in jobs if j.status == JobStatus.PENDING)
    running = sum(1 for j in jobs if j.status == JobStatus.RUNNING)
    completed = sum(1 for j in jobs if j.status == JobStatus.COMPLETED)
    failed = sum(1 for j in jobs if j.status == JobStatus.FAILED)

    with col1:
        st.metric("Pending", pending)
    with col2:
        st.metric("Running", running)
    with col3:
        st.metric("Completed", completed)
    with col4:
        st.metric("Failed", failed)

    # Job list
    st.divider()

    for job in jobs:
        status_emoji = {
            JobStatus.PENDING: "⏳",
            JobStatus.RUNNING: "🔄",
            JobStatus.COMPLETED: "✅",
            JobStatus.FAILED: "❌",
            JobStatus.CANCELLED: "🚫"
        }

        with st.expander(f"{status_emoji[job.status]} {job.job_id} - {job.job_type}"):
            col1, col2 = st.columns(2)

            with col1:
                st.write(f"**Status:** {job.status.value}")
                st.write(f"**Created:** {job.created_at}")
                if job.started_at:
                    st.write(f"**Started:** {job.started_at}")
                if job.completed_at:
                    st.write(f"**Completed:** {job.completed_at}")

            with col2:
                st.write(f"**Progress:** {job.progress * 100:.1f}%")
                st.progress(job.progress)
                if job.message:
                    st.write(f"**Message:** {job.message}")

            if job.error:
                st.error(f"Error: {job.error}")

            # Actions
            col1, col2, col3 = st.columns(3)

            with col1:
                if job.status == JobStatus.PENDING:
                    if st.button("Cancel", key=f"cancel_{job.job_id}"):
                        processor.cancel_job(job.job_id)
                        st.rerun()

            with col2:
                if job.status in [JobStatus.COMPLETED, JobStatus.FAILED, JobStatus.CANCELLED]:
                    if st.button("Delete", key=f"delete_{job.job_id}"):
                        processor.delete_job(job.job_id)
                        st.rerun()

            with col3:
                if job.status == JobStatus.COMPLETED and job.output_files:
                    st.write(f"Output: {len(job.output_files)} files")


def render_batch_results(processor: BatchProcessor):
    """Render batch results view"""
    st.subheader("Batch Results")

    completed_jobs = processor.list_jobs(JobStatus.COMPLETED)

    if not completed_jobs:
        st.info("No completed jobs yet")
        return

    # Select job
    job_options = {j.job_id: j for j in completed_jobs}
    selected_id = st.selectbox(
        "Select Completed Job",
        options=list(job_options.keys())
    )

    job = job_options[selected_id]

    # Show results
    st.markdown(f"### Results for {job.job_id}")

    col1, col2 = st.columns(2)

    with col1:
        st.write(f"**Type:** {job.job_type}")
        st.write(f"**Completed:** {job.completed_at}")

    with col2:
        st.write(f"**Input Files:** {len(job.input_files or [])}")
        st.write(f"**Output Files:** {len(job.output_files or [])}")

    # Output files
    if job.output_files:
        st.markdown("#### Output Files")

        for of in job.output_files:
            of_path = Path(of)
            if of_path.exists():
                col1, col2 = st.columns([3, 1])
                with col1:
                    st.text(of_path.name)
                with col2:
                    with open(of_path, 'rb') as f:
                        st.download_button(
                            "Download",
                            f.read(),
                            of_path.name,
                            key=f"dl_{of}"
                        )

    # Results summary
    if job.results:
        st.markdown("#### Results Summary")
        st.json(job.results)
