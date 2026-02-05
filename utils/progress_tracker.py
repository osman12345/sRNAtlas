"""
Progress Tracker for sRNAtlas Pipeline
Tracks analysis progress across modules
"""
import streamlit as st
from dataclasses import dataclass
from typing import Dict, List, Optional
from datetime import datetime


@dataclass
class PipelineStep:
    """Represents a single pipeline step"""
    name: str
    status: str  # 'pending', 'running', 'completed', 'failed', 'skipped'
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    message: str = ""


PIPELINE_STEPS = [
    ("project", "Project Setup"),
    ("qc", "Quality Control"),
    ("trimming", "Adapter Trimming"),
    ("database", "Reference Database"),
    ("alignment", "Alignment"),
    ("post_qc", "Post-Alignment QC"),
    ("counting", "Read Counting"),
    ("de_analysis", "DE Analysis"),
    ("targets", "Target Prediction"),
    ("enrichment", "GO/Pathway"),
]


def init_progress_tracker():
    """Initialize progress tracker in session state"""
    if 'pipeline_progress' not in st.session_state:
        st.session_state.pipeline_progress = {
            step_id: PipelineStep(name=step_name, status='pending')
            for step_id, step_name in PIPELINE_STEPS
        }
    if 'pipeline_start_time' not in st.session_state:
        st.session_state.pipeline_start_time = None


def update_step_status(step_id: str, status: str, message: str = ""):
    """Update the status of a pipeline step"""
    init_progress_tracker()
    
    if step_id in st.session_state.pipeline_progress:
        step = st.session_state.pipeline_progress[step_id]
        step.status = status
        step.message = message
        
        if status == 'running':
            step.start_time = datetime.now()
        elif status in ['completed', 'failed']:
            step.end_time = datetime.now()


def get_progress_percentage() -> float:
    """Calculate overall pipeline progress percentage"""
    init_progress_tracker()
    
    completed = sum(
        1 for step in st.session_state.pipeline_progress.values()
        if step.status in ['completed', 'skipped']
    )
    total = len(st.session_state.pipeline_progress)
    
    return (completed / total) * 100 if total > 0 else 0


def render_progress_bar():
    """Render the pipeline progress bar"""
    init_progress_tracker()
    
    progress = get_progress_percentage()
    
    # Progress bar
    st.progress(progress / 100)
    st.caption(f"Pipeline Progress: {progress:.0f}%")


def render_progress_overview():
    """Render detailed progress overview"""
    init_progress_tracker()
    
    st.subheader("📊 Pipeline Progress")
    
    # Overall progress
    progress = get_progress_percentage()
    col1, col2, col3 = st.columns(3)
    
    with col1:
        completed = sum(1 for s in st.session_state.pipeline_progress.values() if s.status == 'completed')
        st.metric("Completed", f"{completed}/{len(PIPELINE_STEPS)}")
    
    with col2:
        running = sum(1 for s in st.session_state.pipeline_progress.values() if s.status == 'running')
        st.metric("In Progress", running)
    
    with col3:
        failed = sum(1 for s in st.session_state.pipeline_progress.values() if s.status == 'failed')
        st.metric("Failed", failed, delta_color="inverse" if failed > 0 else "off")
    
    # Progress bar
    st.progress(progress / 100)
    
    # Step details
    st.markdown("---")
    
    status_icons = {
        'pending': '⚪',
        'running': '🔵',
        'completed': '✅',
        'failed': '❌',
        'skipped': '⏭️'
    }
    
    cols = st.columns(5)
    for i, (step_id, step_name) in enumerate(PIPELINE_STEPS):
        step = st.session_state.pipeline_progress[step_id]
        icon = status_icons.get(step.status, '⚪')
        
        with cols[i % 5]:
            st.markdown(f"{icon} **{step_name}**")
            if step.message:
                st.caption(step.message)


def render_mini_progress():
    """Render compact progress indicator for sidebar"""
    init_progress_tracker()
    
    progress = get_progress_percentage()
    completed = sum(1 for s in st.session_state.pipeline_progress.values() if s.status == 'completed')
    total = len(PIPELINE_STEPS)
    
    st.caption(f"Progress: {completed}/{total} steps")
    st.progress(progress / 100)


def reset_progress():
    """Reset all progress tracking"""
    if 'pipeline_progress' in st.session_state:
        del st.session_state.pipeline_progress
    if 'pipeline_start_time' in st.session_state:
        del st.session_state.pipeline_start_time
    init_progress_tracker()
