"""
Tests for progress tracker utility
"""
import pytest
import sys
from pathlib import Path

# Add parent to path for direct module import
sys.path.insert(0, str(Path(__file__).parent.parent / "utils"))


class TestProgressTracker:
    """Tests for pipeline progress tracking"""
    
    def test_pipeline_steps_defined(self):
        """Test that pipeline steps are defined"""
        from progress_tracker import PIPELINE_STEPS
        
        assert len(PIPELINE_STEPS) > 0
        assert all(isinstance(s, tuple) and len(s) == 2 for s in PIPELINE_STEPS)
    
    def test_pipeline_step_class(self):
        """Test PipelineStep dataclass"""
        from progress_tracker import PipelineStep
        
        step = PipelineStep(name="Test Step", status="pending")
        assert step.name == "Test Step"
        assert step.status == "pending"
        assert step.message == ""
    
    def test_pipeline_step_with_message(self):
        """Test PipelineStep with message"""
        from progress_tracker import PipelineStep
        
        step = PipelineStep(name="QC", status="completed", message="Done!")
        assert step.status == "completed"
        assert step.message == "Done!"
    
    def test_pipeline_step_statuses(self):
        """Test valid status values"""
        from progress_tracker import PipelineStep
        
        valid_statuses = ['pending', 'running', 'completed', 'failed', 'skipped']
        for status in valid_statuses:
            step = PipelineStep(name="Test", status=status)
            assert step.status == status
    
    def test_pipeline_steps_have_required_steps(self):
        """Test that key pipeline steps are defined"""
        from progress_tracker import PIPELINE_STEPS
        
        step_ids = [s[0] for s in PIPELINE_STEPS]
        
        required = ['project', 'qc', 'trimming', 'alignment', 'counting']
        for req in required:
            assert req in step_ids, f"Missing required step: {req}"
