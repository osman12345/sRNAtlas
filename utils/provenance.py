"""
Provenance tracking for sRNAtlas
Records pipeline parameters, versions, and checksums for reproducibility
"""
import hashlib
import yaml
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Any, Union
from dataclasses import dataclass, field, asdict
import subprocess
import platform
import sys


@dataclass
class FileRecord:
    """Record of a file with checksum"""
    path: str
    checksum_md5: str
    size_bytes: int
    modified_time: str


@dataclass
class ToolVersion:
    """Version information for a tool"""
    name: str
    version: str
    path: Optional[str] = None


@dataclass
class StepRecord:
    """Record of a pipeline step"""
    step_name: str
    timestamp: str
    parameters: Dict[str, Any]
    input_files: List[FileRecord]
    output_files: List[FileRecord]
    status: str  # 'success', 'warning', 'error'
    duration_seconds: Optional[float] = None
    notes: Optional[str] = None


@dataclass
class ProvenanceRecord:
    """Complete provenance record for a pipeline run"""
    pipeline_name: str
    pipeline_version: str
    run_id: str
    start_time: str
    end_time: Optional[str]
    status: str
    system_info: Dict[str, str]
    tool_versions: List[ToolVersion]
    steps: List[StepRecord]
    parameters: Dict[str, Any]
    metadata: Dict[str, Any] = field(default_factory=dict)


def calculate_file_checksum(file_path: Union[str, Path], algorithm: str = 'md5') -> str:
    """
    Calculate checksum of a file

    Args:
        file_path: Path to file
        algorithm: Hash algorithm ('md5', 'sha256', etc.)

    Returns:
        Hex digest of file hash
    """
    file_path = Path(file_path)

    if not file_path.exists():
        return "FILE_NOT_FOUND"

    hash_obj = hashlib.new(algorithm)

    with open(file_path, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b''):
            hash_obj.update(chunk)

    return hash_obj.hexdigest()


def get_file_record(file_path: Union[str, Path]) -> FileRecord:
    """
    Create a FileRecord for a file

    Args:
        file_path: Path to file

    Returns:
        FileRecord with metadata
    """
    file_path = Path(file_path)

    if not file_path.exists():
        return FileRecord(
            path=str(file_path),
            checksum_md5="FILE_NOT_FOUND",
            size_bytes=0,
            modified_time="N/A"
        )

    stat = file_path.stat()

    return FileRecord(
        path=str(file_path.absolute()),
        checksum_md5=calculate_file_checksum(file_path),
        size_bytes=stat.st_size,
        modified_time=datetime.fromtimestamp(stat.st_mtime).isoformat()
    )


def get_system_info() -> Dict[str, str]:
    """Get system information for provenance"""
    return {
        'platform': platform.platform(),
        'python_version': platform.python_version(),
        'processor': platform.processor(),
        'hostname': platform.node(),
        'os': platform.system(),
        'os_release': platform.release()
    }


def get_tool_version(tool_name: str) -> ToolVersion:
    """
    Get version of an external tool

    Args:
        tool_name: Name of tool (e.g., 'bowtie', 'samtools')

    Returns:
        ToolVersion with version info
    """
    import shutil

    tool_path = shutil.which(tool_name)

    if not tool_path:
        return ToolVersion(
            name=tool_name,
            version="NOT_INSTALLED",
            path=None
        )

    # Tool-specific version commands
    version_commands = {
        'bowtie': [tool_name, '--version'],
        'samtools': [tool_name, '--version'],
        'cutadapt': [tool_name, '--version'],
        'fastqc': [tool_name, '--version'],
        'miranda': [tool_name, '-h'],  # miRanda doesn't have --version
        'python': [tool_name, '--version'],
    }

    cmd = version_commands.get(tool_name, [tool_name, '--version'])

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
        output = result.stdout.strip() or result.stderr.strip()
        # Get first line which usually contains version
        version_line = output.split('\n')[0].strip() if output else "UNKNOWN"
        if not version_line:
            version_line = "UNKNOWN"
        return ToolVersion(
            name=tool_name,
            version=version_line,
            path=tool_path
        )
    except Exception:
        return ToolVersion(
            name=tool_name,
            version="UNKNOWN",
            path=tool_path
        )


def get_python_package_versions() -> Dict[str, str]:
    """Get versions of key Python packages"""
    packages = [
        'streamlit', 'pandas', 'numpy', 'plotly', 'pysam',
        'biopython', 'pydeseq2', 'scipy', 'scikit-learn'
    ]

    versions = {}
    for pkg in packages:
        try:
            import importlib.metadata
            versions[pkg] = importlib.metadata.version(pkg)
        except Exception:
            versions[pkg] = "NOT_INSTALLED"

    return versions


class ProvenanceTracker:
    """
    Track provenance throughout a pipeline run

    Usage:
        tracker = ProvenanceTracker("sRNAtlas", "1.0.0")
        tracker.start_run()

        tracker.start_step("trimming", {"adapter": "AGATCGGAAGAG"})
        tracker.add_input_file("reads.fastq")
        # ... run step ...
        tracker.add_output_file("trimmed.fastq")
        tracker.end_step("success")

        tracker.end_run()
        tracker.save("provenance.yaml")
    """

    def __init__(
        self,
        pipeline_name: str = "sRNAtlas",
        pipeline_version: str = "1.0.0"
    ):
        self.pipeline_name = pipeline_name
        self.pipeline_version = pipeline_version
        self.record: Optional[ProvenanceRecord] = None
        self.current_step: Optional[Dict] = None
        self._step_start_time: Optional[datetime] = None

    def start_run(
        self,
        run_id: Optional[str] = None,
        parameters: Optional[Dict] = None,
        metadata: Optional[Dict] = None
    ) -> str:
        """
        Start a new pipeline run

        Args:
            run_id: Optional custom run ID (auto-generated if None)
            parameters: Global pipeline parameters
            metadata: Additional metadata (project name, etc.)

        Returns:
            Run ID
        """
        if run_id is None:
            run_id = datetime.now().strftime("%Y%m%d_%H%M%S")

        # Get tool versions
        tools = ['bowtie', 'samtools', 'cutadapt', 'fastqc', 'miranda']
        tool_versions = [get_tool_version(t) for t in tools]

        self.record = ProvenanceRecord(
            pipeline_name=self.pipeline_name,
            pipeline_version=self.pipeline_version,
            run_id=run_id,
            start_time=datetime.now().isoformat(),
            end_time=None,
            status='running',
            system_info=get_system_info(),
            tool_versions=tool_versions,
            steps=[],
            parameters=parameters or {},
            metadata=metadata or {}
        )

        return run_id

    def start_step(
        self,
        step_name: str,
        parameters: Optional[Dict] = None
    ):
        """
        Start a new pipeline step

        Args:
            step_name: Name of the step
            parameters: Step-specific parameters
        """
        if self.record is None:
            raise RuntimeError("Call start_run() first")

        self.current_step = {
            'step_name': step_name,
            'timestamp': datetime.now().isoformat(),
            'parameters': parameters or {},
            'input_files': [],
            'output_files': []
        }
        self._step_start_time = datetime.now()

    def add_input_file(self, file_path: Union[str, Path]):
        """Add an input file to the current step"""
        if self.current_step is None:
            raise RuntimeError("Call start_step() first")

        record = get_file_record(file_path)
        self.current_step['input_files'].append(record)

    def add_output_file(self, file_path: Union[str, Path]):
        """Add an output file to the current step"""
        if self.current_step is None:
            raise RuntimeError("Call start_step() first")

        record = get_file_record(file_path)
        self.current_step['output_files'].append(record)

    def end_step(
        self,
        status: str = 'success',
        notes: Optional[str] = None
    ):
        """
        End the current step

        Args:
            status: 'success', 'warning', or 'error'
            notes: Optional notes about the step
        """
        if self.current_step is None:
            raise RuntimeError("No active step")

        duration = None
        if self._step_start_time:
            duration = (datetime.now() - self._step_start_time).total_seconds()

        step_record = StepRecord(
            step_name=self.current_step['step_name'],
            timestamp=self.current_step['timestamp'],
            parameters=self.current_step['parameters'],
            input_files=self.current_step['input_files'],
            output_files=self.current_step['output_files'],
            status=status,
            duration_seconds=duration,
            notes=notes
        )

        self.record.steps.append(step_record)
        self.current_step = None
        self._step_start_time = None

    def end_run(self, status: str = 'success'):
        """
        End the pipeline run

        Args:
            status: Final status ('success', 'warning', 'error')
        """
        if self.record is None:
            raise RuntimeError("No active run")

        self.record.end_time = datetime.now().isoformat()
        self.record.status = status

    def add_metadata(self, key: str, value: Any):
        """Add metadata to the provenance record"""
        if self.record is None:
            raise RuntimeError("Call start_run() first")
        self.record.metadata[key] = value

    def to_dict(self) -> Dict:
        """Convert provenance record to dictionary"""
        if self.record is None:
            return {}

        # Custom conversion to handle dataclasses
        def convert(obj):
            if hasattr(obj, '__dataclass_fields__'):
                return {k: convert(v) for k, v in asdict(obj).items()}
            elif isinstance(obj, list):
                return [convert(item) for item in obj]
            elif isinstance(obj, dict):
                return {k: convert(v) for k, v in obj.items()}
            else:
                return obj

        return convert(self.record)

    def save(
        self,
        output_path: Union[str, Path],
        format: str = 'yaml'
    ):
        """
        Save provenance record to file

        Args:
            output_path: Path for output file
            format: 'yaml' or 'json'
        """
        output_path = Path(output_path)
        data = self.to_dict()

        if format == 'yaml':
            with open(output_path, 'w') as f:
                yaml.dump(data, f, default_flow_style=False, sort_keys=False)
        else:
            with open(output_path, 'w') as f:
                json.dump(data, f, indent=2)

    @classmethod
    def load(cls, input_path: Union[str, Path]) -> 'ProvenanceTracker':
        """
        Load a provenance record from file

        Args:
            input_path: Path to provenance file

        Returns:
            ProvenanceTracker with loaded record
        """
        input_path = Path(input_path)

        with open(input_path, 'r') as f:
            if input_path.suffix in ['.yaml', '.yml']:
                data = yaml.safe_load(f)
            else:
                data = json.load(f)

        tracker = cls(
            pipeline_name=data.get('pipeline_name', 'unknown'),
            pipeline_version=data.get('pipeline_version', 'unknown')
        )

        # Reconstruct record from data
        # This is a simplified reconstruction
        tracker.record = ProvenanceRecord(
            pipeline_name=data.get('pipeline_name', ''),
            pipeline_version=data.get('pipeline_version', ''),
            run_id=data.get('run_id', ''),
            start_time=data.get('start_time', ''),
            end_time=data.get('end_time'),
            status=data.get('status', ''),
            system_info=data.get('system_info', {}),
            tool_versions=[
                ToolVersion(**tv) for tv in data.get('tool_versions', [])
            ],
            steps=[],  # Would need full reconstruction
            parameters=data.get('parameters', {}),
            metadata=data.get('metadata', {})
        )

        return tracker


def generate_provenance_summary(tracker: ProvenanceTracker) -> str:
    """
    Generate a human-readable provenance summary

    Args:
        tracker: ProvenanceTracker with completed run

    Returns:
        Formatted summary string
    """
    if tracker.record is None:
        return "No provenance record available"

    record = tracker.record

    lines = [
        "=" * 60,
        f"PROVENANCE RECORD: {record.pipeline_name} v{record.pipeline_version}",
        "=" * 60,
        "",
        f"Run ID: {record.run_id}",
        f"Status: {record.status}",
        f"Started: {record.start_time}",
        f"Ended: {record.end_time or 'In progress'}",
        "",
        "SYSTEM INFO:",
        f"  Platform: {record.system_info.get('platform', 'unknown')}",
        f"  Python: {record.system_info.get('python_version', 'unknown')}",
        "",
        "TOOL VERSIONS:"
    ]

    for tool in record.tool_versions:
        status = "✓" if tool.version != "NOT_INSTALLED" else "✗"
        lines.append(f"  {status} {tool.name}: {tool.version}")

    lines.extend(["", "PIPELINE STEPS:"])

    for i, step in enumerate(record.steps, 1):
        status_icon = {"success": "✅", "warning": "⚠️", "error": "❌"}.get(step.status, "?")
        lines.append(f"  {i}. {status_icon} {step.step_name}")
        lines.append(f"     Time: {step.timestamp}")

        if step.duration_seconds:
            lines.append(f"     Duration: {step.duration_seconds:.1f}s")

        if step.input_files:
            lines.append(f"     Inputs: {len(step.input_files)} files")
        if step.output_files:
            lines.append(f"     Outputs: {len(step.output_files)} files")

        if step.notes:
            lines.append(f"     Notes: {step.notes}")

    lines.extend(["", "=" * 60])

    return "\n".join(lines)


# Convenience function for Streamlit integration
def get_session_provenance_tracker() -> ProvenanceTracker:
    """
    Get or create a ProvenanceTracker for the current Streamlit session

    Returns:
        ProvenanceTracker instance
    """
    import streamlit as st

    if 'provenance_tracker' not in st.session_state:
        st.session_state.provenance_tracker = ProvenanceTracker()

    return st.session_state.provenance_tracker


def render_provenance_ui():
    """Render provenance tracking UI in Streamlit"""
    import streamlit as st

    st.subheader("📋 Analysis Provenance")

    tracker = get_session_provenance_tracker()

    if tracker.record is None:
        st.info("No analysis has been run yet. Provenance will be recorded as you run analyses.")
        return

    # Summary
    record = tracker.record
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric("Run ID", record.run_id[:8] + "...")
    with col2:
        st.metric("Steps Completed", len(record.steps))
    with col3:
        status_icon = {"success": "✅", "running": "🔄", "error": "❌"}.get(record.status, "?")
        st.metric("Status", f"{status_icon} {record.status}")
    with col4:
        if record.steps:
            total_time = sum(s.duration_seconds or 0 for s in record.steps)
            st.metric("Total Time", f"{total_time:.1f}s")

    # Steps timeline
    st.markdown("#### Pipeline Steps")

    for step in record.steps:
        status_icon = {"success": "✅", "warning": "⚠️", "error": "❌"}.get(step.status, "?")

        with st.expander(f"{status_icon} {step.step_name}"):
            st.markdown(f"**Timestamp:** {step.timestamp}")

            if step.duration_seconds:
                st.markdown(f"**Duration:** {step.duration_seconds:.1f} seconds")

            if step.parameters:
                st.markdown("**Parameters:**")
                st.json(step.parameters)

            if step.input_files:
                st.markdown(f"**Input Files:** ({len(step.input_files)})")
                for f in step.input_files[:5]:
                    st.caption(f"📄 {Path(f.path).name} ({f.size_bytes:,} bytes)")
                    st.caption(f"   MD5: {f.checksum_md5[:16]}...")

            if step.output_files:
                st.markdown(f"**Output Files:** ({len(step.output_files)})")
                for f in step.output_files[:5]:
                    st.caption(f"📄 {Path(f.path).name} ({f.size_bytes:,} bytes)")

            if step.notes:
                st.info(step.notes)

    # Download options
    st.divider()
    col1, col2 = st.columns(2)

    with col1:
        yaml_content = yaml.dump(tracker.to_dict(), default_flow_style=False)
        st.download_button(
            "📥 Download Provenance (YAML)",
            yaml_content,
            f"provenance_{record.run_id}.yaml",
            "text/yaml"
        )

    with col2:
        json_content = json.dumps(tracker.to_dict(), indent=2)
        st.download_button(
            "📥 Download Provenance (JSON)",
            json_content,
            f"provenance_{record.run_id}.json",
            "application/json"
        )

    # Summary view
    with st.expander("📄 View Full Summary"):
        st.text(generate_provenance_summary(tracker))
