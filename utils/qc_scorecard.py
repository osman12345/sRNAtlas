"""
QC Scorecard utilities for sRNAtlas
Provides summary metrics and traffic-light flags for quality assessment
"""
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from enum import Enum


class QCStatus(Enum):
    """QC status levels"""
    OK = "ok"
    WARNING = "warning"
    CRITICAL = "critical"
    UNKNOWN = "unknown"


@dataclass
class QCMetric:
    """A single QC metric with thresholds"""
    name: str
    value: float
    unit: str
    status: QCStatus
    description: str
    threshold_warning: Optional[float] = None
    threshold_critical: Optional[float] = None


@dataclass
class SampleScorecard:
    """QC scorecard for a single sample"""
    sample_name: str
    metrics: List[QCMetric]
    overall_status: QCStatus
    summary: str
    recommendations: List[str]


# Default thresholds for small RNA-seq QC
QC_THRESHOLDS = {
    'total_reads': {
        'warning': 1_000_000,  # Less than 1M reads
        'critical': 100_000,   # Less than 100k reads
        'direction': 'min'     # Higher is better
    },
    'adapter_pct': {
        'warning': 30,         # More than 30% adapters
        'critical': 60,        # More than 60% adapters
        'direction': 'max'     # Lower is better
    },
    'trimming_pass_rate': {
        'warning': 70,         # Less than 70% pass
        'critical': 50,        # Less than 50% pass
        'direction': 'min'     # Higher is better
    },
    'alignment_rate': {
        'warning': 50,         # Less than 50% alignment
        'critical': 20,        # Less than 20% alignment
        'direction': 'min'     # Higher is better
    },
    'rrna_pct': {
        'warning': 20,         # More than 20% rRNA
        'critical': 50,        # More than 50% rRNA
        'direction': 'max'     # Lower is better
    },
    'mean_quality': {
        'warning': 25,         # Less than Q25
        'critical': 20,        # Less than Q20
        'direction': 'min'     # Higher is better
    },
    'mirna_fraction': {
        'warning': 30,         # Less than 30% miRNA
        'critical': 10,        # Less than 10% miRNA
        'direction': 'min'     # Higher is better (for miRNA studies)
    }
}


def evaluate_metric(
    name: str,
    value: float,
    unit: str = "",
    description: str = "",
    thresholds: Optional[Dict] = None
) -> QCMetric:
    """
    Evaluate a single metric against thresholds

    Args:
        name: Metric name
        value: Metric value
        unit: Unit of measurement
        description: Description of the metric
        thresholds: Custom thresholds (optional, uses defaults if not provided)

    Returns:
        QCMetric with status
    """
    if thresholds is None:
        thresholds = QC_THRESHOLDS.get(name, {})

    threshold_warning = thresholds.get('warning')
    threshold_critical = thresholds.get('critical')
    direction = thresholds.get('direction', 'min')

    # Determine status
    if value is None or np.isnan(value):
        status = QCStatus.UNKNOWN
    elif direction == 'min':
        # Higher is better
        if threshold_critical and value < threshold_critical:
            status = QCStatus.CRITICAL
        elif threshold_warning and value < threshold_warning:
            status = QCStatus.WARNING
        else:
            status = QCStatus.OK
    else:
        # Lower is better (direction == 'max')
        if threshold_critical and value > threshold_critical:
            status = QCStatus.CRITICAL
        elif threshold_warning and value > threshold_warning:
            status = QCStatus.WARNING
        else:
            status = QCStatus.OK

    return QCMetric(
        name=name,
        value=value,
        unit=unit,
        status=status,
        description=description,
        threshold_warning=threshold_warning,
        threshold_critical=threshold_critical
    )


def classify_size_distribution(
    lengths: List[int],
    counts: Optional[Dict[int, int]] = None
) -> Tuple[str, str]:
    """
    Classify the size distribution pattern

    Args:
        lengths: List of read lengths
        counts: Optional pre-computed length counts

    Returns:
        Tuple of (pattern_name, description)
    """
    if not lengths and not counts:
        return "unknown", "No length data available"

    if counts is None:
        from collections import Counter
        counts = Counter(lengths)

    total = sum(counts.values())
    if total == 0:
        return "unknown", "No reads"

    # Calculate fractions for key ranges
    mirna_range = sum(counts.get(i, 0) for i in range(20, 25))  # 20-24 nt
    sirna_range = sum(counts.get(i, 0) for i in range(21, 25))  # 21-24 nt
    pirna_range = sum(counts.get(i, 0) for i in range(26, 32))  # 26-31 nt
    trf_range = sum(counts.get(i, 0) for i in range(30, 41))    # 30-40 nt
    hc_sirna_range = sum(counts.get(24, 0) for i in [24])       # 24 nt (plant hc-siRNA)

    mirna_pct = 100 * mirna_range / total
    pirna_pct = 100 * pirna_range / total
    trf_pct = 100 * trf_range / total
    hc_sirna_pct = 100 * hc_sirna_range / total

    # Classify based on dominant pattern
    patterns = [
        (mirna_pct, "miRNA-dominant", f"Strong 20-24 nt peak ({mirna_pct:.1f}%)"),
        (pirna_pct, "piRNA-dominant", f"Strong 26-31 nt peak ({pirna_pct:.1f}%)"),
        (trf_pct, "tRF-dominant", f"Strong 30-40 nt peak ({trf_pct:.1f}%)"),
    ]

    # Find dominant pattern
    best_pattern = max(patterns, key=lambda x: x[0])

    if best_pattern[0] > 40:
        return best_pattern[1], best_pattern[2]
    elif best_pattern[0] < 20:
        return "degraded", "No clear peak - possibly degraded sample"
    else:
        return "mixed", f"Mixed pattern (miRNA: {mirna_pct:.1f}%, piRNA: {pirna_pct:.1f}%)"


def classify_plant_sirna(counts: Dict[int, int]) -> Tuple[str, str]:
    """
    For plant data: classify 21-nt vs 24-nt dominance

    Args:
        counts: Length counts dictionary

    Returns:
        Tuple of (pattern, implications)
    """
    total = sum(counts.values())
    if total == 0:
        return "unknown", "No data"

    nt21 = counts.get(21, 0)
    nt24 = counts.get(24, 0)

    nt21_pct = 100 * nt21 / total
    nt24_pct = 100 * nt24 / total

    if nt21_pct > nt24_pct * 1.5:
        return "21-nt dominant", f"miRNA/ta-siRNA enriched ({nt21_pct:.1f}% at 21 nt)"
    elif nt24_pct > nt21_pct * 1.5:
        return "24-nt dominant", f"hc-siRNA/phasiRNA enriched ({nt24_pct:.1f}% at 24 nt)"
    else:
        return "balanced", f"Balanced 21/24 nt ({nt21_pct:.1f}% / {nt24_pct:.1f}%)"


def generate_sample_scorecard(
    sample_name: str,
    qc_results: Dict,
    trim_results: Optional[Dict] = None,
    alignment_results: Optional[Dict] = None,
    contamination_results: Optional[Dict] = None,
    is_plant: bool = False
) -> SampleScorecard:
    """
    Generate a comprehensive QC scorecard for a sample

    Args:
        sample_name: Sample identifier
        qc_results: Basic QC results from analyze_fastq_basic
        trim_results: Trimming results (optional)
        alignment_results: Alignment results (optional)
        contamination_results: Contamination check results (optional)
        is_plant: Whether this is plant data (affects 21/24 nt classification)

    Returns:
        SampleScorecard with all metrics and recommendations
    """
    metrics = []
    recommendations = []

    # === Read Count ===
    total_reads = qc_results.get('total_reads', 0)
    metrics.append(evaluate_metric(
        'total_reads',
        total_reads,
        'reads',
        'Total sequenced reads'
    ))

    if total_reads < 1_000_000:
        recommendations.append("Consider deeper sequencing for better coverage")

    # === Read Quality ===
    mean_quality = qc_results.get('mean_quality', 0)
    metrics.append(evaluate_metric(
        'mean_quality',
        mean_quality,
        'Phred',
        'Mean quality score'
    ))

    if mean_quality < 25:
        recommendations.append("Quality scores are low - check sequencing run quality")

    # === Size Distribution ===
    lengths = qc_results.get('lengths', [])
    length_dist = qc_results.get('length_distribution', {})

    if lengths or length_dist:
        pattern, pattern_desc = classify_size_distribution(lengths, length_dist)
        metrics.append(QCMetric(
            name='size_pattern',
            value=0,  # Not numeric
            unit=pattern,
            status=QCStatus.OK if pattern != 'degraded' else QCStatus.WARNING,
            description=pattern_desc
        ))

        if pattern == 'degraded':
            recommendations.append("Sample may be degraded - check RNA integrity")

        # Plant-specific classification
        if is_plant and length_dist:
            plant_pattern, plant_desc = classify_plant_sirna(length_dist)
            metrics.append(QCMetric(
                name='plant_sirna_pattern',
                value=0,
                unit=plant_pattern,
                status=QCStatus.OK,
                description=plant_desc
            ))

    # === Adapter Content (from contamination check) ===
    if contamination_results and contamination_results.get('status') == 'success':
        adapters_detected = contamination_results.get('adapters_detected', [])
        adapter_pct = max(contamination_results.get('adapter_content', {}).values() or [0])

        metrics.append(evaluate_metric(
            'adapter_pct',
            adapter_pct,
            '%',
            f"Adapter content ({', '.join(adapters_detected) if adapters_detected else 'none detected'})"
        ))

        if adapter_pct > 30:
            recommendations.append("High adapter content - ensure trimming is applied")

    # === Trimming Results ===
    if trim_results:
        pass_rate = trim_results.get('pass_rate', 100)
        metrics.append(evaluate_metric(
            'trimming_pass_rate',
            pass_rate,
            '%',
            'Reads passing trimming filters'
        ))

        if pass_rate < 70:
            recommendations.append("Many reads removed during trimming - check adapter selection")

    # === Alignment Results ===
    if alignment_results and alignment_results.get('status') == 'success':
        stats = alignment_results.get('stats', {})
        alignment_rate = stats.get('alignment_rate', 0)

        metrics.append(evaluate_metric(
            'alignment_rate',
            alignment_rate,
            '%',
            'Reads aligned to reference'
        ))

        if alignment_rate < 50:
            recommendations.append("Low alignment rate - verify reference organism or check adapters")

    # === rRNA Contamination ===
    if contamination_results and 'rrna_pct' in contamination_results:
        rrna_pct = contamination_results.get('rrna_pct', 0)
        metrics.append(evaluate_metric(
            'rrna_pct',
            rrna_pct,
            '%',
            'rRNA contamination'
        ))

        if rrna_pct > 20:
            recommendations.append("High rRNA contamination - consider rRNA depletion")

    # === Determine Overall Status ===
    statuses = [m.status for m in metrics]
    if QCStatus.CRITICAL in statuses:
        overall_status = QCStatus.CRITICAL
        summary = "Sample has critical QC issues that need attention"
    elif QCStatus.WARNING in statuses:
        overall_status = QCStatus.WARNING
        summary = "Sample has some QC warnings - review recommended"
    elif QCStatus.UNKNOWN in statuses:
        overall_status = QCStatus.OK
        summary = "Sample passes QC (some metrics unavailable)"
    else:
        overall_status = QCStatus.OK
        summary = "Sample passes all QC checks"

    if not recommendations:
        recommendations.append("No issues detected - sample looks good!")

    return SampleScorecard(
        sample_name=sample_name,
        metrics=metrics,
        overall_status=overall_status,
        summary=summary,
        recommendations=recommendations
    )


def scorecard_to_dataframe(scorecards: List[SampleScorecard]) -> pd.DataFrame:
    """
    Convert a list of scorecards to a summary DataFrame

    Args:
        scorecards: List of SampleScorecard objects

    Returns:
        DataFrame with samples as rows and metrics as columns
    """
    rows = []

    for sc in scorecards:
        row = {
            'Sample': sc.sample_name,
            'Status': sc.overall_status.value.upper(),
            'Summary': sc.summary
        }

        for metric in sc.metrics:
            if metric.unit and metric.unit not in ['%', 'reads', 'Phred']:
                # Non-numeric metric (pattern)
                row[metric.name] = metric.unit
            else:
                row[metric.name] = f"{metric.value:.1f}" if isinstance(metric.value, float) else metric.value

            row[f"{metric.name}_status"] = metric.status.value

        rows.append(row)

    return pd.DataFrame(rows)


def get_status_emoji(status: QCStatus) -> str:
    """Get emoji for QC status"""
    return {
        QCStatus.OK: "✅",
        QCStatus.WARNING: "⚠️",
        QCStatus.CRITICAL: "❌",
        QCStatus.UNKNOWN: "❓"
    }.get(status, "❓")


def get_status_color(status: QCStatus) -> str:
    """Get color for QC status"""
    return {
        QCStatus.OK: "green",
        QCStatus.WARNING: "orange",
        QCStatus.CRITICAL: "red",
        QCStatus.UNKNOWN: "gray"
    }.get(status, "gray")
