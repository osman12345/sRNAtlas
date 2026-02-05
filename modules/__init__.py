"""
Modules package for sRNAtlas
"""
from .qc_module import render_qc_page
from .alignment_module import render_alignment_page
from .counting_module import render_counting_page
from .de_module import render_de_page
from .enrichment_module import render_enrichment_page
from .reports_module import render_reports_page
from .settings_module import render_settings_page
from .help_module import render_help_page

__all__ = [
    'render_qc_page',
    'render_alignment_page',
    'render_counting_page',
    'render_de_page',
    'render_enrichment_page',
    'render_reports_page',
    'render_settings_page',
    'render_help_page',
]
