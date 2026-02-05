# Contributing to sRNAtlas

Thank you for your interest in contributing to sRNAtlas! This document provides guidelines and information for contributors.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [How to Contribute](#how-to-contribute)
- [Development Setup](#development-setup)
- [Coding Standards](#coding-standards)
- [Testing](#testing)
- [Pull Request Process](#pull-request-process)

---

## Code of Conduct

We are committed to providing a welcoming and inclusive environment. Please:

- Be respectful and considerate
- Welcome newcomers and help them get started
- Focus on constructive feedback
- Accept different viewpoints gracefully

---

## Getting Started

### Prerequisites

- Python 3.9+
- Git
- Conda (recommended) or pip

### Fork and Clone

1. Fork the repository on GitHub
2. Clone your fork:
   ```bash
   git clone https://github.com/osman12345/sRNAtlas.git
   cd sRNAtlas
   ```
3. Add upstream remote:
   ```bash
   git remote add upstream https://github.com/osman12345/sRNAtlas.git
   ```

---

## How to Contribute

### Reporting Bugs

1. Check existing issues to avoid duplicates
2. Use the bug report template
3. Include:
   - Clear description of the bug
   - Steps to reproduce
   - Expected vs actual behavior
   - System information (OS, Python version)
   - Error messages or screenshots

### Suggesting Features

1. Check existing issues and discussions
2. Use the feature request template
3. Describe the use case and benefits
4. Consider implementation complexity

### Contributing Code

1. Find an issue to work on (or create one)
2. Comment on the issue to claim it
3. Create a feature branch
4. Write code and tests
5. Submit a pull request

---

## Development Setup

### Create Development Environment

```bash
# Clone repository
git clone https://github.com/osman12345/sRNAtlas.git
cd sRNAtlas

# Create conda environment
conda create -n srnatlas-dev python=3.10 -y
conda activate srnatlas-dev

# Install dependencies
pip install -r requirements.txt

# Install development tools
pip install pytest pytest-cov black flake8 mypy

# Install bioinformatics tools
conda install -c bioconda bowtie samtools -y
```

### Run the Application

```bash
streamlit run app/main.py
```

### Run Tests

```bash
# All tests
python -m pytest

# With coverage
python -m pytest --cov=utils --cov=modules

# Specific file
python -m pytest tests/test_file_handlers.py -v
```

---

## Coding Standards

### Python Style

- Follow PEP 8 guidelines
- Use Black for code formatting
- Maximum line length: 100 characters
- Use meaningful variable and function names

### Code Formatting

```bash
# Format code with Black
black .

# Check style with flake8
flake8 .
```

### Type Hints

Use type hints for function signatures:

```python
def process_reads(
    fastq_file: Path,
    min_length: int = 18,
    max_length: int = 35
) -> pd.DataFrame:
    """Process FASTQ reads and return statistics."""
    ...
```

### Docstrings

Use Google-style docstrings:

```python
def align_reads(fastq_file: Path, index_prefix: Path) -> dict:
    """
    Align reads to reference using Bowtie.
    
    Args:
        fastq_file: Path to input FASTQ file
        index_prefix: Path prefix for Bowtie index
        
    Returns:
        Dictionary containing alignment statistics
        
    Raises:
        FileNotFoundError: If input file doesn't exist
        subprocess.CalledProcessError: If alignment fails
    """
    ...
```

### Imports

Organize imports in this order:
1. Standard library
2. Third-party packages
3. Local modules

```python
import os
import sys
from pathlib import Path

import pandas as pd
import numpy as np
import streamlit as st

from utils.file_handlers import read_fastq_stats
from config import settings
```

---

## Testing

### Writing Tests

- Place tests in the `tests/` directory
- Name test files `test_*.py`
- Name test functions `test_*`
- Use pytest fixtures for common setup

```python
import pytest
from utils.file_handlers import read_fastq_stats

class TestFastqReader:
    def test_read_valid_file(self, temp_fastq_file):
        """Test reading a valid FASTQ file."""
        result = read_fastq_stats(temp_fastq_file)
        assert result['total_reads'] > 0
        
    def test_read_nonexistent_file(self):
        """Test handling of nonexistent file."""
        with pytest.raises(FileNotFoundError):
            read_fastq_stats(Path('/nonexistent.fastq'))
```

### Test Coverage

Aim for at least 80% code coverage for new features.

---

## Pull Request Process

### Before Submitting

1. **Sync with upstream:**
   ```bash
   git fetch upstream
   git rebase upstream/main
   ```

2. **Run tests:**
   ```bash
   python -m pytest
   ```

3. **Check code style:**
   ```bash
   black --check .
   flake8 .
   ```

4. **Update documentation** if needed

### PR Guidelines

- Create a descriptive title
- Reference related issues
- Describe changes made
- Include screenshots for UI changes
- Ensure all tests pass

### PR Template

```markdown
## Description
Brief description of changes

## Related Issue
Fixes #123

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Documentation update
- [ ] Refactoring

## Checklist
- [ ] Tests pass locally
- [ ] Code follows style guidelines
- [ ] Documentation updated
- [ ] No new warnings
```

### Review Process

1. Maintainers will review your PR
2. Address any requested changes
3. Once approved, PR will be merged

---

## Project Structure

```
sRNAtlas/
├── app/               # Main application
├── modules/           # Analysis modules
├── utils/             # Utility functions
├── config/            # Configuration
├── tests/             # Test files
├── docs/              # Documentation
└── assets/            # Static assets
```

### Adding a New Module

1. Create module file in `modules/`
2. Add to navigation in `app/main.py`
3. Create tests in `tests/`
4. Update documentation

---

## Questions?

- Open a GitHub Discussion
- Check existing issues
- Review documentation

Thank you for contributing to sRNAtlas!
