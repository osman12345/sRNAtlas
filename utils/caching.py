"""
Caching utilities for sRNAtlas
Provides helper functions for Streamlit caching with complex types
"""
import hashlib
from pathlib import Path
from typing import Any, List, Dict, Optional, Tuple, Union
import pandas as pd
import numpy as np


def hash_dataframe(df: pd.DataFrame) -> str:
    """Create a hash string from a DataFrame for caching purposes"""
    if df is None:
        return "none"
    # Use shape + column names + sample of data for efficient hashing
    hash_str = f"{df.shape}_{list(df.columns)}_{df.index.tolist()[:10]}"
    # Add sample of actual values for collision avoidance
    if len(df) > 0:
        sample = df.head(5).to_string()
        hash_str += f"_{sample}"
    return hashlib.md5(hash_str.encode()).hexdigest()


def hash_file(file_path: Union[str, Path]) -> str:
    """Create a hash based on file path and modification time"""
    path = Path(file_path)
    if not path.exists():
        return f"missing_{str(path)}"
    # Use path + mtime + size for efficient change detection
    stat = path.stat()
    return f"{path.absolute()}_{stat.st_mtime}_{stat.st_size}"


def hash_list(items: List[Any]) -> str:
    """Convert a list to a hashable tuple string"""
    if items is None:
        return "none"
    return str(tuple(sorted(str(x) for x in items)))


def hash_dict(d: Dict) -> str:
    """Create a hash from a dictionary"""
    if d is None:
        return "none"
    # Sort keys for consistent hashing
    sorted_items = sorted((str(k), str(v)) for k, v in d.items())
    return hashlib.md5(str(sorted_items).encode()).hexdigest()


def make_hashable(obj: Any) -> Any:
    """
    Convert an object to a hashable form for Streamlit caching.

    Handles:
    - Path -> str
    - DataFrame -> hash string
    - List -> tuple
    - Dict -> frozenset of items
    - None -> "none"
    - Other -> str
    """
    if obj is None:
        return "none"
    elif isinstance(obj, Path):
        return hash_file(obj)
    elif isinstance(obj, pd.DataFrame):
        return hash_dataframe(obj)
    elif isinstance(obj, list):
        return tuple(make_hashable(x) for x in obj)
    elif isinstance(obj, dict):
        return tuple(sorted((k, make_hashable(v)) for k, v in obj.items()))
    elif isinstance(obj, (str, int, float, bool)):
        return obj
    else:
        return str(obj)
