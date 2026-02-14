"""
Standardized Result type for sRNAtlas
Provides consistent success/error returns across all modules
"""
from typing import Any, Dict, Optional, TypeVar, Generic
from dataclasses import dataclass, field


T = TypeVar('T')


@dataclass
class Result:
    """
    Standardized result container for function returns.

    Usage:
        # Success
        return Result.ok(data={'counts': count_dict}, message='Counted 1000 reads')

        # Error
        return Result.error('File not found: input.bam', category='file_error')

        # Check result
        result = some_function()
        if result.success:
            process(result.data)
        else:
            show_error(result.error_message)
    """
    success: bool
    data: Optional[Dict[str, Any]] = field(default_factory=dict)
    error_message: Optional[str] = None
    error_category: Optional[str] = None
    warnings: list = field(default_factory=list)
    message: Optional[str] = None

    @classmethod
    def ok(cls, data: Optional[Dict[str, Any]] = None, message: str = '', warnings: list = None) -> 'Result':
        """Create a success result"""
        return cls(
            success=True,
            data=data or {},
            message=message,
            warnings=warnings or []
        )

    @classmethod
    def error(cls, error_message: str, category: str = 'unknown', data: Optional[Dict] = None) -> 'Result':
        """Create an error result"""
        return cls(
            success=False,
            error_message=error_message,
            error_category=category,
            data=data or {}
        )

    def add_warning(self, warning: str):
        """Add a warning to the result"""
        self.warnings.append(warning)

    def to_dict(self) -> Dict:
        """Convert to legacy dict format for backward compatibility"""
        if self.success:
            result = {'status': 'success'}
            result.update(self.data)
            if self.message:
                result['message'] = self.message
            if self.warnings:
                result['warnings'] = self.warnings
        else:
            result = {
                'status': 'error',
                'error': self.error_message,
                'category': self.error_category
            }
        return result

    @classmethod
    def from_dict(cls, d: Dict) -> 'Result':
        """Create Result from legacy dict format"""
        if d.get('status') == 'success':
            return cls.ok(
                data={k: v for k, v in d.items() if k != 'status'},
                message=d.get('message', '')
            )
        else:
            return cls.error(
                error_message=d.get('error', 'Unknown error'),
                category=d.get('category', 'unknown')
            )

    def __bool__(self):
        """Allow using Result in boolean context"""
        return self.success
