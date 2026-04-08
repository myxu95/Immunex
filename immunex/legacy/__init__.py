"""
Legacy compatibility surface for retired Immunex interfaces.

This package is intentionally outside the default public API. Import from here
only when maintaining older integrations that have not yet migrated to the
manifest-based discovery and batch execution contracts.
"""

from .task_discovery import TaskDiscovery
from .batch_processor import BatchProcessor

__all__ = [
    "TaskDiscovery",
    "BatchProcessor",
]
