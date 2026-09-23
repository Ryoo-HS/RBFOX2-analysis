"""rbpc: a rule-based peak caller for CLIP-seq / RBP binding data."""
__version__ = "0.1.0"
from .peaks import Peak  # noqa: E402,F401
__all__ = ["Peak", "__version__"]
