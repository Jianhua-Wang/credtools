"""Web pages for CREDTOOLS visualization."""

# Import pages to register them with Dash
try:
    from . import home, locus
except ImportError:
    # Pages require web dependencies which may not be installed
    pass

__all__ = ["home", "locus"]
