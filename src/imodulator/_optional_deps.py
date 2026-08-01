"""Helpers for optional dependencies (femwell, solcore).

These packages are declared as pip extras rather than hard requirements, so a
plain ``import`` of imodulator (or of a simulator that does not need them)
succeeds without them installed. Call :func:`require` at the point of use to
turn a missing package into an actionable message.
"""

import importlib


def require(package, extra):
    """Import ``package``, or raise with the install command for its extra.

    Returns the imported module so callers can also use it directly.
    """
    try:
        return importlib.import_module(package)
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            f"'{package}' is required for this feature but is not installed. "
            f"Install it with:  pip install imodulator[{extra}]"
        ) from exc
