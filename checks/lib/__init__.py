"""Shared helpers for the verification checks."""
import os


class SkipCheck(Exception):
    """Raised by a check that is intentionally not run in the current mode."""


def skip(reason):
    """Skip the current check.

    Under pytest this defers to pytest's own skip; under run_all.py it raises
    SkipCheck, which the runner reports as SKIP.
    """
    if "PYTEST_CURRENT_TEST" in os.environ:
        import pytest

        pytest.skip(reason)
    raise SkipCheck(reason)
