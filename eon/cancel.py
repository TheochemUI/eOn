"""Cooperative cancel for in-process jobs.

The token is polled before each job and before each in-process
dispatch. A force-call inside a compiled relax still runs to the next
poll. ``cancel_state`` sets the flag.
"""

from __future__ import annotations

import threading


class Cancelled(Exception):
    """Raised when a job sees a cancelled token."""


class CancelToken:
    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._cancelled = False

    def cancel(self) -> None:
        with self._lock:
            self._cancelled = True

    @property
    def cancelled(self) -> bool:
        with self._lock:
            return self._cancelled

    def raise_if_cancelled(self) -> None:
        if self.cancelled:
            raise Cancelled("job cancelled")
