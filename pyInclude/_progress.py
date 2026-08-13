# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Private terminal progress rendering for compiled HASE simulations."""

from __future__ import annotations

import io
import os
import time


def _human_duration(seconds):
    total_seconds = max(0, int(seconds))
    hours, remainder = divmod(total_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    parts = []
    if hours:
        parts.append(f"{hours}h")
    if minutes:
        parts.append(f"{minutes}m")
    if seconds or not parts:
        parts.append(f"{seconds}s")
    return " ".join(parts)


class _ProgressBar:
    """Render backend step completion on one persistent terminal line."""

    _symbols = (
        "ø",
        "¤",
        "º",
        "°",
        "`",
        "°",
        "°",
        "¤",
        "ø",
        ",",
        "¸",
        ",",
    )

    def __init__(self, stream, *, clock=time.monotonic, terminal_columns=None):
        self._stream = stream
        self._clock = clock
        self._terminal_columns = terminal_columns
        self._buffer = io.StringIO()
        self._started = None
        self._tic = 0
        self._max_total = 0
        self._last_line = ""
        self._line_width = 0
        self._drawn = False
        self._complete = False

    def _bar(self, progress, length):
        if os.name == "nt":
            return "#" * progress + " " * (length - progress)
        return "".join(self._symbols[(self._tic + index) % len(self._symbols)] for index in range(progress)) + (
            " " * (length - progress)
        )

    def _available_columns(self):
        if self._terminal_columns is not None:
            return int(self._terminal_columns)
        try:
            return os.get_terminal_size(self._stream.fileno()).columns
        except (AttributeError, OSError, io.UnsupportedOperation):
            return None

    def _format(self, current, total, elapsed):
        percentage = current / total
        estimated_total = elapsed / percentage if percentage > 0.0 else 0.0
        remaining = max(0.0, estimated_total - elapsed)
        width = len(str(total))

        status = (
            f"] {int(percentage * 100):3d}%"
            f" ({current:{width}d}/{total})"
            f" after {_human_duration(elapsed)}"
            f" ({_human_duration(estimated_total)} total, {_human_duration(remaining)} remaining)"
        )
        maximum_length = 16 if os.name == "nt" else 50
        columns = self._available_columns()
        if columns is None:
            length = maximum_length
        else:
            fixed_width = len("[PROGRESS] [") + len(status)
            length = max(1, min(maximum_length, columns - fixed_width - 1))
        progress = min(length, int(percentage * length))

        self._buffer.seek(0)
        self._buffer.truncate()
        self._buffer.write("[PROGRESS] [")
        self._buffer.write(self._bar(progress, length))
        self._buffer.write(status)
        return self._buffer.getvalue()

    def update(self, current, total):
        current = int(current)
        total = int(total)
        if total <= 0:
            return
        if self._started is None:
            self._started = self._clock()
        self._max_total = max(self._max_total, total)
        total = self._max_total
        current = min(current, total)
        elapsed = max(0.0, self._clock() - self._started)
        if self._tic and elapsed <= 0.035 * self._tic and current != total:
            return

        self._tic += 1
        self._last_line = self._format(current, total, elapsed)
        self._line_width = max(self._line_width, len(self._last_line))
        self._stream.write("\r" + self._last_line.ljust(self._line_width))
        if current == total:
            self._stream.write("\n")
            self._drawn = False
            self._complete = True
        else:
            self._drawn = True
        self._stream.flush()

    def clear(self):
        if not self._drawn:
            return
        self._stream.write("\r" + " " * self._line_width + "\r")
        self._stream.flush()
        self._drawn = False

    def redraw(self):
        if not self._last_line or self._complete:
            return
        self._stream.write("\r" + self._last_line.ljust(self._line_width))
        self._stream.flush()
        self._drawn = True

    def finish(self):
        if not self._drawn:
            return
        self._stream.write("\n")
        self._stream.flush()
        self._drawn = False
