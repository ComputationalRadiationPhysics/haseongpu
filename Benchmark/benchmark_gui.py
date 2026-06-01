#!/usr/bin/env python3
# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Small GUI for plotting HASE benchmark CSV files."""

from __future__ import annotations

import argparse
import csv
import statistics
import tkinter as tk
from collections import defaultdict
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure


NO_LABEL = "<no label>"
DEFAULT_Y_COLUMN = "DurationNs"
TIMESTAMP_ALIASES = {
    "Timestamp": "Timestamp",
    "GlobalTimestampUnixNs": "Timestamp",
    "EventUnixNs": "Timestamp",
    "RunUnixNs": "Timestamp",
}


def normalize_header(name: str) -> str:
    return TIMESTAMP_ALIASES.get(name, name)


def natural_key(value: str) -> tuple[int, object]:
    try:
        return (0, float(value))
    except ValueError:
        return (1, value)


class BenchmarkData:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.rows: list[dict[str, str]] = []
        self.columns: list[str] = []
        self.reload()

    def reload(self) -> None:
        if not self.path.exists():
            raise FileNotFoundError(self.path)

        with self.path.open(newline="") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames is None:
                raise ValueError(f"{self.path} has no CSV header")

            fieldnames = [normalize_header(name) for name in reader.fieldnames]
            rows: list[dict[str, str]] = []
            for raw_row in reader:
                row: dict[str, str] = {}
                for old_name, new_name in zip(reader.fieldnames, fieldnames):
                    row[new_name] = raw_row.get(old_name, "")
                rows.append(row)

        if DEFAULT_Y_COLUMN not in fieldnames:
            raise ValueError(f"{self.path} does not contain required column {DEFAULT_Y_COLUMN}")

        self.columns = list(dict.fromkeys(fieldnames))
        self.rows = rows

    def plottable_columns(self) -> list[str]:
        excluded = {DEFAULT_Y_COLUMN, "EndUnixNs"}
        return [column for column in self.columns if column not in excluded]

    def grouped_average(self, x_column: str, label_column: str | None) -> dict[str, dict[str, float]]:
        buckets: dict[tuple[str, str], list[float]] = defaultdict(list)

        for row in self.rows:
            try:
                duration = float(row[DEFAULT_Y_COLUMN])
            except (TypeError, ValueError):
                continue

            x_value = row.get(x_column, "")
            label = row.get(label_column, "") if label_column else NO_LABEL
            buckets[(label, x_value)].append(duration)

        series: dict[str, dict[str, float]] = defaultdict(dict)
        for (label, x_value), values in buckets.items():
            series[label][x_value] = statistics.fmean(values)
        return dict(series)


class BenchmarkGui(tk.Tk):
    def __init__(self, initial_csv: Path | None) -> None:
        super().__init__()
        self.title("HASE Benchmark Viewer")
        self.geometry("1100x760")

        self.data: BenchmarkData | None = None
        self.csv_path = tk.StringVar(value=str(initial_csv) if initial_csv else "")
        self.x_column = tk.StringVar()
        self.label_column = tk.StringVar(value=NO_LABEL)
        self.unit = tk.StringVar(value="ms")
        self.status = tk.StringVar(value="Select a benchmark CSV file.")
        self.label_filter_vars: dict[str, tk.BooleanVar] = {}
        self.label_filter_menu: tk.Menu | None = None
        self.label_filter_button: ttk.Menubutton | None = None
        self._updating_label_filter = False

        self._build_controls()
        self._build_plot()

        if initial_csv is not None:
            self.load_csv(initial_csv)

    def _build_controls(self) -> None:
        frame = ttk.Frame(self, padding=8)
        frame.pack(side=tk.TOP, fill=tk.X)

        ttk.Label(frame, text="CSV").grid(row=0, column=0, sticky=tk.W)
        csv_entry = ttk.Entry(frame, textvariable=self.csv_path)
        csv_entry.grid(row=0, column=1, columnspan=5, sticky=tk.EW, padx=6)
        ttk.Button(frame, text="Open", command=self.choose_csv).grid(row=0, column=6, padx=(0, 6))
        ttk.Button(frame, text="Reload", command=self.reload_csv).grid(row=0, column=7)

        ttk.Label(frame, text="X").grid(row=1, column=0, sticky=tk.W, pady=(8, 0))
        self.x_combo = ttk.Combobox(frame, textvariable=self.x_column, state="readonly")
        self.x_combo.grid(row=1, column=1, sticky=tk.EW, padx=6, pady=(8, 0))

        ttk.Label(frame, text="Label").grid(row=1, column=2, sticky=tk.W, pady=(8, 0))
        self.label_combo = ttk.Combobox(frame, textvariable=self.label_column, state="readonly")
        self.label_combo.grid(row=1, column=3, sticky=tk.EW, padx=6, pady=(8, 0))

        ttk.Label(frame, text="Y unit").grid(row=1, column=4, sticky=tk.W, pady=(8, 0))
        self.unit_combo = ttk.Combobox(
            frame,
            textvariable=self.unit,
            state="readonly",
            values=("ns", "us", "ms", "s"),
            width=8,
        )
        self.unit_combo.grid(row=1, column=5, sticky=tk.W, padx=6, pady=(8, 0))

        ttk.Button(frame, text="Plot", command=self.plot).grid(row=1, column=6, columnspan=2, sticky=tk.EW, pady=(8, 0))

        ttk.Label(frame, text="Label filter").grid(row=2, column=0, sticky=tk.W, pady=(8, 0))
        self.label_filter_button = ttk.Menubutton(frame, text="All label values", state=tk.DISABLED)
        self.label_filter_button.grid(row=2, column=1, columnspan=5, sticky=tk.EW, padx=6, pady=(8, 0))
        self.label_filter_menu = tk.Menu(self.label_filter_button, tearoff=False)
        self.label_filter_button.configure(menu=self.label_filter_menu)

        ttk.Label(frame, textvariable=self.status).grid(row=3, column=0, columnspan=8, sticky=tk.W, pady=(8, 0))

        frame.columnconfigure(1, weight=1)
        frame.columnconfigure(3, weight=1)

        self.x_column.trace_add("write", lambda *_: self.plot())
        self.unit.trace_add("write", lambda *_: self.plot())
        self.label_column.trace_add("write", lambda *_: self._on_label_column_changed())

    def _build_plot(self) -> None:
        self.figure = Figure(figsize=(10, 6), dpi=100)
        self.axis = self.figure.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        toolbar = NavigationToolbar2Tk(self.canvas, self, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(side=tk.TOP, fill=tk.X)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def choose_csv(self) -> None:
        filename = filedialog.askopenfilename(
            title="Open benchmark CSV",
            filetypes=(("CSV files", "*.csv"), ("All files", "*.*")),
        )
        if filename:
            self.load_csv(Path(filename))

    def reload_csv(self) -> None:
        path = Path(self.csv_path.get())
        if not path:
            return
        self.load_csv(path)

    def load_csv(self, path: Path) -> None:
        try:
            self.data = BenchmarkData(path)
        except Exception as exc:
            messagebox.showerror("Could not load benchmark CSV", str(exc))
            return

        self.csv_path.set(str(path))
        columns = self.data.plottable_columns()
        labels = [NO_LABEL] + columns
        self.x_combo.configure(values=columns)
        self.label_combo.configure(values=labels)

        if not self.x_column.get() or self.x_column.get() not in columns:
            preferred = "Timestamp" if "Timestamp" in columns else columns[0] if columns else ""
            self.x_column.set(preferred)
        if self.label_column.get() not in labels:
            self.label_column.set(NO_LABEL)
        else:
            self._refresh_label_filter_options()

        self.status.set(f"Loaded {len(self.data.rows)} rows from {path}")
        self.plot()

    def _on_label_column_changed(self) -> None:
        self._refresh_label_filter_options()
        self.plot()

    def _refresh_label_filter_options(self) -> None:
        if self.data is None or self.label_filter_menu is None or self.label_filter_button is None:
            return

        had_existing_filter = bool(self.label_filter_vars)
        previous_selected = self.selected_label_values()
        self.label_filter_menu.delete(0, tk.END)
        self.label_filter_vars = {}

        label_column = self.label_column.get()
        if label_column == NO_LABEL:
            self.label_filter_button.configure(text="All label values", state=tk.DISABLED)
            return

        values = sorted({row.get(label_column, "") for row in self.data.rows}, key=natural_key)
        if not values:
            self.label_filter_button.configure(text="No label values", state=tk.DISABLED)
            return

        self._updating_label_filter = True
        for value in values:
            is_selected = value in previous_selected if had_existing_filter else True
            self.label_filter_vars[value] = tk.BooleanVar(value=is_selected)
        self._updating_label_filter = False

        self.label_filter_menu.add_command(label="Select all", command=self.select_all_label_values)
        self.label_filter_menu.add_command(label="Clear", command=self.clear_label_values)
        self.label_filter_menu.add_separator()
        for value in values:
            self.label_filter_menu.add_checkbutton(
                label=value if value else "<empty>",
                variable=self.label_filter_vars[value],
                command=self._on_label_filter_changed,
            )

        self.label_filter_button.configure(state=tk.NORMAL)
        self._update_label_filter_button_text()

    def _on_label_filter_changed(self) -> None:
        if self._updating_label_filter:
            return
        self._update_label_filter_button_text()
        self.plot()

    def selected_label_values(self) -> set[str]:
        return {value for value, variable in self.label_filter_vars.items() if variable.get()}

    def select_all_label_values(self) -> None:
        self._updating_label_filter = True
        for variable in self.label_filter_vars.values():
            variable.set(True)
        self._updating_label_filter = False
        self._on_label_filter_changed()

    def clear_label_values(self) -> None:
        self._updating_label_filter = True
        for variable in self.label_filter_vars.values():
            variable.set(False)
        self._updating_label_filter = False
        self._on_label_filter_changed()

    def _update_label_filter_button_text(self) -> None:
        if self.label_filter_button is None:
            return

        selected = self.selected_label_values()
        total = len(self.label_filter_vars)
        if total == 0:
            text = "All label values"
        elif len(selected) == total:
            text = "All label values"
        elif len(selected) == 0:
            text = "No label values"
        elif len(selected) == 1:
            text = next(iter(selected)) or "<empty>"
        else:
            text = f"{len(selected)} label values"
        self.label_filter_button.configure(text=text)

    def plot(self) -> None:
        if self.data is None:
            return

        x_column = self.x_column.get()
        label_column = self.label_column.get()
        if not x_column:
            return
        if label_column == NO_LABEL:
            label_column = None

        try:
            self.data.reload()
        except Exception as exc:
            self.status.set(f"Reload failed: {exc}")
            return

        factor = {"ns": 1.0, "us": 1_000.0, "ms": 1_000_000.0, "s": 1_000_000_000.0}[self.unit.get()]
        self._refresh_label_filter_options()
        grouped = self.data.grouped_average(x_column, label_column)
        selected_labels = self.selected_label_values() if label_column else set()
        if selected_labels:
            grouped = {label: values for label, values in grouped.items() if label in selected_labels}
        elif label_column and self.label_filter_vars:
            grouped = {}
        x_values = sorted({x for values in grouped.values() for x in values}, key=natural_key)

        self.axis.clear()
        if not x_values:
            self.axis.set_title("No plottable benchmark rows")
            self.canvas.draw_idle()
            return

        positions = list(range(len(x_values)))
        for label in sorted(grouped, key=natural_key):
            y_values = [grouped[label].get(x_value) for x_value in x_values]
            plot_x = [position for position, value in zip(positions, y_values) if value is not None]
            plot_y = [value / factor for value in y_values if value is not None]
            self.axis.plot(plot_x, plot_y, marker="o", label=label)

        self.axis.set_xticks(positions)
        self.axis.set_xticklabels(x_values, rotation=35, ha="right")
        self.axis.set_xlabel(x_column)
        self.axis.set_ylabel(f"Average duration [{self.unit.get()}]")
        self.axis.grid(True, axis="y", alpha=0.3)
        if label_column:
            self.axis.legend(title=label_column)
        self.figure.tight_layout()
        self.canvas.draw_idle()
        self.status.set(f"Loaded {len(self.data.rows)} rows from {self.data.path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Visualize HASE benchmark CSV files.")
    parser.add_argument("csv", nargs="?", type=Path, help="Path to hase_benchmark.csv")
    args = parser.parse_args()

    app = BenchmarkGui(args.csv)
    app.mainloop()


if __name__ == "__main__":
    main()
