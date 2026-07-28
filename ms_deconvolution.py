#!/usr/bin/env python3
"""
AutoChem / residual-gas-analyzer mass-spectrometry deconvolution.

Fits each time-resolved mass spectrum as a nonnegative linear combination of
user-editable reference fragmentation patterns. This is intended for unit-mass
EI/RGA data (for example, an MKS mass spectrometer coupled to an AutoChem III).

Run:
    python ms_deconvolution.py path/to/export.csv
    python ms_deconvolution.py --make-demo demo.csv
    python ms_deconvolution.py demo.csv --output demo_results
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
import tempfile
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "autochem_ms_matplotlib"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.signal import savgol_filter


VERSION = "1.1.0"
SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_LIBRARY = SCRIPT_DIR / "analytes.csv"


@dataclass
class ParsedMS:
    time_seconds: np.ndarray
    masses: np.ndarray
    signals: np.ndarray  # time x mass
    time_label: str
    source_format: str
    original_time: np.ndarray | None = None
    scan_ids: np.ndarray | None = None
    source_delimiter: str | None = None


def _normalized_name(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


def _mass_from_column(value: object) -> float | None:
    """Recognize common wide-format mass-channel headings."""
    raw = str(value).strip()
    normalized = raw.lower().replace("−", "-")
    patterns = (
        r"^\s*(\d+(?:\.\d+)?)\s*$",
        r"^\s*(?:m\s*/?\s*z|mz|mass|amu)\s*[:=_\-\[\( ]*\s*(\d+(?:\.\d+)?)",
        r"^\s*(\d+(?:\.\d+)?)\s*(?:amu|u|m/z|mz)\s*$",
        r"^\s*(?:signal|current|pressure|intensity)\s*[:=_\-\[\( ]*\s*(\d+(?:\.\d+)?)\s*(?:amu|u|m/z|mz)?",
    )
    for pattern in patterns:
        match = re.search(pattern, normalized)
        if match:
            mass = float(match.group(1))
            if 0 < mass <= 1000:
                return mass
    return None


def _find_column(columns: Iterable[object], explicit: str | None, candidates: tuple[str, ...]) -> object | None:
    columns = list(columns)
    if explicit:
        for col in columns:
            if str(col) == explicit or _normalized_name(col) == _normalized_name(explicit):
                return col
        raise ValueError(f"Requested column {explicit!r} was not found. Columns: {columns}")
    normalized_candidates = tuple(_normalized_name(x) for x in candidates)
    for col in columns:
        normalized = _normalized_name(col)
        if normalized in normalized_candidates:
            return col
    for col in columns:
        normalized = _normalized_name(col)
        if any(candidate in normalized for candidate in normalized_candidates):
            return col
    return None


def read_table(path: Path, sheet: str | int | None = None) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        selected_sheet: str | int = 0 if sheet is None else sheet
        frame = pd.read_excel(path, sheet_name=selected_sheet)
        frame.attrs["source_delimiter"] = "Excel worksheet"
        return frame

    likely_separator = "\t" if suffix in {".txt", ".tsv"} else ","
    separators = [likely_separator, None, "\t", ",", ";"]
    attempts = []
    seen: set[str | None] = set()
    for separator in separators:
        if separator not in seen:
            attempts.append({"sep": separator, "engine": "python", "comment": "#"})
            seen.add(separator)
    errors = []
    for options in attempts:
        try:
            frame = pd.read_csv(path, **options)
            if frame.shape[1] >= 2:
                delimiter_names = {"\t": "tab", ",": "comma", ";": "semicolon", None: "auto-detected"}
                frame.attrs["source_delimiter"] = delimiter_names.get(options["sep"], repr(options["sep"]))
                return frame
        except Exception as exc:  # pragma: no cover - message collected for user
            errors.append(str(exc))
    raise ValueError("Could not parse the input as CSV/TSV. " + " | ".join(errors[:2]))


def _time_to_seconds(series: pd.Series, name: object, unit: str) -> tuple[np.ndarray, str]:
    name_text = str(name)
    if pd.api.types.is_datetime64_any_dtype(series):
        parsed = pd.to_datetime(series)
        seconds = (parsed - parsed.iloc[0]).dt.total_seconds().to_numpy(float)
        return seconds, name_text

    numeric = pd.to_numeric(series, errors="coerce")
    numeric_fraction = float(numeric.notna().mean())
    if numeric_fraction >= 0.80:
        values = numeric.to_numpy(float)
        if not np.isfinite(values).all():
            values = pd.Series(values).interpolate(limit_direction="both").to_numpy()
        inferred = unit
        normalized = _normalized_name(name)
        if unit == "auto":
            if "hour" in normalized or normalized.endswith("hr") or normalized.endswith("h"):
                inferred = "hours"
            elif "min" in normalized:
                inferred = "minutes"
            else:
                inferred = "seconds"
        factor = {"seconds": 1.0, "minutes": 60.0, "hours": 3600.0}[inferred]
        return (values - values[0]) * factor, f"{name_text} ({inferred})"

    timedelta = pd.to_timedelta(series.astype(str), errors="coerce")
    if timedelta.notna().mean() >= 0.80:
        seconds = timedelta.dt.total_seconds().to_numpy(float)
        return seconds - seconds[0], name_text

    datetimes = pd.to_datetime(series.astype(str), errors="coerce", format="mixed")
    if datetimes.notna().mean() >= 0.80:
        seconds = (datetimes - datetimes.iloc[0]).dt.total_seconds().to_numpy(float)
        return seconds, name_text

    raise ValueError(
        f"Time column {name_text!r} could not be interpreted. "
        "Use a numeric elapsed-time column or a datetime/timedelta column."
    )


def parse_input(
    frame: pd.DataFrame,
    data_format: str,
    time_column: str | None,
    mass_column: str | None,
    signal_column: str | None,
    time_unit: str,
) -> ParsedMS:
    source_delimiter = frame.attrs.get("source_delimiter")
    frame = frame.dropna(axis=0, how="all").dropna(axis=1, how="all")
    if frame.empty:
        raise ValueError("The input table is empty.")

    time_col = _find_column(
        frame.columns,
        time_column,
        ("time", "elapsedtime", "elapsed", "seconds", "second", "minutes", "minute", "timestamp", "datetime"),
    )
    mass_col = _find_column(frame.columns, mass_column, ("mz", "m/z", "mass", "amu", "atomicmass"))
    signal_col = _find_column(
        frame.columns,
        signal_column,
        ("signal", "intensity", "current", "ioncurrent", "partialpressure", "pressure", "value"),
    )
    scan_col = _find_column(frame.columns, None, ("scan", "scannumber", "scanindex"))

    if data_format == "auto":
        data_format = "long" if mass_col is not None and signal_col is not None else "wide"

    if data_format == "long":
        if mass_col is None or signal_col is None:
            raise ValueError(
                "Long format needs mass and signal columns. Use --mass-column and --signal-column "
                f"to name them explicitly. Columns: {list(frame.columns)}"
            )
        if time_col is None:
            raise ValueError("Long format needs a time column; use --time-column.")
        working = frame[[time_col, mass_col, signal_col]].copy()
        working[mass_col] = pd.to_numeric(working[mass_col], errors="coerce")
        working[signal_col] = pd.to_numeric(working[signal_col], errors="coerce")
        working = working.dropna()
        if working.empty:
            raise ValueError("No numeric mass/signal rows remained after parsing.")
        time_values, time_label = _time_to_seconds(working[time_col], time_col, time_unit)
        working["_time_seconds"] = time_values
        pivoted = (
            working.groupby(["_time_seconds", mass_col], as_index=False)[signal_col]
            .mean()
            .pivot(index="_time_seconds", columns=mass_col, values=signal_col)
            .sort_index()
            .sort_index(axis=1)
        )
        return ParsedMS(
            time_seconds=pivoted.index.to_numpy(float),
            masses=pivoted.columns.to_numpy(float),
            signals=pivoted.to_numpy(float),
            time_label=time_label,
            source_format="long",
            source_delimiter=source_delimiter,
        )

    mass_columns: list[tuple[object, float]] = []
    for col in frame.columns:
        if col == time_col:
            continue
        mass = _mass_from_column(col)
        if mass is not None:
            mass_columns.append((col, mass))
    if len(mass_columns) < 2:
        raise ValueError(
            "Could not find at least two mass channels in wide format. Expected headings such as "
            "'18', 'm/z 18', 'Mass 18', or '18 amu'. Use long format if the file has one mass per row. "
            f"Columns: {list(frame.columns)}"
        )

    if time_col is None:
        time_values = np.arange(len(frame), dtype=float)
        time_label = "row index (assumed 1 second per row)"
        original_time = None
    else:
        time_values, time_label = _time_to_seconds(frame[time_col], time_col, time_unit)
        original_time = frame[time_col].astype(str).to_numpy()
    scan_ids = (
        pd.to_numeric(frame[scan_col], errors="coerce").to_numpy()
        if scan_col is not None
        else None
    )

    numeric = pd.DataFrame(index=frame.index)
    for col, mass in mass_columns:
        numeric[str(col)] = pd.to_numeric(frame[col], errors="coerce")
    valid_rows = np.isfinite(time_values) & (numeric.notna().mean(axis=1).to_numpy() >= 0.50)
    numeric = numeric.loc[valid_rows].interpolate(limit_direction="both")
    time_values = time_values[valid_rows]
    if original_time is not None:
        original_time = original_time[valid_rows]
    if scan_ids is not None:
        scan_ids = scan_ids[valid_rows]
    if numeric.empty:
        raise ValueError("No usable numeric signal rows remained after parsing.")

    masses = np.array([mass for _, mass in mass_columns], dtype=float)
    values = numeric.to_numpy(float)
    # Average duplicate channels representing the same nominal mass.
    rounded = np.rint(masses).astype(int)
    unique_masses = np.unique(rounded)
    combined = np.column_stack([np.nanmean(values[:, rounded == mass], axis=1) for mass in unique_masses])
    order = np.argsort(time_values)
    is_mks_layout = (
        scan_col is not None
        and time_col is not None
        and all(str(col).strip().lower().startswith("mass ") for col, _ in mass_columns)
    )
    return ParsedMS(
        time_seconds=time_values[order],
        masses=unique_masses.astype(float),
        signals=combined[order],
        time_label=time_label,
        source_format="wide MKS scan export" if is_mks_layout else "wide",
        original_time=original_time[order] if original_time is not None else None,
        scan_ids=scan_ids[order] if scan_ids is not None else None,
        source_delimiter=source_delimiter,
    )


def load_library(path: Path) -> tuple[pd.DataFrame, list[str], np.ndarray, dict[str, float | None]]:
    library = pd.read_csv(path, comment="#")
    required = {"species", "mz", "relative_intensity"}
    missing = required - set(library.columns)
    if missing:
        raise ValueError(f"Fragment library is missing columns: {sorted(missing)}")
    if "enabled" in library.columns:
        enabled = library["enabled"].astype(str).str.strip().str.lower().isin({"1", "true", "yes", "y"})
        library = library.loc[enabled].copy()
    library["species"] = library["species"].astype(str).str.strip()
    library["mz"] = pd.to_numeric(library["mz"], errors="raise")
    library["relative_intensity"] = pd.to_numeric(library["relative_intensity"], errors="raise")
    library = library[library["relative_intensity"] > 0].copy()
    if library.empty:
        raise ValueError("No enabled, positive-intensity library rows were found.")
    species = list(dict.fromkeys(library["species"]))
    library_masses = np.sort(np.unique(np.rint(library["mz"]).astype(int)))
    sensitivity: dict[str, float | None] = {}
    for name in species:
        rows = library[library["species"] == name]
        if "sensitivity_signal_per_ppm" in library.columns:
            values = pd.to_numeric(rows["sensitivity_signal_per_ppm"], errors="coerce").dropna().unique()
            sensitivity[name] = float(values[0]) if len(values) else None
        else:
            sensitivity[name] = None
    return library, species, library_masses, sensitivity


def build_pattern_matrix(
    library: pd.DataFrame, species: list[str], measured_masses: np.ndarray
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    measured_nominal = np.rint(measured_masses).astype(int)
    matrix = np.zeros((len(measured_masses), len(species)), dtype=float)
    warnings: list[str] = []
    for j, name in enumerate(species):
        rows = library[library["species"] == name]
        for _, row in rows.iterrows():
            nominal = int(round(float(row["mz"])))
            indices = np.where(measured_nominal == nominal)[0]
            if len(indices):
                matrix[indices[0], j] += float(row["relative_intensity"])
        column_sum = matrix[:, j].sum()
        if column_sum <= 0:
            warnings.append(f"{name}: none of its library masses are present in the measured data; component removed.")
        else:
            matrix[:, j] /= column_sum
    keep = matrix.sum(axis=0) > 0
    return matrix[:, keep], keep, warnings


def parse_mass_list(text: str | None) -> set[int]:
    if not text:
        return set()
    result: set[int] = set()
    for token in re.split(r"[,;\s]+", text.strip()):
        if not token:
            continue
        if "-" in token:
            start_text, end_text = token.split("-", 1)
            start, end = int(float(start_text)), int(float(end_text))
            result.update(range(min(start, end), max(start, end) + 1))
        else:
            result.add(int(round(float(token))))
    return result


def baseline_correct(
    signals: np.ndarray,
    time_seconds: np.ndarray,
    mode: str,
    baseline_points: int | None,
    baseline_start: float | None,
    baseline_end: float | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    if mode == "none":
        zeros = np.zeros(signals.shape[1], dtype=float)
        noise = np.nanstd(signals[: max(3, min(20, len(signals)))], axis=0)
        return np.clip(signals, 0, None), zeros, noise, "No baseline subtraction"

    if baseline_start is not None or baseline_end is not None:
        start = -np.inf if baseline_start is None else baseline_start
        end = np.inf if baseline_end is None else baseline_end
        mask = (time_seconds >= start) & (time_seconds <= end)
        description = f"Median over {start:g}–{end:g} s"
    else:
        n_points = baseline_points
        if n_points is None:
            n_points = max(3, min(100, int(math.ceil(len(signals) * 0.05))))
        n_points = min(max(1, n_points), len(signals))
        mask = np.zeros(len(signals), dtype=bool)
        mask[:n_points] = True
        description = f"Median of first {n_points} points"
    if mask.sum() < 1:
        raise ValueError("The selected baseline interval contains no data points.")
    baseline = np.nanmedian(signals[mask], axis=0)
    noise = 1.4826 * np.nanmedian(np.abs(signals[mask] - baseline), axis=0)
    corrected = signals - baseline
    # Negative residuals cannot be represented by a nonnegative mixture.
    corrected = np.clip(corrected, 0, None)
    return corrected, baseline, noise, description


def choose_weights(noise: np.ndarray, use_noise_weights: bool) -> np.ndarray:
    if not use_noise_weights:
        return np.ones_like(noise)
    positive = noise[np.isfinite(noise) & (noise > 0)]
    floor = float(np.median(positive) * 0.25) if len(positive) else 1.0
    floor = max(floor, np.finfo(float).eps)
    return 1.0 / np.maximum(np.nan_to_num(noise, nan=floor), floor)


def fit_nnls(signals: np.ndarray, patterns: np.ndarray, weights: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    weighted_patterns = patterns * weights[:, None]
    coefficients = np.zeros((signals.shape[0], patterns.shape[1]), dtype=float)
    fitted = np.zeros_like(signals)
    for i, spectrum in enumerate(signals):
        coefficients[i], _ = nnls(weighted_patterns, spectrum * weights)
        fitted[i] = patterns @ coefficients[i]
    return coefficients, fitted


def smooth_coefficients(coefficients: np.ndarray, window: int) -> np.ndarray:
    if window <= 1 or len(coefficients) < 3:
        return coefficients
    if window % 2 == 0:
        window += 1
    window = min(window, len(coefficients) if len(coefficients) % 2 else len(coefficients) - 1)
    if window < 3:
        return coefficients
    polyorder = min(2, window - 1)
    smoothed = savgol_filter(coefficients, window_length=window, polyorder=polyorder, axis=0, mode="interp")
    return np.clip(smoothed, 0, None)


def pairwise_similarity(patterns: np.ndarray, species: list[str]) -> list[dict[str, float | str]]:
    norms = np.linalg.norm(patterns, axis=0)
    records: list[dict[str, float | str]] = []
    for i in range(len(species)):
        for j in range(i + 1, len(species)):
            denominator = norms[i] * norms[j]
            similarity = float(patterns[:, i] @ patterns[:, j] / denominator) if denominator else 0.0
            records.append({"species_1": species[i], "species_2": species[j], "cosine_similarity": similarity})
    return sorted(records, key=lambda item: float(item["cosine_similarity"]), reverse=True)


def estimate_unknown_profiles(
    positive_residual: np.ndarray, n_components: int, iterations: int = 500
) -> tuple[np.ndarray, np.ndarray]:
    """Small dependency-free nonnegative matrix factorization for residual discovery."""
    if n_components <= 0:
        return np.empty((len(positive_residual), 0)), np.empty((0, positive_residual.shape[1]))
    rng = np.random.default_rng(42)
    scale = max(float(np.nanmean(positive_residual)), np.finfo(float).eps)
    w = rng.random((positive_residual.shape[0], n_components)) * math.sqrt(scale)
    h = rng.random((n_components, positive_residual.shape[1])) * math.sqrt(scale)
    eps = np.finfo(float).eps
    for _ in range(iterations):
        h *= (w.T @ positive_residual) / (w.T @ w @ h + eps)
        w *= (positive_residual @ h.T) / (w @ h @ h.T + eps)
    # Normalize spectra to sum one, keeping scale in temporal profiles.
    row_sums = h.sum(axis=1)
    row_sums[row_sums <= 0] = 1.0
    h /= row_sums[:, None]
    w *= row_sums[None, :]
    return w, h


def safe_filename(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_") or "component"


def save_overview_plot(
    output: Path,
    time_seconds: np.ndarray,
    coefficients: np.ndarray,
    relative_pct: np.ndarray,
    species: list[str],
    calibrated: dict[str, np.ndarray],
) -> None:
    x = time_seconds / 60.0
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True, constrained_layout=True)
    for j, name in enumerate(species):
        axes[0].plot(x, coefficients[:, j], label=name, linewidth=1.6)
    axes[0].set_ylabel("Deconvolved signal coefficient")
    axes[0].set_title("MS component signals")
    axes[0].legend(ncol=min(4, max(1, len(species))), fontsize=8)
    axes[0].grid(alpha=0.2)

    axes[1].stackplot(x, relative_pct.T, labels=species, alpha=0.85)
    axes[1].set_ylabel("Relative fitted signal (%)")
    axes[1].set_xlabel("Elapsed time (min)")
    axes[1].set_ylim(0, 100)
    axes[1].grid(alpha=0.2)
    fig.savefig(output / "component_overview.png", dpi=180)
    plt.close(fig)

    if calibrated:
        fig, ax = plt.subplots(figsize=(12, 5), constrained_layout=True)
        for name, values in calibrated.items():
            ax.plot(x, values, label=name, linewidth=1.6)
        ax.set_xlabel("Elapsed time (min)")
        ax.set_ylabel("Calibrated concentration (ppm)")
        ax.set_title("Calibrated concentrations")
        ax.legend()
        ax.grid(alpha=0.2)
        fig.savefig(output / "calibrated_concentrations.png", dpi=180)
        plt.close(fig)


def save_fit_plots(
    output: Path,
    time_seconds: np.ndarray,
    masses: np.ndarray,
    corrected: np.ndarray,
    fitted: np.ndarray,
    r2: np.ndarray,
) -> None:
    x = time_seconds / 60.0
    residual = corrected - fitted
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True, constrained_layout=True)
    selected_indices = np.argsort(np.nanmax(corrected, axis=0))[-min(10, corrected.shape[1]) :]
    for index in selected_indices:
        axes[0].plot(x, corrected[:, index], label=f"m/z {masses[index]:g}", linewidth=1.0)
    axes[0].set_ylabel("Baseline-corrected signal")
    axes[0].set_title("Most intense measured channels")
    axes[0].legend(ncol=5, fontsize=8)
    axes[0].grid(alpha=0.2)
    axes[1].plot(x, r2, color="black", linewidth=1.2)
    axes[1].axhline(0.90, color="tab:red", linestyle="--", linewidth=1)
    axes[1].set_ylabel("Spectrum fit R²")
    axes[1].set_xlabel("Elapsed time (min)")
    axes[1].set_ylim(min(-0.1, float(np.nanmin(r2))), 1.03)
    axes[1].grid(alpha=0.2)
    fig.savefig(output / "fit_quality.png", dpi=180)
    plt.close(fig)

    if residual.size:
        limit = np.nanpercentile(np.abs(residual), 99)
        limit = max(float(limit), np.finfo(float).eps)
        fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)
        image = ax.imshow(
            residual.T,
            aspect="auto",
            origin="lower",
            extent=[x[0], x[-1], masses[0], masses[-1]],
            cmap="coolwarm",
            vmin=-limit,
            vmax=limit,
            interpolation="nearest",
        )
        ax.set_xlabel("Elapsed time (min)")
        ax.set_ylabel("m/z")
        ax.set_title("Fit residuals (measured − reconstructed)")
        fig.colorbar(image, ax=ax, label="Residual signal")
        fig.savefig(output / "residual_heatmap.png", dpi=180)
        plt.close(fig)


def make_demo(path: Path, library_path: Path, points: int = 300) -> None:
    library, species, library_masses, _ = load_library(library_path)
    masses = np.arange(2, 101, dtype=float)
    patterns, keep, _ = build_pattern_matrix(library, species, masses)
    species = [name for name, include in zip(species, keep) if include]
    t = np.arange(points, dtype=float) * 10.0
    c = np.zeros((points, len(species)), dtype=float)

    def add(name: str, values: np.ndarray) -> None:
        if name in species:
            c[:, species.index(name)] = values

    reaction = 1.0 / (1.0 + np.exp(-(t - 600.0) / 45.0))
    add("2-Butanone", reaction * (350.0 + 20.0 * np.sin(t / 280.0)))
    add("Hydrogen", 35.0 * reaction)
    add("Methane", 30.0 * reaction)
    add("Ethylene", 110.0 * reaction * np.exp(-t / 5000.0))
    add("Ethane", 60.0 * reaction)
    add("Propylene", 45.0 * reaction)
    add("Propane", 70.0 * reaction)
    add("Carbon monoxide", 45.0 * reaction)
    add("Carbon dioxide", 12.0 * reaction)
    add("Water", 10.0 * reaction)
    add("Methyl vinyl ketone", 55.0 * np.exp(-0.5 * ((t - 1600.0) / 450.0) ** 2))
    rng = np.random.default_rng(7)
    baseline = 2.0 + 0.02 * masses
    signals = c @ patterns.T + baseline[None, :] + rng.normal(0, 0.35, (points, len(masses)))
    signals = np.clip(signals, 0, None)
    frame = pd.DataFrame({"Elapsed Time (s)": t})
    for j, mass in enumerate(masses.astype(int)):
        frame[f"m/z {mass}"] = signals[:, j]
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)
    truth = pd.DataFrame({"time_seconds": t})
    for j, name in enumerate(species):
        truth[f"true_{safe_filename(name)}"] = c[:, j]
    truth.to_csv(path.with_name(path.stem + "_truth.csv"), index=False)


def analyze(args: argparse.Namespace) -> Path:
    input_path = Path(args.input).expanduser().resolve()
    library_path = Path(args.library).expanduser().resolve()
    output = Path(args.output).expanduser().resolve() if args.output else input_path.with_name(input_path.stem + "_results")
    output.mkdir(parents=True, exist_ok=True)

    frame = read_table(input_path, sheet=args.sheet)
    parsed = parse_input(
        frame,
        data_format=args.format,
        time_column=args.time_column,
        mass_column=args.mass_column,
        signal_column=args.signal_column,
        time_unit=args.time_unit,
    )
    library, all_species, _, sensitivity = load_library(library_path)

    include_masses = parse_mass_list(args.include_masses)
    exclude_masses = parse_mass_list(args.exclude_masses)
    nominal = np.rint(parsed.masses).astype(int)
    channel_mask = np.ones(len(parsed.masses), dtype=bool)
    if include_masses:
        channel_mask &= np.array([mass in include_masses for mass in nominal])
    if exclude_masses:
        channel_mask &= np.array([mass not in exclude_masses for mass in nominal])
    if channel_mask.sum() < 2:
        raise ValueError("Fewer than two mass channels remain after include/exclude filtering.")

    masses = parsed.masses[channel_mask]
    raw = parsed.signals[:, channel_mask]
    finite_columns = np.isfinite(raw).mean(axis=0) >= 0.80
    masses = masses[finite_columns]
    raw = pd.DataFrame(raw[:, finite_columns]).interpolate(limit_direction="both").to_numpy(float)

    corrected, baseline, noise, baseline_description = baseline_correct(
        raw,
        parsed.time_seconds,
        mode=args.baseline,
        baseline_points=args.baseline_points,
        baseline_start=args.baseline_start,
        baseline_end=args.baseline_end,
    )
    patterns, keep, warnings = build_pattern_matrix(library, all_species, masses)
    species = [name for name, include in zip(all_species, keep) if include]
    if patterns.shape[1] == 0:
        raise ValueError("No analyte has a reference ion present in the measured channels.")

    weights = choose_weights(noise, args.noise_weighting)
    coefficients, fitted = fit_nnls(corrected, patterns, weights)
    coefficients = smooth_coefficients(coefficients, args.smooth_window)
    fitted = coefficients @ patterns.T

    total = coefficients.sum(axis=1)
    detection_threshold = 3.0 * float(np.linalg.norm(noise)) if args.baseline != "none" else 0.0
    detected = total > detection_threshold
    relative_pct = np.divide(
        coefficients,
        total[:, None],
        out=np.zeros_like(coefficients),
        where=detected[:, None],
    ) * 100.0

    residual = corrected - fitted
    sse = np.sum(residual**2, axis=1)
    centered = corrected - np.mean(corrected, axis=1, keepdims=True)
    sst = np.sum(centered**2, axis=1)
    r2 = np.divide(sst - sse, sst, out=np.full_like(sst, np.nan), where=sst > 0)
    rmse = np.sqrt(np.mean(residual**2, axis=1))

    calibrated: dict[str, np.ndarray] = {}
    for j, name in enumerate(species):
        factor = sensitivity.get(name)
        if factor is not None and np.isfinite(factor) and factor > 0:
            calibrated[name] = coefficients[:, j] / factor

    results = pd.DataFrame()
    if parsed.original_time is not None:
        results["timestamp"] = parsed.original_time
    if parsed.scan_ids is not None:
        results["scan"] = parsed.scan_ids
    results["time_seconds"] = parsed.time_seconds
    results["time_minutes"] = parsed.time_seconds / 60.0
    results["total_fitted_signal"] = total
    results["signal_detected"] = detected
    for j, name in enumerate(species):
        slug = safe_filename(name)
        results[f"signal_{slug}"] = coefficients[:, j]
        results[f"relative_pct_{slug}"] = relative_pct[:, j]
        if name in calibrated:
            results[f"ppm_{slug}"] = calibrated[name]
    results["fit_r2"] = r2
    results["fit_rmse"] = rmse
    results.to_csv(output / "deconvolved_components.csv", index=False)

    reconstructed_columns: dict[str, np.ndarray] = {}
    if parsed.original_time is not None:
        reconstructed_columns["timestamp"] = np.repeat(parsed.original_time, len(masses))
    if parsed.scan_ids is not None:
        reconstructed_columns["scan"] = np.repeat(parsed.scan_ids, len(masses))
    reconstructed_columns.update(
        {
            "time_seconds": np.repeat(parsed.time_seconds, len(masses)),
            "mz": np.tile(masses, len(parsed.time_seconds)),
            "raw_signal": raw.ravel(),
            "baseline_corrected_signal": corrected.ravel(),
            "fitted_signal": fitted.ravel(),
            "residual": residual.ravel(),
        }
    )
    reconstructed = pd.DataFrame(reconstructed_columns)
    reconstructed.to_csv(output / "reconstructed_spectra_long.csv", index=False)

    baseline_table = pd.DataFrame({"mz": masses, "baseline": baseline, "baseline_noise_MAD": noise})
    baseline_table.to_csv(output / "baseline_and_noise.csv", index=False)

    similarities = pairwise_similarity(patterns, species)
    pd.DataFrame(similarities).to_csv(output / "component_similarity.csv", index=False)

    unknown_profiles, unknown_spectra = estimate_unknown_profiles(
        np.clip(residual, 0, None), args.unknown_components
    )
    if args.unknown_components > 0:
        unknown_time = pd.DataFrame(
            {"time_seconds": parsed.time_seconds, "time_minutes": parsed.time_seconds / 60.0}
        )
        unknown_library_rows = []
        for component in range(args.unknown_components):
            label = f"unknown_{component + 1}"
            unknown_time[label] = unknown_profiles[:, component]
            ranked = np.argsort(unknown_spectra[component])[::-1]
            for index in ranked:
                if unknown_spectra[component, index] > 0:
                    unknown_library_rows.append(
                        {"component": label, "mz": masses[index], "relative_intensity": unknown_spectra[component, index]}
                    )
        unknown_time.to_csv(output / "unknown_component_profiles.csv", index=False)
        pd.DataFrame(unknown_library_rows).to_csv(output / "unknown_component_spectra.csv", index=False)

    save_overview_plot(output, parsed.time_seconds, coefficients, relative_pct, species, calibrated)
    save_fit_plots(output, parsed.time_seconds, masses, corrected, fitted, r2)

    high_similarity = [item for item in similarities if float(item["cosine_similarity"]) >= 0.90]
    for item in high_similarity:
        warnings.append(
            f"Strong overlap: {item['species_1']} vs {item['species_2']} "
            f"(cosine similarity {float(item['cosine_similarity']):.3f})."
        )
    if 28 in np.rint(masses).astype(int):
        warnings.append(
            "m/z 28 contains contributions from N2, CO, ethylene, ethane, and CO2 fragments. "
            "Do not interpret that channel alone."
        )
    if 28 not in np.rint(masses).astype(int):
        warnings.append("m/z 28 is absent/excluded; CO and C2 hydrocarbon estimates may be weak.")
    if 43 in np.rint(masses).astype(int):
        warnings.append(
            "m/z 43 is shared by 2-butanone and several hydrocarbon/oxygenate fragments. "
            "Use the 2-butanone molecular ion at m/z 72 and companion ions for feed quantification."
        )
    nominal_masses = np.rint(masses).astype(int)
    if 40 in nominal_masses:
        mass40_index = int(np.where(nominal_masses == 40)[0][0])
        other_baselines = np.abs(np.delete(baseline, mass40_index))
        typical_baseline = float(np.nanmedian(other_baselines[np.isfinite(other_baselines)]))
        if typical_baseline > 0 and abs(float(baseline[mass40_index])) > 100.0 * typical_baseline:
            warnings.append(
                "m/z 40 is exceptionally large relative to other channels, consistent with an Ar carrier/background peak. "
                "If Ar is not analytically relevant, consider --exclude-masses 35-40 or --noise-weighting; "
                "inspect m/z 35-39 before excluding the full range."
            )
    median_r2 = float(np.nanmedian(r2)) if np.isfinite(r2).any() else float("nan")
    if np.isfinite(median_r2) and median_r2 < 0.90:
        warnings.append(
            f"Median fit R² is {median_r2:.3f}; the library may be missing products, "
            "the baseline may be inappropriate, or instrument-specific fragmentation differs."
        )

    report = {
        "program_version": VERSION,
        "input_file": str(input_path),
        "input_format": parsed.source_format,
        "input_delimiter": parsed.source_delimiter,
        "time_interpretation": parsed.time_label,
        "scan_ids_preserved": parsed.scan_ids is not None,
        "n_time_points": int(len(parsed.time_seconds)),
        "masses_used": [float(value) for value in masses],
        "components_fitted": species,
        "baseline": baseline_description,
        "noise_weighting": bool(args.noise_weighting),
        "savgol_smoothing_window_points": int(args.smooth_window),
        "relative_fraction_reporting_threshold": detection_threshold,
        "median_fit_r2": median_r2,
        "median_fit_rmse": float(np.nanmedian(rmse)),
        "calibrated_components": list(calibrated),
        "warnings": warnings,
    }
    (output / "analysis_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    readable_warnings = "\n".join(f"- {warning}" for warning in warnings) or "- None"
    report_text = f"""\
AutoChem MS deconvolution report
================================
Input: {input_path.name}
Detected format: {parsed.source_format}
Detected delimiter: {parsed.source_delimiter}
Time interpretation: {parsed.time_label}
Scan IDs preserved: {parsed.scan_ids is not None}
Time points: {len(parsed.time_seconds)}
Mass channels used: {", ".join(f"{x:g}" for x in masses)}
Components: {", ".join(species)}
Baseline: {baseline_description}
Median fit R²: {median_r2:.4f}
Median fit RMSE: {float(np.nanmedian(rmse)):.6g}

Interpretation
--------------
"signal_*" columns are deconvolved, pattern-scaled MS signals.
"relative_pct_*" columns are each component's fraction of the fitted signal,
not automatically a gas-phase mole percentage. Relative percentages are set
to zero below a total fitted-signal threshold of {detection_threshold:.6g}.

Absolute ppm values are only exported when a positive
sensitivity_signal_per_ppm is supplied for that species in analytes.csv.
Calibrate using your own instrument, ion-source settings, capillary, pressure,
and flow conditions.

Warnings / identifiability checks
---------------------------------
{readable_warnings}
"""
    (output / "analysis_report.txt").write_text(report_text, encoding="utf-8")
    return output


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Deconvolute time-resolved AutoChem/MKS mass-spectrometer data with nonnegative least squares.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """
            Examples:
              python ms_deconvolution.py run1.csv
              python ms_deconvolution.py run1.xlsx --sheet Data
              python ms_deconvolution.py run1.csv --baseline-start 0 --baseline-end 300
              python ms_deconvolution.py run1.csv --exclude-masses 28,32,40
              python ms_deconvolution.py --make-demo demo.csv
            """
        ),
    )
    parser.add_argument("input", nargs="?", help="CSV, TSV, TXT, XLSX, or XLS export")
    parser.add_argument("--output", help="Output folder (default: INPUT_results)")
    parser.add_argument("--library", default=str(DEFAULT_LIBRARY), help="Analyte fragment-pattern CSV")
    parser.add_argument("--format", choices=("auto", "wide", "long"), default="auto")
    parser.add_argument("--sheet", help="Excel sheet name")
    parser.add_argument("--time-column", help="Exact time-column heading")
    parser.add_argument("--mass-column", help="Exact mass-column heading for long format")
    parser.add_argument("--signal-column", help="Exact signal-column heading for long format")
    parser.add_argument("--time-unit", choices=("auto", "seconds", "minutes", "hours"), default="auto")
    parser.add_argument("--baseline", choices=("median", "none"), default="median")
    parser.add_argument("--baseline-points", type=int, help="Number of initial points used as baseline")
    parser.add_argument("--baseline-start", type=float, help="Baseline interval start, in elapsed seconds")
    parser.add_argument("--baseline-end", type=float, help="Baseline interval end, in elapsed seconds")
    parser.add_argument("--include-masses", help="Only fit these nominal masses, e.g. '12,16-18,26-30,44'")
    parser.add_argument("--exclude-masses", help="Exclude nominal masses, e.g. '32,40'")
    parser.add_argument(
        "--noise-weighting",
        action="store_true",
        help="Weight channels by inverse baseline noise (use after checking the baseline interval)",
    )
    parser.add_argument(
        "--smooth-window",
        type=int,
        default=1,
        help="Optional odd Savitzky-Golay window in points; 1 disables smoothing",
    )
    parser.add_argument(
        "--unknown-components",
        type=int,
        default=0,
        help="Extract this many unlabeled nonnegative components from positive residuals",
    )
    parser.add_argument("--make-demo", metavar="PATH", help="Write synthetic wide-format demo CSV and truth CSV")
    parser.add_argument("--version", action="version", version=VERSION)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    try:
        if args.make_demo:
            demo_path = Path(args.make_demo).expanduser().resolve()
            make_demo(demo_path, Path(args.library).expanduser().resolve())
            print(f"Demo written to: {demo_path}")
            print(f"Truth written to: {demo_path.with_name(demo_path.stem + '_truth.csv')}")
            return 0
        if not args.input:
            parser.error("Provide an input file, or use --make-demo PATH.")
        if args.smooth_window < 1:
            parser.error("--smooth-window must be at least 1.")
        if args.unknown_components < 0:
            parser.error("--unknown-components cannot be negative.")
        output = analyze(args)
        print(f"Analysis complete: {output}")
        print(f"Open: {output / 'component_overview.png'}")
        print(f"Read: {output / 'analysis_report.txt'}")
        return 0
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
