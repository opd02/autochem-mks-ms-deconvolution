# AutoChem / MKS MS Deconvolution

This ready-to-run Python tool decomposes each time-resolved mass spectrum into a nonnegative mixture of expected species. It was preconfigured for:

- methyl propionate feed over Pt/C at 350 °C
- N₂ carrier
- ethane, ethylene, CO, CO₂, water, and methyl acrylate as candidate products

It accepts wide data (one column per m/z) or long data (time, m/z, signal), from CSV/TSV/TXT or Excel.

## What it reports

- Deconvolved signal coefficient versus time for every component
- Relative fitted-signal percentage versus time
- Absolute ppm only when you enter calibration sensitivities
- Reconstructed spectra, residuals, per-time fit R²/RMSE, and diagnostic plots
- Pairwise component-pattern similarity, which exposes components that cannot be separated reliably
- Optional unlabeled components extracted from residual signal

The model is nonnegative least squares (NNLS), not a black-box classifier:

`measured spectrum ≈ reference-pattern matrix × nonnegative component signals`

This is generally the appropriate starting model when expected species are known and unit-mass EI fragmentation is approximately linear.

## 1. Install

Install Python 3.10 or newer from <https://www.python.org/downloads/>. On Windows, check **Add Python to PATH** during installation.

Open Command Prompt or PowerShell in this folder:

```bash
python -m venv .venv
```

Windows:

```bat
.venv\Scripts\activate
python -m pip install -r requirements.txt
```

macOS/Linux:

```bash
source .venv/bin/activate
python -m pip install -r requirements.txt
```

## 2. Run

The easiest Windows method is to drag your CSV/XLSX export onto `run_windows.bat`.

From a terminal:

```bash
python ms_deconvolution.py "C:\path\to\run.csv"
```

Results go into `run_results` beside the input file.

For Excel:

```bash
python ms_deconvolution.py run.xlsx --sheet Data
```

## Accepted input layouts

### Wide format

Headings such as `18`, `m/z 18`, `Mass 18`, and `18 amu` are recognized:

```text
Elapsed Time (s),m/z 2,m/z 3,...,m/z 100
0,1.2,0.8,...,3.1
10,1.1,0.9,...,3.2
```

If time is in minutes/hours and the heading does not say so:

```bash
python ms_deconvolution.py run.csv --time-unit minutes
```

### Long format

```text
Time,Mass,Signal
0,18,1.2
0,28,20.1
10,18,1.3
```

If headings are unusual:

```bash
python ms_deconvolution.py run.csv --format long ^
  --time-column "Cycle Time" --mass-column "AMU" --signal-column "SEM Current"
```

## 3. Set the baseline correctly

By default, the program subtracts the median of the first 5% of points (at least 3, at most 100). Those points should be a carrier-only or pre-reaction baseline.

Specify an elapsed-time interval in seconds:

```bash
python ms_deconvolution.py run.csv --baseline-start 0 --baseline-end 300
```

Or use exactly 50 initial scans:

```bash
python ms_deconvolution.py run.csv --baseline-points 50
```

To analyze an already baseline-corrected file:

```bash
python ms_deconvolution.py run.csv --baseline none
```

## 4. Replace starter spectra with your instrument spectra

Edit `analytes.csv`. Each row contains:

| Column | Meaning |
|---|---|
| `species` | Component name |
| `mz` | Nominal mass channel |
| `relative_intensity` | Fragment intensity on any consistent relative scale |
| `sensitivity_signal_per_ppm` | Optional fitted coefficient per ppm |
| `enabled` | `true` or `false` |
| `notes` | Free text |

The included small-molecule EI patterns are reasonable starter patterns. The ester patterns are deliberately labeled approximate. Fragment ratios depend on electron energy, ion source, tuning, pressure, capillary, and detector. For publishable selectivity:

1. Record N₂-only baseline at the same flow/pressure.
2. Record methyl propionate feed over an inert bed or bypass at the same inlet conditions.
3. Run known mixtures or individual standards for CO, CO₂, ethane, ethylene, and methyl acrylate.
4. Replace relative intensities in `analytes.csv` with the observed ratios.
5. Fit sensitivity versus known ppm for each gas and enter the slope in `sensitivity_signal_per_ppm`.

Repeat the same sensitivity value on all rows for one species, or place it on one row.

## Critical identifiability issue for this experiment

N₂, CO, ethylene, ethane, and a CO₂ fragment all contribute at m/z 28. Because N₂ is the carrier, small changes in total flow or pressure can dwarf the product signals there.

- Do not quantify CO or C₂ products from m/z 28 alone.
- Retain m/z 26, 27, 29, and 30 for C₂ discrimination.
- Use m/z 44 for CO₂, m/z 18 for water, and 55/59/86/88 for the esters.
- A carrier-only baseline helps, but instrument-specific calibration is still necessary.
- `component_similarity.csv` and warnings in `analysis_report.txt` show when patterns are too similar.

Nitrogen is included but disabled in `analytes.csv`. Enable it only if you need to model carrier drift and have stable instrument-specific N₂ fragments. Doing so can make CO/ethylene estimates less identifiable.

## Useful options

Exclude background or unreliable channels:

```bash
python ms_deconvolution.py run.csv --exclude-masses 28,32,40
```

Use only selected masses:

```bash
python ms_deconvolution.py run.csv --include-masses 12,16-18,25-31,43-45,55,57-60,85-89
```

Weight by inverse noise measured during the baseline:

```bash
python ms_deconvolution.py run.csv --noise-weighting
```

Smooth the fitted time profiles with an 11-point Savitzky–Golay filter:

```bash
python ms_deconvolution.py run.csv --smooth-window 11
```

Find two unlabeled patterns remaining in the positive residuals:

```bash
python ms_deconvolution.py run.csv --unknown-components 2
```

Unknown components are exploratory NMF factors, not chemical identifications. Inspect their major m/z peaks and add chemically plausible candidates to `analytes.csv`.

## Test it before using real data

```bash
python ms_deconvolution.py --make-demo demo.csv
python ms_deconvolution.py demo.csv --output demo_results
```

Compare `demo_results/deconvolved_components.csv` with `demo_truth.csv`.

## Output files

| File | Contents |
|---|---|
| `deconvolved_components.csv` | Signals, relative fitted percentages, optional ppm, R², RMSE |
| `component_overview.png` | Main time-profile plots |
| `fit_quality.png` | Intense raw channels and spectrum R² |
| `residual_heatmap.png` | Missing/poorly fitted signal by time and m/z |
| `reconstructed_spectra_long.csv` | Raw, corrected, fitted, and residual values |
| `component_similarity.csv` | Pairwise spectral overlap |
| `baseline_and_noise.csv` | Subtracted baselines and measured noise |
| `analysis_report.txt/.json` | Parsing, fitting, and warning summary |

## Interpretation limits

`relative_pct_*` is the share of the **fitted MS signal**, not automatically mol% or carbon selectivity. Different molecules have different ionization cross-sections, fragmentation yields, transfer efficiencies, and detector responses. True gas concentration requires calibration. Carbon selectivity additionally requires molar flow and carbon-number accounting.

To avoid meaningless percentages when only baseline noise is present, relative percentages are set to zero until the total fitted signal exceeds three times the combined baseline-noise norm. The exact threshold is recorded in `analysis_report.txt`.

Reference spectra are a starting point for deconvolution, not proof that a component is present. Treat a product assignment as strong only when multiple characteristic ions rise together, the fitted coefficient exceeds baseline noise, fit residuals improve, and a standard or orthogonal method supports it.

## Reference-data notes

The default species and molecular masses were checked against the NIST Chemistry WebBook. NIST provides EI spectra but notes licensing restrictions on spectrum downloads for some compounds, so this package does not redistribute NIST spectral files. Open reference spectra can also be obtained from MassBank, whose records are openly versioned.

- NIST Chemistry WebBook: <https://webbook.nist.gov/>
- MassBank: <https://massbank.eu/MassBank/>
- MassBank open-data repository: <https://github.com/MassBank/MassBank-data>
