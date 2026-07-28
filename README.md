# AutoChem / MKS MS Deconvolution

This Python tool decomposes each time-resolved, unit-mass EI spectrum into a nonnegative mixture of expected species. Version 1.1 is configured for:

- 2-butanone vapor over Pt/C at 400 °C
- N₂ carrier, with optional Ar-background handling
- scans from m/z 2 through at least m/z 75
- quoted MKS tab-delimited exports as well as ordinary CSV and Excel files

The model is nonnegative least squares (NNLS):

`measured spectrum ≈ reference-pattern matrix × nonnegative component signals`

NNLS is transparent and appropriate when the expected species are known and signal/fragmentation is approximately linear. It cannot turn an ambiguous unit-mass spectrum into a unique chemical identification.

## What it reports

- Deconvolved component signals and relative fitted-signal percentages over time
- Original MKS timestamps and scan numbers
- Absolute ppm when instrument calibration sensitivities are provided
- Reconstructed spectra, residuals, R²/RMSE, and diagnostic plots
- Pairwise component-pattern similarity to expose non-identifiable combinations
- Optional unlabeled components extracted from positive residual signal

## Install

Install Python 3.10 or newer from <https://www.python.org/downloads/>. On Windows, check **Add Python to PATH**.

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

## Run an MKS export

The supplied MKS format is accepted directly, even when the filename ends in `.txt`:

```text
"Time"  "Scan"  "Mass 2"  "Mass 3" ... "Mass 75"
"12/4/2025 5:41:48 PM"  1  1.8118e+02 ... 9.4256e+00
```

The actual delimiter is a tab. Quoted headings, scientific notation, negative detector values, the trailing blank column, timestamps, and scan numbers are handled automatically:

```bash
python ms_deconvolution.py "225124ptcustab_000001 - Copy.txt"
```

Results are written beside the input in a folder ending in `_results`. On Windows, you can instead drag the export onto `run_windows.bat`.

Other supported formats:

- CSV/TSV/TXT wide format: one column per mass
- Long format: one row per time/mass/signal observation
- XLSX/XLS Excel files

For Excel:

```bash
python ms_deconvolution.py run.xlsx --sheet Data
```

For unusual long-format headings:

```bash
python ms_deconvolution.py run.csv --format long ^
  --time-column "Cycle Time" --mass-column "AMU" --signal-column "SEM Current"
```

## Candidate chemistry for 2-butanone on Pt at 400 °C

The default library is intentionally hierarchical:

| Status | Component | Chemical rationale / useful channels |
|---|---|---|
| Enabled | 2-butanone | Reactant; m/z 43 is intense but nonspecific, so retain molecular ion m/z 72 |
| Enabled | H₂ + methyl vinyl ketone | Direct dehydrogenation pair; MVK has distinguishing m/z 55 and 70 |
| Enabled | CO + propane | Stoichiometric decarbonylation pair: C₄H₈O → C₃H₈ + CO |
| Enabled | Methane, ethane/ethylene, propane/propylene | C–C scission followed by surface hydrogen transfer or dehydrogenation |
| Enabled | H₂O and CO₂ | Track as background or secondary oxidative products; they are not automatically primary products in oxygen-free N₂ |
| Disabled | 2-butanol | Hydrogenation or transfer-hydrogenation product; less likely without a hydrogen source |
| Disabled | Acetone and acetaldehyde | Possible secondary C–C redistribution/scission oxygenates |
| Disabled | Butenes and butane | Possible secondary C₄ deoxygenation products |
| Disabled | N₂, O₂, and Ar | Carrier/background components; enable only when deliberately modeling their drift |

The disabled candidates remain in `analytes.csv` so they can be enabled after their distinguishing ions or standards support them. Enabling every chemically possible species at once is not recommended: the patterns become too correlated for a stable solution.

Carbon deposited on the catalyst is another plausible high-temperature sink, but it is not directly observable by the online MS. A carbon balance or post-run catalyst characterization is needed to assess it.

## Important overlaps

### m/z 28

N₂, CO, ethylene, ethane, propane, and a CO₂ fragment all contribute at m/z 28.

- Do not quantify CO or C₂ products from m/z 28 alone.
- Retain m/z 26, 27, 29, and 30.
- A carrier-only baseline helps but does not replace calibration.

### m/z 43

2-butanone, propane, acetone, butane, and other fragments can contribute at m/z 43.

- Use m/z 72 with companion 2-butanone fragments for the reactant.
- Use m/z 70 and 55 for methyl vinyl ketone.
- Use m/z 44 with the propane fragment pattern, while accounting for CO₂.

### Large m/z 40 background

The short example export has an exceptionally large m/z 40 signal with m/z 36/38 companions, consistent with argon. This may be a carrier from that particular experiment, a residual background, or a different run than the intended N₂ experiment.

The code warns when m/z 40 dominates. To remove only the main Ar channel:

```bash
python ms_deconvolution.py run.txt --exclude-masses 40
```

If m/z 35–39 are also dominated by Ar isotope/tailing/background signal, exclude the range after inspecting those channels:

```bash
python ms_deconvolution.py run.txt --exclude-masses 35-40
```

Or down-weight noisy background channels using the carrier-only baseline:

```bash
python ms_deconvolution.py run.txt --noise-weighting
```

## Choose the baseline correctly

By default, the median of the first 5% of scans is subtracted (at least 3 and at most 100 points). These scans must represent an appropriate pre-reaction or carrier-only period.

Specify an elapsed-time range in seconds:

```bash
python ms_deconvolution.py run.txt --baseline-start 0 --baseline-end 300
```

Or specify a number of initial scans:

```bash
python ms_deconvolution.py run.txt --baseline-points 50
```

For data already baseline-corrected:

```bash
python ms_deconvolution.py run.txt --baseline none
```

Negative values after subtraction are clipped to zero for NNLS. The original signed values remain in `reconstructed_spectra_long.csv`.

## Calibrate before reporting concentrations or selectivity

Edit `analytes.csv`:

| Column | Meaning |
|---|---|
| `species` | Component name |
| `mz` | Nominal mass |
| `relative_intensity` | Fragment intensity on a consistent scale |
| `sensitivity_signal_per_ppm` | Optional fitted coefficient per ppm |
| `enabled` | `true` or `false` |
| `notes` | Role and limitations |

For defensible concentration or selectivity:

1. Record carrier-only baseline at the same flow, pressure, temperature, and capillary settings.
2. Record pure 2-butanone feed through an inert bed or bypass.
3. Run known mixtures or standards for likely gases.
4. Replace starter fragment ratios with those observed on this instrument.
5. Fit sensitivity against known ppm and enter the slope in `sensitivity_signal_per_ppm`.

`relative_pct_*` is the share of the **fitted MS signal**, not mol%, conversion, carbon selectivity, or product selectivity. True selectivity additionally requires response-factor correction, molar flows, reactant conversion, and carbon-number accounting.

## Useful options

Use selected masses:

```bash
python ms_deconvolution.py run.txt --include-masses 2,12-18,25-30,41-45,55-58,69-74
```

Smooth the fitted time profiles with an 11-point Savitzky–Golay filter:

```bash
python ms_deconvolution.py run.txt --smooth-window 11
```

Extract two unlabeled patterns from positive residuals:

```bash
python ms_deconvolution.py run.txt --unknown-components 2
```

Unknown components are exploratory NMF factors, not chemical identifications.

## Test the installation

Generate and fit a synthetic 2-butanone experiment:

```bash
python ms_deconvolution.py --make-demo demo.csv
python ms_deconvolution.py demo.csv --output demo_results
```

Run the MKS-format parser regression test:

```bash
python -m unittest discover -s tests -v
```

## Output files

| File | Contents |
|---|---|
| `deconvolved_components.csv` | Timestamp, scan, fitted signals, relative percentages, optional ppm, R², RMSE |
| `component_overview.png` | Main time-profile plots |
| `fit_quality.png` | Intense corrected channels and spectrum R² |
| `residual_heatmap.png` | Missing or poorly fitted signal by time and m/z |
| `reconstructed_spectra_long.csv` | Original signed signal, corrected signal, fit, and residual |
| `component_similarity.csv` | Pairwise spectral overlap |
| `baseline_and_noise.csv` | Subtracted baseline and baseline noise |
| `analysis_report.txt/.json` | Parsing, fitting, and warning summary |

Relative percentages are set to zero until the total fitted signal exceeds three times the combined baseline-noise norm.

## Reference-data notes

The starter EI patterns are for initial model construction, not publication-ready calibration. The 2-butanone and methyl-vinyl-ketone identities and spectra were checked against open MassBank records and the NIST Chemistry WebBook. NIST restricts redistribution of some downloadable spectra, so no NIST spectral file is bundled.

- 2-butanone MassBank record: <https://massbank.eu/MassBank/RecordDisplay?id=MSBNK-Fac_Eng_Univ_Tokyo-JP001803>
- Methyl vinyl ketone MassBank record: <https://massbank.eu/MassBank/RecordDisplay?id=MSBNK-Fac_Eng_Univ_Tokyo-JP002135>
- NIST Chemistry WebBook: <https://webbook.nist.gov/>
- MassBank open-data repository: <https://github.com/MassBank/MassBank-data>
