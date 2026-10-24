# Flux function models: execution and provenance

English paper and figures; German handover. Existing k-factor paper and production
controller files stay unchanged. Branch: `codex/flux-function-models`.
Final outputs: 11 manuscript pages and 377 atlas pages (376 group/model figures).
Both PDFs compiled successfully and were rendered for visual inspection; no unresolved
references or overfull boxes remain.
The supplied `Flusskennfelder_Funktionsmodelle_Entwurf.pdf` is conceptual background;
its schematic figures and proposed criteria are not empirical results for this study.

Implementation sequence:
1. Reuse the fresh grouped Set1--3 measured-current/voltage-derived flux exports.
   Inspect the optional processed temperature fixture and the controller signal paths.
2. Fit shared Gaussian bases and cubic tensor splines using group-weighted linear
   least squares. Select size, width and regularization on separate validation data.
3. Refit selected models on training plus validation, then compare identical physical
   test points with a newly fitted full LUT and existing k variants A--G.
4. Check analytical derivatives numerically, report their disagreement with local
   training-data slope estimates, and time the same two-flux/six-derivative package.
5. Generate per-group surfaces/error maps, summary tables, English paper and PDF.

No test framework, production integration or new k family. The few numerical checks
are direct calculations. No external bibliography, following the user's preference.

Start from this repository in MATLAB R2025b:
```matlab
run_flux_function_study
```
Individual comparison scripts load saved models without refitting.

## Inputs and reused calculations

- `prepare_reduction_data.m`: raw acquisition channels, arithmetic repeat means,
  measured currents, and stationary voltage-derived fluxes. All three fresh compact
  inputs are `Resources/DataSets/Set*/extracted/evaluation/research_flux_measurements.mat`.
  Their `Total.<speed>.<excitation>.Meas.PE_psi_d_meas/PE_psi_q_meas` fields precede
  the interpolation used here. These are not smoothed shared LUT values.
- `prepare_flux_fit_data.m`: reads these exports, records requested speed and stator
  temperature ranges, and makes the common split. Nominal excitation defines groups;
  pointwise measured `I_exc_meas` is the third function argument. Temperatures were
  not controlled, so these are not experimentally isothermal current maps.
- Optional temperature reference:
  `Resources/GUI_VICE-MeasurementEvalKit/SourceCode/tests/Artifacts/json_api_headless_real_run/psm_temperature/flux_psm_temperature_30C_70C.mat_result.mat`,
  `Total.pos_2500_rpm.pos_0_A.Meas`. Its run log refers to an absent acquisition fixture.
  Raw voltage/current fields are retained in a processed artifact, whose physical
  provenance cannot be independently established. The coordinate is rotor temperature
  in **degrees C**, while temperature derivatives can be expressed per K.
  Both temperature bands remain in training; the holdout is by stator-current tile.
- The current evaluation kit's `flux_calculation.m` supplies the measured-voltage
  stationary convention. `fit_calculation.m` fits flux surfaces and differentiates
  them; those older fitted derivatives are not imported as independent measurements.
  `validate_winding_factor.m` likewise provides training-correlated local-regression
  derivative diagnostics. The new comparison recomputes local slopes from its own
  final training rows using 40/80/120 neighbors.
- `fit_reduction_map.m`, `reduction_coordinates.m`, `reduction_lookup.m` are reused
  unchanged for normalized source triangulation, raster storage and physical queries.
  `fit_flux_k_models.m` is the original `fit_methods`/`objective` from
  `analyze_flux_map_reduction.m`, with only the entry-point name changed.
  It retains A--G, the internal group fitting, optimizer budgets and G restriction.
- `eval_flux_model.m` differentiates the stored interpolants, applying the exact k
  chain rule. The explicit output indexing is `J(point,flux,input)`.
  `flux_bspline_basis.m` and `flux_tensor_basis.m` provide cubic values/partials.
- `plot_flux_function_study.m` reuses the blue-reference/orange-mesh, rotated d-flux
  view and physical signed-map layout from `plot_reduction_study.m`. The original
  plotting script is not called: it would regenerate the separate k-factor paper.
  The shared LUT mask can suppress otherwise finite coefficient predictions, including
  training boundary slices. Empty maps mean no common support, not an invented zero error.

## Controller inspection

Source read only: `Resources/VICE800_Software_Standard_003_004_000_002/a_OPEN_Subsystems/Control_P.slx`,
its `simulink/blockdiagram.xml`, path
`Control_P/Current_Controller/d/q-Stromregler_DCI/LUT1` (SID 2664).

- Inputs: `Id_haupt_soll` (2665), `Iq_haupt_soll` (2666), `Id_DCI_add` (2667).
  They come through tags `Id_Haupt_soll`, `Iq_Haupt_soll`, `I_d_DCI_add`, from
  current-controller request inlets. The local lookup-port names `Id_filt`/`Iq_filt`
  do not mean measured currents in this path.
- `LUT_Flux_MTPA/Psi_d_LUT` (2672) and `Psi_q_LUT` (2673) use main requests.
  `Add2` (2668) adds the DCI term before `LUT_Flux_SUM` (2682) and
  `LUT_Inductivities` (2695). Flux SUM SIDs are 2685/2686; Ldd is 2698.
- Outputs: two main-request fluxes, two summed-request fluxes and Ldd/Lqq/Ldq/Lqd
  at the summed request. No independent temperature/excitation input enters this
  local wrapper. The research package (two fluxes plus six partials at one point)
  is therefore not the entire existing wrapper package.
- This establishes the tap at the current-controller inlet before the local DCI
  addition. It does not establish that an unchanged upstream torque request leaves
  every downstream request unchanged. No cache hit-rate trace or target task timing
  was acquired. Exact reuse needs all actual query arguments, parameters and mode;
  tolerance reuse would be a distinct approximation. Both could also apply to LUTs.

## Saved outputs and replay

All numeric files are under `results/flux_fits/` and are supplied in `model_results.zip`.
Extract the ZIP into that directory to evaluate/plot saved models without refitting.
Fresh rebuilding additionally needs the three compact inputs above; the optional
temperature case needs its original processed artifact. Replay already saved data
if that optional source is absent; do not reinterpret it as a missing machine test.

- `data.mat`: values, physical inputs, active units/conditions, split masks.
- `rbf_models.mat`, `bspline_models.mat`, `baseline_models.mat`: complete evaluators,
  basis/scaling/bounds, selected coefficients and training identities. RBF centers
  are shared; spline coefficients use first-axis-fastest MATLAB tensor ordering.
- `comparison.mat`: point predictions, Jacobians, support masks and local slopes.
- `metrics.csv`: modes 1 reconstruction, 2 independent inner test, 3 upper boundary;
  group 0 pools points, positive groups are source-group indices.
- Methods 1--10 are RBF, Spline, LUT, kA, kB, kC, kD, kE, kF, kG.
  In `derivatives.csv`, reference methods 11/12 mean 40/120 versus 80 neighbors.
  First three sets use mH throughout; dataset 4's third-input differences use mVs/K.
- `differences_from_LUT.csv`: fresh training-only LUT comparison;
  `existing_LUT_differences.csv`: pre-existing all-data LUT comparison, which is not
  independent validation. The latter comes from the previous study's archived
  `flux_map_reduction.mat`; the script extracts that archive if needed.
- `model_sizes.csv`: numeric online arrays versus research MATLAB structure bytes.
  `runtime.csv`: seven batch-average measurements of the same scalar output package.
  Neither includes a target-hardware demonstration or an additional geometric guard.
- `numerical_checks.csv`, `lookup_consistency.csv`: centered differences, saved
  training-index overlap, and comparison with the original LUT evaluator. Direct
  checks only; no test framework or assertions.

Saved-model replay from the repository root:
```matlab
compare_flux_fits
compare_fits_with_k_factor
plot_flux_function_study
```
Compile the main paper from `function_model_study/manuscript` and the atlas from
`function_model_study`, running each command twice for cross-references:
```powershell
Set-Location C:/Git/Winding_Factor/function_model_study/manuscript
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
Set-Location ..
pdflatex -interaction=nonstopmode -halt-on-error atlas.tex
pdflatex -interaction=nonstopmode -halt-on-error atlas.tex
```
