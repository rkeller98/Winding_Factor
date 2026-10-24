# Research execution, 24 October 2026

Scope: stator flux linkages only. English manuscript, German handover.
Work is on `codex/flux-reduction-study`; production software and raw files are read-only.
The pre-existing manuscript/scripts were copied to `%TEMP%/winding_factor_before_20261024`
and preserved in baseline commit `fa08911`.

Execution order (all commands from `C:/Git/Winding_Factor`):
1. `prepare_reduction_data`: average acquisition repetitions by file and operating-point index,
   retain measured currents, calculate stationary flux from measured dq voltages.
2. `analyze_winding_factors`, `validate_winding_factor`: conventional contours, local derivative
   diagnostics and three analytic numerical examples.
3. `audit_reduction_sources`: fresh vs historical preprocessing, original saved coordinates and constants.
4. `analyze_flux_map_reduction`: A--H, reconstruction, nested group holdout,
   controlled raster and scaling comparisons, generated tables and complete group figures.
5. `plot_reduction_study`: independent surfaces, folding/holdout Jacobian checks, boundary tables,
   all group figures and generated manuscript tables.
6. Compile `main.tex` and `atlas.tex`, render and inspect the PDFs.

Final build: 25 manuscript pages and 297 atlas pages. Both PDFs compiled twice;
all pages were rendered and visually surveyed, with enlarged checks of the changed figures.
No unresolved references or overfull boxes remain.

```powershell
matlab -batch "restoredefaultpath; addpath(pwd); prepare_reduction_data; analyze_winding_factors; validate_winding_factor; audit_reduction_sources; analyze_flux_map_reduction; plot_reduction_study"
```

Run LaTeX from `winding_factor_study` to resolve relative inputs; the exact final build
commands below use that working directory. Run twice for references:

```powershell
Set-Location C:/Git/Winding_Factor/winding_factor_study
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error atlas.tex
pdflatex -interaction=nonstopmode -halt-on-error atlas.tex
```

The three independent sets were fitted concurrently in three MATLAB processes, each with
`maxNumCompThreads(2); analyze_flux_map_reduction(1)` (respectively 2 and 3).
The last completed set collects all three completed result files. Serial execution with no argument
produces the same deterministic fitting procedure. Logs: `results/reduction_set1.log` through `set3.log`.

Source provenance:
- Original acquisition files remain under `Resources/DataSets/Set*/raw`.
- Fresh compact grouped flux exports: each set's `extracted/evaluation/research_flux_measurements.mat`.
- `results/reduction_data.mat` and `Set*_flux_points.csv` contain measured currents, recomputed
  fluxes, nominal labels, file/index identity and repetition counts. `data_inventory.csv` reports
  current tracking, excitation scatter, speed and recorded temperatures per group.
- The local corrected numerical source inspected was the measurement-evaluation checkout at
  `a66300e3`. `flux_calculation.m`, `fit_calculation.m` and `Torque_LUT_Calc.m` were read only.
  The research implementation uses the same measured-voltage stationary flux formula and current
  corrected output/input derivative convention. It deliberately rebuilds all sets with arithmetic
  averaging and no fresh outlier rejection; it does not claim to rerun the complete GUI pipeline.
- Set1/2 historical exports are only a provenance/coordinate control. Their flux difference is
  evaluated at matched physical coordinates and exported in `historical_preprocessing_comparison.csv`.
  These differences are not subtracted from reconstruction RMSE as a causal partition.
- Set1 and Set2 record different current sensor delays (0 and -24 microseconds). No delay causality
  is inferred. Set3 lacks machine-ID/phase-selector/delay/RMS comparison channels.
- `main 2.pdf` was not present in the project or supplied attachments. The supplied review requirements
  were therefore used directly. References to the handwritten working note were removed.

Settings are literal constants in the scripts: 81 nodes per reduced axis, 40/80/120 derivative neighbors,
at most 120 query records for each of three inner groups, equal component weights, regularization
`1e-8 Vs^2` on added coefficients, and bounded finite-budget parameter searches. Final outer-group
errors use every common supported point, not the inner fitting subsample. D--G use a nonzero search
coordinate origin so zero physical coefficients receive the same initial simplex step as nonzero ones.
Optimizer exit flags are retained; a budget stop is not a convergence claim.

Core results:
- `results/flux_map_reduction.mat`: all maps, physical-point predictions, masks and fold parameters.
- `results/winding_validation.mat`: derivative estimates and analytic controls.
- `method_metrics.csv`, `group_metrics.csv`, `boundary_groups.csv`: error/coverage/complexity.
- `surface_pairs.csv`, `projected_folding.csv`, `holdout_jacobians.csv`: independent geometry checks.
- `Set*_reduction.mat`: independent complete dataset results, including scaling/raster controls.
- Main figures and complete atlas figures are generated under `winding_factor_study/figures`.

Atlas maintenance: `atlas_pages.tex` contains three `AtlasDataset` calls with each set's excitation
list. The shared page layout and nested `foreach` loops are defined once in `atlas.tex`.
The d-flux 3D views exchange displayed d/q axes; heatmaps retain their d/q order.
All figure axes use the requested `quantity in unit` notation. These are presentation changes;
the fitted models, errors and numerical scaling are unchanged.
Figure 1 uses open meshes so the individual excitation levels remain visible.
`research_results.zip` contains the current numerical MAT/CSV outputs and execution logs;
extract it into `results/` to inspect the saved results without repeating the fits.

CSV conventions: method numbers 1--8 mean A--H; mode 1 is all-data reconstruction and mode 2
is interior whole-group prediction. The two endpoint groups are excluded from the all-model
interior mask, so their mode-2 common count in `group_metrics.csv` is deliberately zero; use
`boundary_groups.csv` for their actual prediction errors and support. H's own-support boundary
metrics must not be presented as a matched-support comparison with reduced models.

Review mapping:

| Requested issue | Implementation / manuscript location |
| --- | --- |
| Current source, Set3, measured vs nominal currents | preparation + source audit; data section |
| Contours, moving iq, centroid, quotient, secants | contour script + estimator section |
| Common exactness, local/global, determinant, LS | derivative script + exact-reduction section |
| Actual shared model and return reconstruction | fit/lookup helpers + experiment section |
| Independent transformed originals | plotting postprocessor + independent-surface section |
| A--H and integrated shift, two axes, bounded all-current model | comparison script + model section |
| No outer-group leakage | inner fitting inside each outer fold; prediction section |
| Identical support, signed errors, per-group coverage | metric masks/CSV; complete atlas |
| Raster and coordinate scaling | controls in comparison script; numerical-control results |
| English narrative, checked references, no working-note citation | rewritten main manuscript |
| Executed controls, compiled PDFs, recommendation, review fixes | generated outputs and final handover |

The detailed authorized specification is the user's attached research request. No test framework,
assertions, production GUI changes, scalar magnetic model or deployment are part of this work.
