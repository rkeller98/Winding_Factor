# Running and maintaining the k-factor study

The numerical reference is `matlab/k_factor/estimate_contour_factor.m`.
Its original 100 levels, arithmetic edge-intersection means, edge cutoff,
absolute numerator, signed denominator and aggregation are retained.
The 25 October 2026 cleanup reproduced all six contour fields and the saved
40/80/120-neighbour derivative studies exactly for Sets 1–3.

From the repository root in MATLAB R2025b:

```matlab
addpath(pwd);
% Only if canonical data need preparing again:
run_matlab_analysis([], 'raw', 'prepare');
addpath(fullfile(pwd,'matlab','k_factor'));
analyze_winding_factors;       % Cached canonical data; tabbed figures, no images
validate_winding_factor;       % Local derivatives; not required if cache exists
```

`results/reduction_data.mat` stores the canonical structs; `reduction_inputs`
extracts measured current/flux arrays and groups by `ie_ref`. For raw inputs,
`ie_ref` comes from the requested excitation; measured excitation is retained
for all numerical coordinates. Eval imports use their stored block references.
No copy of these MAT files is required by the follow-up analysis.

For a complete new A–H fitting run, use `run_k_factor_analysis`. It includes
`analyze_flux_map_reduction` and can be expensive. Do not run it merely to
inspect contour plots or reuse local derivatives. The historical source audit
`audit_reduction_sources` is separate and requires full original FluxEval inputs.
No function-model study is necessary for these k-factor investigations.

## Plotting and explicit exports

Each contour-analysis dataset opens one figure with tabs for measured, shifted,
reduced and reconstructed flux, a measured/reconstructed overlay, per-group
errors, contour statistics and the 3-D centroid diagnostic. The field tabs share
nominal-current checkboxes and `Alle`/`Keine`. The overlay uses a translucent
measured surface, reconstruction mesh and signed-error-coloured points.
Error tabs include density histograms and point statistics.

Image export is opt-in:

```matlab
analyze_winding_factors('png');
analyze_winding_factors('pdf');
% Rebuild only existing reduction displays/tables, without refitting:
plot_reduction_study;          % Many figures; no image exports by default
```

An explicit `analyze_winding_factors('png')` updates contour PNG copies used by
the manuscript. `plot_reduction_study('png')` regenerates all A–H atlas images
and tables; use it only when those data or displays actually need updating.
Individual tab layouts can be exported directly with `exportgraphics`.
Default MAT/CSV result writes are independent of optional image export.

## Existing results and atlas

- `results/winding_factors.mat`: contour factors and diagnostic data.
- `results/winding_validation.mat`: derivative estimates, masks and analytic examples.
- `results/flux_map_reduction.mat` and `Set*_reduction.mat`: existing A–H fits,
  reconstructions, holdouts and numerical controls.
- Existing `research_results.zip`: archived numerical outputs and logs; extract
  into `results/` when local caches are absent. Check provenance against current
  canonical data before combining archives from different runs.

`atlas_pages.tex` contains three excitation lists: 13, 14 and 10 groups.
The atlas remains 296 existing A–H comparison pages plus its introduction.
It uses 41×41 display queries and the previously stored models; it is distinct
from interactive contour tabs and the new local-factor plots. No atlas images
or group lists were regenerated during this cleanup.

Run LaTeX twice from this folder to resolve relative inputs and references:

```powershell
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error atlas.tex
pdflatex -interaction=nonstopmode -halt-on-error atlas.tex
```

The follow-up lives in [k_factor_investigation](../k_factor_investigation/README.md).
It changes no reference estimator and performs no new global optimization.
Local differential factors, contour departure, reconstruction RMSE, prediction
RMSE and derivative residuals remain different statistics. No symmetry or
causal error decomposition is assumed.
