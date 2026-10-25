# Winding_Factor

Two MATLAB studies of stator flux and coordinate reduction in an EESM.
The Python version is a learning project. Meas/raw data preparation and the
conventional contour k-factor estimator are implemented; other analyses remain
documented placeholders.

```text
run_matlab_analysis.m       Run both MATLAB studies
matlab/k_factor/           k-factor and coordinate-reduction study
matlab/function_model/     Function models, fits, and their comparisons
matlab/common/             Shared preparation, geometry, and path helpers
run_python_analysis.py     Run the implemented Python contour study
python/k_factor/           Importer, contour estimator, other analysis placeholders
python/function_model/     Function-model module and study-entry placeholders
Resources/                 Shared measurement data and MeasEval resources
results/                   Generated MAT/CSV results (ignored by Git)
winding_factor_study/      First study: manuscripts, figures, results archive
function_model_study/      Second study: manuscripts, figures, results archive
README.md, AGENTS.md        Project description and working guidelines
```

MATLAB R2025b with base MATLAB is the environment used for the existing studies.
Parameters, data selection, group splits, interpolation methods, and error
definitions are preserved. Raw data and study folders remain in their existing locations.

## Running MATLAB

All three options work from any working directory. In a fresh MATLAB session,
choose the required entry point:

```matlab
project='/absolute/path/to/Winding_Factor';

% Only the k-factor study:
addpath(fullfile(project,'matlab','k_factor'));
run_k_factor_analysis

% Only the function-model study:
addpath(fullfile(project,'matlab','function_model'));
run_function_model_analysis

% Both studies, k-factor first:
addpath(project);
run_matlab_analysis
```

Each study resolves its own common helpers and MeasEval paths. The root entry
point calls only the two study entry points; it contains no calculation logic.
`run_flux_function_study` remains a compatibility alias for
`run_function_model_analysis` in `matlab/function_model/`.

The k-factor study runs `prepare_reduction_data`, `analyze_winding_factors`,
`validate_winding_factor`, `audit_reduction_sources`, `analyze_flux_map_reduction`,
and `plot_reduction_study`, in that order.
The function-model study prepares missing grouped measurements only if needed,
then runs `prepare_flux_fit_data`, `fit_flux_rbf`, `fit_flux_bspline`,
`fit_flux_baselines`, `compare_flux_fits`, `compare_fits_with_k_factor`, and
`plot_flux_function_study`. The plot step also calls `write_flux_function_tables`.

The function-model comparison uses its own training/validation fits for the LUT
and k models A–G. These are different from the k-factor study's all-data and
holdout fits and cannot replace them. For the additional all-data reference
comparison, it reuses `results/flux_map_reduction.mat`; if missing, it imports
only that file from `winding_factor_study/research_results.zip`, preserving other
results. If neither source exists, the error gives the exact path and commands
for `prepare_reduction_data; analyze_flux_map_reduction`. It never starts the
other study's audit or plots. The combined run reuses the freshly generated
reference and grouped inputs without repeating their preparation or fitting.

All regular figures, including atlas figures, are created visibly and remain
open after export. Each call creates new figures without closing previous ones.
Temporary font/interpreter defaults are restored. LaTeX compilation remains separate.

```sh
matlab -batch "addpath('/absolute/path/to/Winding_Factor'); run_matlab_analysis"
```

For a single-study batch run, add its folder and call its entry point instead.
The batch commands produce the same exports; windows close when MATLAB exits.

## File responsibilities

| Folder | Files and responsibilities |
| --- | --- |
| `matlab/k_factor/` | `run_k_factor_analysis`, `analyze_winding_factors`, `validate_winding_factor`, `audit_reduction_sources`, `analyze_flux_map_reduction`, `plot_reduction_study` |
| `matlab/function_model/` | `run_function_model_analysis` (legacy alias `run_flux_function_study`), `prepare_flux_fit_data`, `fit_flux_rbf`, `fit_flux_bspline`, `fit_flux_baselines`, `compare_flux_fits`, `compare_fits_with_k_factor`, `plot_flux_function_study`, `write_flux_function_tables` |
| `matlab/function_model/` helpers | `fit_flux_k_models` fits this study's k baselines; `eval_flux_model`, `flux_bspline_basis`, `flux_tensor_basis`, `flux_group_loss`, and `flux_model_size` support function-model fitting and evaluation |
| `matlab/common/` | `prepare_reduction_data` is used by both studies; `estimate_contour_factor`, `fit_reduction_map`, `reduction_coordinates`, and `reduction_lookup` supply shared numerical operations; `setup_matlab_paths` and `require_analysis_file` handle paths and missing inputs |
| `python/k_factor/` | Named Meas/raw importer, conventional contour estimator and study starter; derivative analysis, source audit, other coordinate models and plots remain placeholders |
| `python/function_model/` | Placeholders for data splits, RBF/spline/baseline fits, both comparisons, plots/tables, and `run_function_model_analysis.py` |

## Individual entry points

Add the folder containing the desired entry point, then call it directly.
For example, in a fresh session:

```matlab
addpath('/absolute/path/to/Winding_Factor/matlab/function_model');
fit_flux_rbf   % Requires saved data.mat; the error names its preparation step.
```

For shared preparation alone, add `matlab/common/` and call
`prepare_reduction_data`. For a k-factor step, add `matlab/k_factor/` and call,
for example, `analyze_flux_map_reduction(2)`. Each analysis entry point bootstraps
its common dependencies, has its own function workspace, and derives resource
and output paths from its source location. Missing inputs are reported with
full paths and concrete preparation commands.

| Entry point | Required inputs / preparation |
| --- | --- |
| `prepare_reduction_data` | Original Set1–3 files in `Resources/DataSets/Set*/raw/` |
| `analyze_winding_factors`, `validate_winding_factor` | `results/reduction_data.mat` from `prepare_reduction_data` |
| `audit_reduction_sources` | Reduction data and complete historical FluxEval files; see below |
| `analyze_flux_map_reduction` | Reduction data; optional `(1)`, `(2)`, or `(3)` selects one set; collects all completed saved sets |
| `plot_reduction_study` | `reduction_data.mat`, `winding_validation.mat`, `flux_map_reduction.mat` |
| `prepare_flux_fit_data` | Three `extracted/evaluation/research_flux_measurements.mat` files from `prepare_reduction_data` |
| `fit_flux_rbf`, `fit_flux_bspline`, `fit_flux_baselines` | `results/flux_fits/data.mat` from `prepare_flux_fit_data` |
| `compare_flux_fits` | `data.mat`, `rbf_models.mat`, `bspline_models.mat`, `baseline_models.mat` in `results/flux_fits/` |
| `compare_fits_with_k_factor` | `data.mat`, `comparison.mat`, and the A–H reference; extracts the winding study archive if the reference is missing |
| `plot_flux_function_study` | Function-model data, comparison results, candidate CSVs from both fits, and comparison tables |
| `write_flux_function_tables` | Saved function-model data and CSVs from both comparisons |
| `run_k_factor_analysis` | Original raw data and complete historical audit sources; runs the first study |
| `run_function_model_analysis` / `run_flux_function_study` | Grouped measurement data, or raw data to generate missing inputs; saved all-data reference or supplied archive |

Helpers requiring arguments are `fit_reduction_map`, `reduction_coordinates`,
`reduction_lookup`, `estimate_contour_factor`, `fit_flux_k_models`, `eval_flux_model`,
`flux_bspline_basis`, `flux_tensor_basis`, `flux_group_loss`, and `flux_model_size`.
They are not additional analysis steps. `setup_matlab_paths` and
`require_analysis_file` only handle paths and missing-input messages.
The historical extraction scripts and `Resources/calculate_winding_factor.m`
are not part of either study's regular workflow.

## Resources and results

The historical audit additionally requires these locally held files, which are
not included in the repository:

- `Resources/DataSets/Set1/FluxEval/Matlab_200_201_202_FluxEval.mat`
- `Resources/DataSets/Set2/FluxEval/eval_flux_eesm_excitation_0A_to_40A.mat`

Compact Meas exports cannot replace their saved LUT/Calc coordinates.
Without these originals, the shared entry point stops at the audit and reports
the missing path. The second study can run separately using the existing grouped data.
The optional temperature artifact is included only if it exists at
`Resources/GUI_VICE-MeasurementEvalKit/SourceCode/tests/Artifacts/json_api_headless_real_run/psm_temperature/flux_psm_temperature_30C_70C.mat_result.mat`;
otherwise Set1–3 are processed as before. It is not independent temperature
validation. Replaying the archive retains the previously saved fourth case.

MAT/CSV outputs remain in `results/` and `results/flux_fits/`;
fresh grouped data remain in `Resources/DataSets/Set*/extracted/evaluation/`.
PNG exports remain in `<study>/figures/` and `figures/atlas/`; generated LaTeX
tables and atlas lists remain in the corresponding study folders.
Saved models can be reused without refitting (extraction replaces results with the same names):

```matlab
project='/absolute/path/to/Winding_Factor';
unzip(fullfile(project,'winding_factor_study','research_results.zip'),fullfile(project,'results'));
unzip(fullfile(project,'function_model_study','model_results.zip'),fullfile(project,'results','flux_fits'));
addpath(fullfile(project,'matlab','function_model'));
addpath(fullfile(project,'matlab','k_factor'));
compare_flux_fits
compare_fits_with_k_factor
plot_flux_function_study
% First study: analyze_winding_factors; plot_reduction_study
```

MeasEval remains unchanged in [Resources/MeasEval](Resources/MeasEval/README.md).
The actual Python classes are `PointSet` (point geometry; NumPy, hashlib),
`SampledField` (field values on a `PointSet`; NumPy), and `EdgeIntersector`
(SciPy Delaunay triangulation, filtered edges, and level intersections; both model
classes, NumPy, SciPy). `EdgeIntersections` is the additional result dataclass in
the same module. The optional demo uses Matplotlib and `himmelblau`.
These resources are existing library code.

The MATLAB studies directly use `delaunay`, `delaunayTriangulation`,
`pointLocation`, `scatteredInterpolant`, and `griddedInterpolant`; they do not call
these Python classes. The four supplied MATLAB geometry functions
`backend.numeric.geometry.{unique_edges,valid_edges,edge_intersections,triangle_intersections}`
are supplementary resources. Each study adds only its own folder, `matlab/common/`, and
`Resources/MeasEval/matlab/`. The combined entry point adds both study folders.
No GUI or repository paths are added recursively.

Manuscripts and further provenance: [winding study](winding_factor_study/RUN.md),
[function models](function_model_study/RUN.md). Run LaTeX twice in the respective
manuscript/atlas directories listed there. MAT/JSON data: [DataSets](Resources/DataSets/README.md).

## Python development

The local `.venv` uses Python 3.14, NumPy and SciPy and is ignored by Git.
From the project root:

```sh
source .venv/bin/activate
python python/k_factor/prepare_reduction_data.py
```

For a new checkout: `python3.14 -m venv .venv`, activate it as above, then
`python -m pip install numpy scipy`.

`ReductionDataPreparation` reads one grouped Meas export or raw measurement
file using the existing `.json.gz` copies. It detects the format automatically;
the default is Set1's `extracted/evaluation/research_flux_measurements.json.gz`.
Grouped `Total.<speed>.<excitation>.Meas` data use saved fluxes directly.
Raw data are averaged and their stationary fluxes calculated as before.

`run()`, `from_meas()` and `from_raw()` return dictionaries with the same twelve
names. Each value is a one-dimensional floating-point NumPy array of length N;
the same index identifies the same measurement point in every array.

| Name | Unit | Meas channel / block field | Raw source |
| --- | --- | --- | --- |
| `id_meas` | A | `PE_Id_meas_avg` | same channel |
| `iq_meas` | A | `PE_Iq_meas_avg` | same channel |
| `ie_meas` | A | `I_exc_meas` | same channel |
| `rpm_meas` | rpm, mechanical | `CAN_RPM_meas` | same channel; no conversion |
| `id_req` | A | `PE_Id_req_meas` | same channel |
| `iq_req` | A | `PE_Iq_req_meas` | same channel |
| `ie_req` | A | `I_exc_req` | same channel |
| `rpm_req` | rpm, mechanical | `PE_Rpm_req` | same channel |
| `psi_d` | Vs | `PE_psi_d_meas` | `(PE_Uq_meas_avg - Const_Machine_R_s * iq_meas) / PE_w_el_rad_avg` |
| `psi_q` | Vs | `PE_psi_q_meas` | `(Const_Machine_R_s * id_meas - PE_Ud_meas_avg) / PE_w_el_rad_avg` |
| `ie_ref` | A | block `Ref_exc_current` | nominal `I_exc_req` group key |
| `rpm_ref` | rpm, mechanical | block `Ref_Speed` | nominal `PE_Rpm_req` group key |

Meas vectors are copied point by point; only scalar block references are repeated.
Raw samples retain arithmetic averaging by recorded speed request, excitation
request and `PE_Index_meas`, followed by stationary flux calculation. The raw
references come from those grouping keys, matching
`matlab/common/prepare_reduction_data.m`; requested and measured arrays remain
separate. Raw preprocessing still excludes index <= 0, nonfinite measured
currents/flux inputs and averaged electrical angular speeds with absolute value
<= 1 rad/s. Missing optional channels do not exclude samples.

`run()` applies one common mask to every array: points with nonfinite `id_meas`,
`iq_meas`, `ie_meas`, `psi_d` or `psi_q` are removed. Those five quantities are
required. Missing additional quantities become arrays of `NaN` and are reported
in one warning per import. Missing requested values are never filled from group
references. Raw data additionally require measured dq voltages, electrical
angular speed, stator resistance and the measurement index to compute fluxes and
average samples. If raw grouping keys are unknown, affected records are retained
individually rather than averaging potentially unrelated points; unknown
references remain `NaN`. Channel length errors identify the source block/channel.

Within each Meas speed group, complete excitation blocks are sorted numerically
and stably by `Ref_exc_current`, following `Resources/calculate_winding_factor.m`.
Equal-reference blocks stay separate, unknown references follow known ones, and
point order within each block is retained. Raw groups are also sorted by nominal
excitation within each speed, retaining speed input order and operating points
in first-occurrence order. Every named array follows the same ordering.
Measured currents remain the physical map coordinates.

```python
from python.k_factor.prepare_reduction_data import ReductionDataPreparation

preparation = ReductionDataPreparation()
data = preparation.run()
mask = data["ie_req"] == 6.0
id_meas = data["id_meas"][mask]
iq_meas = data["iq_meas"][mask]
psi_d = data["psi_d"][mask]
```

`preparation.save(data)` writes the common `results/python/reduction_data.npz`;
the standalone command does both. Downstream code uses
`data = numpy.load("results/python/reduction_data.npz")` and accesses the same
fields regardless of the source. Pass a raw or evaluated JSON path to
`ReductionDataPreparation(filename)` to choose another input. Historical Meas
exports can differ numerically from freshly prepared raw data; normalization
preserves their existing values.

The conventional contour analysis is implemented in
`python/k_factor/analyze_winding_factors.py`. These entry points run that analysis
only, defaulting to the three compact research exports:

```sh
.venv/bin/python python/k_factor/analyze_winding_factors.py
.venv/bin/python python/k_factor/run_k_factor_analysis.py
.venv/bin/python run_python_analysis.py
```

All accept optional measurement `.json.gz` paths, for example:

```sh
.venv/bin/python run_python_analysis.py Resources/DataSets/Set1/extracted/evaluation/research_flux_measurements.json.gz
```

For inspecting intermediate arrays without printing a summary:

```python
from python.k_factor.analyze_winding_factors import analyze_winding_factors

data = ReductionDataPreparation().run()
results = analyze_winding_factors(data)
contours = results[3000.0]  # keyed by nominal mechanical rpm
print(contours["k"], contours["mean_ie"], contours["exclusion_reason"])
```

Each speed is processed separately, with excitation groups sorted numerically
by `ie_ref`. `estimate_contour_factor(data)` handles one known speed. Missing
`rpm_ref` prevents dataset grouping; missing `ie_ref` prevents estimation for
that speed. No unknown groups are pooled or skipped. Results include `k`,
`status`, `reason`, sorted `ie_ref`, 100 `levels`, (group, level)
`intersection_counts`, `mean_id`, `mean_iq`, `mean_ie`, and (adjacent pair, level)
`delta_ie`, `ratios`, `valid`, `exclusion_reason`. `edges` and `group_errors`
make geometry failures inspectable. An undeterminable factor is `NaN` with
`status="not_determinable"`; a failed group triangulation blocks the whole speed.

The reference is `matlab/common/estimate_contour_factor.m` and section 3.2 of
`winding_factor_study/main.tex`. Both specify the same estimator. Implementation
details supplied by MATLAB are: inclusive 100 levels at 25--75% of the pooled
flux **value range**, unique undirected Delaunay edges with squared length
strictly below nine times the mean squared unique-edge length, inclusive
edge interpolation weights in [0,1], and equal weight per edge intersection.
Shared-vertex intersections remain repeated; constant-flux edges give no unique
intersection. The existing `batched_edge_intersections` kernel is used with
`numeric_policy="matlab-ieee"`, keeping measured ie on the same weights as id/iq.
No fixed-iq condition, extra input averaging, coordinate normalization or smoothing
is added. SciPy Delaunay can choose different diagonals from MATLAB in ambiguous
geometries; no triangulation options are tuned to match a final factor.

Adjacent factors are `abs(delta_mean_id) / delta_mean_ie`. Only nonfinite factors
and `abs(factor) < eps` are excluded, followed by one arithmetic mean over all
remaining pair/level entries. MATLAB imposes no small-denominator cutoff and
retains negative denominators. Missing contours leave NaN means; neighboring
pairs are never replaced by nonadjacent pairs. The MATLAB centroid-flux-error
diagnostic does not affect k and is not computed here.

`Resources/calculate_winding_factor.m` confirms nominal group ordering and
channel mapping, but currently ends after a measurement plot; it contains no
factor calculation. `Resources/MeasEval/matlab/.../unique_edges.m` and
`edge_intersections.m` match the estimator's edge/intersection policy.
`valid_edges.m` instead weights edge occurrences per simplex when computing
its threshold and is not the policy used here. No `Torque_LUT_Calc.m` is
present in this checkout. Importing any Python entry point starts no analysis.
Other k models, function models, fits and plots remain unimplemented.
