# Results and execution

The study is written in English. MATLAB code, raw-data preparation,
reconstruction, and prediction evaluation were run for Set1, Set2, and Set3.
The production evaluation software and all raw files were left unchanged.

Set1 contains 7,476 operating points, Set2 8,092, and the newly generated Set3
evaluation 9,022. Set3 uses a different machine configuration (four rather than
three pole pairs). Nominal excitation currents define the groups; calculations
use measured currents.

Joint RMSE of both flux components for completely omitted interior excitation
current groups, in mVs:

| Set | A: existing | C: joint constant fit | E: q-/excitation-dependent | F: two constant axis factors | G: three currents | H: three-current reference |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 5.998 | 5.518 | 5.453 | 5.490 | 5.513 | 1.787 |
| 2 | 6.309 | 5.736 | 5.470 | 5.653 | 5.443 | 1.780 |
| 3 | 6.970 | 6.996 | 6.839 | 6.529 | 6.425 | 1.480 |

Each row uses the same supported physical query points for all eight methods.
Coverage is 89.15%, 89.03%, and 90.75% of the interior groups, respectively.
Boundary groups are reported separately in `results/boundary_groups.csv`.
The joint metric is the square root of the mean squared error across both
components; it is not an arithmetic mean of individual group RMSEs.

Recommendation: If reduction to two map axes is required and both axes may be
shifted, F is the strongest simple compromise. Its two constant parameters
improve joint error over A by approximately 8.5%, 10.4%, and 6.3% across the three
datasets. E is a useful alternative for Set1/2 when the physical q-axis must
remain unchanged. G does not consistently outperform simpler variants; its small
additional improvement does not generally justify the extra parameters and
invertibility conditions. If accuracy takes priority, H is preferable:
approximately 21 times as many map nodes, but substantially smaller prediction errors.

Exact constant d-reduction requires both `Lde = k Ldd` and `Lqe = k Lqd`.
The derivative diagnostics retain residuals, and the shifted original surfaces
sometimes require different flux values at the same reduced position. This does
not establish an isolated physical error floor: measurement, temperature,
preprocessing, and interpolation remain included. Varying iq alone does not
explain the results.

Limitations: There is no application-based acceptance threshold or independent
absolute flux calibration. Set3 lacks phase/RMS metadata, among other channels;
its rotor-temperature channel cannot be used as a calibrated temperature trace.
Four G fits reached the finite search budget; no global optimum is claimed.
Extrapolation of the interpolants was disabled. Supported reduced queries for a
boundary group must still be distinguished from interpolation within the
measured three-current domain.

The small analytic controls produce the expected results: constant and integrated
variable reduction agree up to rounding; the incorrect product of excitation
current and local factor differs; the deliberately incompatible case has
different q-fluxes at identical reduced coordinates. No test harness or
assertions were added.

Serial reproduction in PowerShell from the project directory, using the current
`matlab/k_factor/` layout:

```powershell
Set-Location C:/Git/Winding_Factor
matlab -batch "restoredefaultpath; addpath('C:/Git/Winding_Factor/matlab/k_factor'); run_k_factor_analysis"
Set-Location C:/Git/Winding_Factor/winding_factor_study
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error atlas.tex
pdflatex -interaction=nonstopmode -halt-on-error atlas.tex
```

In the reported run, the three independent datasets were processed concurrently
in separate MATLAB processes with `maxNumCompThreads(2)` and
`analyze_flux_map_reduction(1)`, `(2)`, and `(3)`, respectively. The corrected runs
took approximately 1,000, 1,226, and 943 seconds. Logs are in
`results/reduction_set*.log`. MATLAB version: R2025b. External citations and the
bibliography were removed from the paper at the user's request.

Data provenance, detailed outputs, and the mapping of review requirements are
in [RUN.md](RUN.md). Complete individual-group comparisons are in the atlas;
numerical MAT/CSV files are also supplied as a results archive.
Axis labels use "in A" and "in mVs". The 3D views of ψd exchange the displayed
d/q axes; heatmaps retain their existing arrangement. Figure 1 uses open meshes
to keep individual excitation levels visible. A shared LaTeX template generates
atlas pages by looping over three dataset definitions instead of repeating page blocks.
The work remains on `codex/flux-reduction-study`; no merge or deployment was performed.
