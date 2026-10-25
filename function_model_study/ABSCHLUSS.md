# Function-model study results

The separate MATLAB evaluations for Set1–3 were completed. Gaussian RBFs and
cubic tensor-product B-splines share basis values for both fluxes. Comparisons
use a full LUT rebuilt from the same training data and the existing k variants
A–G. The previous k-factor paper and production control software were not changed.

Joint flux RMSE on completely held-out interior excitation-current groups, in mVs:

| Set | RBF | Spline | Full LUT | kA | kF | Common points |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 3.074 | 1.450 | 1.744 | 4.521 | 4.570 | 1,034 / 1,156 |
| 2 | 3.200 | 1.552 | 1.823 | 5.726 | 5.070 | 1,551 / 1,734 |
| 3 | 1.853 | 1.210 | 1.421 | 6.115 | 5.706 | 1,654 / 1,819 |

All ten methods are evaluated on the same supported physical points within each
set. Common coverage is approximately 89.45%, 89.45%, and 90.93%. Overall errors
are pooled from squared point errors and point counts. The CSVs additionally
contain d/q RMSE, separate 95% quantiles and maxima, all groups, reconstruction
errors, own support coverage, and differences from the fresh and existing LUTs.
The existing LUT had seen all groups and is explicitly not an independent
evaluation source.

For both fluxes together, the splines require 972 / 972 / 490 coefficients and
approximately 7.98 / 7.98 / 4.17 KiB of numerical online storage. The full LUT
requires approximately 2,154 / 2,154 / 1,539 KiB. The selected RBFs each have 128
centers, 264 flux coefficients including the affine part, and approximately
5.16 KiB of online data. Documentation, MATLAB structure overhead, workspace,
and an additional geometric validity guard are excluded from these online figures.

Recommendation: Prototype a local cubic spline evaluator on the target hardware
as the next accuracy/storage candidate, keeping the RBF as a runtime comparison.
The generic MATLAB RBF evaluator is faster here: approximately 14–27 µs per
complete query versus 69–137 µs for splines. Every method returns two fluxes,
six partial derivatives, and validity status. These are PC batch averages,
not individual-call timing bounds or evidence of 20 kHz operation.

The inspected existing control block processes requested currents and an
additional d-current component. It requires two flux evaluations at different
current points and four d/q inductances. Replacing it fully therefore involves
more than a single research query. Exact reuse requires unchanged complete
inputs, parameters, and modes; LUTs could also use this approach. No reuse rate
was measured.

There is no application-based threshold for flux/derivative errors or target
hardware demonstration including available task budget, required outputs,
numeric format, and validity checking. Independently measured derivatives are
unavailable. Local derivative references come from training flux data; their
disagreements are diagnostics. Temperatures were not constant in Set1–3.
The separately evaluated temperature test file contains 680 points, has no
available raw source, and includes only two temperature bands. Its 44 common
evaluation points concern spatially held-out current regions; they establish
neither physical temperature generalization nor a 4D demonstration. The RBF
performs best there, which also argues against a universal ranking.

The upper boundary group is reported separately and rejected outside the
training range. Even within the current box, the conservative common LUT mask
can hide boundary surfaces. Empty error maps indicate missing common support,
not zero error. The atlas contains all 376 group/method pages and a cover page.
The English paper has 11 pages; both PDFs were compiled, rendered, and visually inspected.

Direct numerical checks: the largest disagreement between analytical derivatives
and centered differences was approximately 2.22e-11 in the respective physical
units; LUT values agree with the existing helper to approximately 5.55e-17 Vs.
No saved training indices overlap test or boundary points. All k optimizers
reported positive exit flags in this run. This does not establish a global
optimum or independent physical validation.

Start from `C:/Git/Winding_Factor` using the current `matlab/function_model/` layout:

```powershell
matlab -batch "restoredefaultpath; maxNumCompThreads(2); addpath('C:/Git/Winding_Factor/matlab/function_model'); run_function_model_analysis"
```

Models, input points, and tables are archived in `model_results.zip`.
To evaluate without refitting, extract it into `results/flux_fits/` and run
`compare_flux_fits; compare_fits_with_k_factor; plot_flux_function_study`.
File references and complete LaTeX commands are in [RUN.md](RUN.md).
