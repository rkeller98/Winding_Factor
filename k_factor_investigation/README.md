# k-factor investigation

Follow-up to `winding_factor_study`, using the same canonical measurements and
saved derivatives. The conventional implementation remains the numerical reference.
This folder investigates identification geometry and local differential factors;
it introduces no global flux fit, RBF, spline or new parameter optimization.

## Run with existing results

From the repository root in MATLAB R2025b:

```matlab
addpath(fullfile(pwd,'k_factor_investigation','matlab'));
[averaging, figures] = investigate_contour_averaging;
local = investigate_local_factors;
fig = plot_local_factor(local{1});
plot_local_factor(local{1},local{1}.kd,'k_d',local{1}.valid_pair);
plot_local_factor(local{1},local{1}.kq,'k_q',local{1}.valid_pair);
```

Inputs are `../results/reduction_data.mat`, `winding_validation.mat` and
`flux_map_reduction.mat`. No MAT files are copied. The existing results archive
is described in `../winding_factor_study/RUN.md`; preparation and derivative
calculations remain in the original MATLAB folders.

Figures stay open. Local-factor plots have nominal-current checkboxes and
`Alle`/`Keine`, controlling each group's mesh and points together. Each group is
triangulated separately using its original points; triangles touching invalid
values or exceeding the existing normalized long-edge rule are removed.
Missing areas are not retriangulated or extrapolated. Axes are physical id, iq
and the local factor; colour is measured excitation. Negative supported factors
are retained. Rotation and zoom use MATLAB's ordinary figure controls.

There is no default image export. Examples of explicit exports:

```matlab
investigate_contour_averaging(1,[6 7],50,'png');
plot_local_factor(local{1},[],'k_{loc}',[],'pdf');
```

Exports go to `figures/`; a small representative set accompanies the manuscript.

## A: contour averaging

For a d-flux level lambda, the reference averages every retained edge
intersection in measured `(id,iq,ie)` coordinates. Its global factor is the
unweighted mean of `abs(delta mean(id))/delta mean(ie)` over adjacent nominal
groups and the original 100 levels in the middle 25–75% of the flux value range.
The denominator retains its sign. Additional display levels do not enter k.

`investigate_contour_averaging(set,group_indices,level_index,export_format)`
defaults to Set1, groups 6/7 (17/20 A), original level 50. A scalar group index
selects that group and its successor. Other levels and pairs use the same code.

The first figure overlays the available contour segments, retained intersections,
open endpoints and four representatives: arithmetic edge mean, exact physical
arclength centroid, componentwise median and componentwise 10%-per-tail trimmed
mean. Segment length is measured in physical `(id,iq)` coordinates; excitation
is interpolated along the segments. The length centroid integrates the segment
midpoints exactly, not a sampled approximation. It need not lie on the contour.

The second figure compares representatives and fixed-iq roots. Equal-arclength
and equal-per-segment resampling keep the interpolant fixed, while deterministic
half-density, hull-preserving remeshing changes both interpolation and sampling.
Median and trimmed statistics are coordinatewise; neither is automatically a
geometric contour centre. Isolated retained crossings are part of the reference
mean but have no segment length and therefore do not contribute to the length
centroid. These support differences must be considered when present.

Fixed-iq comparisons retain all detected roots and form a factor only when
both groups have exactly one root. Ambiguous branches are left unresolved.
Open endpoints describe available mesh support, not proof of physical truncation.
The mirror diagnostic queries the full original interpolant at `(id,-iq)`:
its support fraction describes the convex hull, not complete edge-filtered
contour support. Missing mirrors remain missing. Symmetry is not imposed.

Initial observation at this one level/pair: the exact reference mean is reproduced;
edge-mean k is 13.0010684, length-centroid k is 13.3751163 A/A (+2.877%).
The q-coordinates change from -22.846/7.989 A to 1.376/0.393 A. Both centroids
remain off the selected flux level. These results demonstrate sensitivity in
this example, not a globally superior estimator or a physical causal partition.

## B/C: local scalar and two-axis factors

`investigate_local_factors(set_ids,min_rcond)` reuses saved 80-neighbour local
affine derivatives, with output index before input index. The scalar value is
the saved differential least-squares optimum

`k_loc = (a' b)/(a' a)`, `a=[Ldd;Lqd]`, `b=[Lde;Lqe]`.

It minimizes `norm(b-k*a)^2` locally with equal component weights. It is not a
global flux-error optimum. Its validity mask retains the reference neighbourhood
screen and `norm(a)>0.1*median(norm(a))`; no additional value clipping is applied.

The two-axis factors solve `J*[kd;kq]=b`, with
`J=[Ldd Ldq;Lqd Lqq]`. They require finite saved derivatives and independently
screen `rcond(J)>=1e-3`; this is a numerical choice, not a calibrated uncertainty
bound. Thresholds 1e-4, 1e-3 and 1e-2 are reported together. Ldq and Lqd are not
assumed equal. The stored global Method F constants are read for comparison,
without refitting. Printed summaries report variation and group medians.

For these saved derivatives, all 24,590 points pass the stated screens; the
minimum per-set rcond is approximately 0.170, 0.169 and 0.210. This does not
exclude sensitivity to derivative estimation or experimental uncertainty.
The 5th–95th percentiles of k_loc span approximately 4.27–15.43, 4.37–15.48
and 4.19–16.27 A/A. Local kd/kq variation is compared with global F; no equality
of the identification objectives is assumed.

## Existing work and open questions

| Topic | Status / next question |
|---|---|
| A: conventional contour factor | Existing reference; geometry sensitivity investigated here at one level/pair |
| B/C in the original paper | Existing fixed-iq d-flux and joint d/q-flux prediction criteria; no refit here |
| D/E/G | Existing excitation-/q-/all-current-dependent shifts; retain their results, no new optimization |
| F | Existing common two-axis constant shift; new local differential comparison reuses it |
| H | Existing three-current interpolation reference; no new global model |
| Determinant and eta | Existing nonparallel-sensitivity diagnostics; neither is a flux-error percentage |
| Local versus global k | Different weighting and objectives; a global LS factor is not the arithmetic mean of local factors |
| Small denominators | Separate sensitivity of excitation differences, scalar sensitivities and stator-Jacobian conditioning |
| Variable differential factors | A coefficient field does not justify multiplying it by ie to obtain a finite shift |
| Symmetry and cross-saturation | Approximate symmetry must be checked against both measured values and available support |
| Mesh/normalization/boundaries | Extend A over additional levels/pairs and alternative mesh choices before attributing a cause |

Two-axis `(id+kd*ie,iq+kq*ie)` is distinct from assigning different d-axis shifts
to psi_d and psi_q. With spatially varying kd/kq, local compensation directions
must be integrated and assessed for globally consistent coordinates; substituting
local coefficients into the constant-shift formula omits derivative terms.

Open hypotheses include intersection-density bias, incomplete-contour effects
and asymmetric support as contributors to centroid shifts. The first example
does not rank them or show that they explain all observed flux errors. Extend
the same diagnostics over the original levels/groups, examine perturbations and
fixed-branch comparisons, and compare common support before recommending a method.
Global flux approximation and the later `function_model_study` remain out of scope.
