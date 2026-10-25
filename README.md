# EESM stator-flux coordinate reduction

The separate coefficient-model comparison is in
[function_model_study](function_model_study/RUN.md). Start it with
`run_flux_function_study`: shared Gaussian RBFs, cubic tensor splines,
full LUT and the existing k variants, including Set3. The original reduction
paper below is preserved.

Research comparison of seven k-based two-coordinate models and a complete
three-current reference, using Set1, Set2 and the freshly processed Set3.
The paper and figure labels remain in English.

The existing MeasEval Delaunay and edge-intersection implementations are included
in [Resources/MeasEval](Resources/MeasEval/README.md), including `PointSet`,
`SampledField` and the MATLAB geometry helpers from `modernization-baseline`.

From this repository in MATLAB R2025b:

```matlab
prepare_reduction_data
analyze_winding_factors
validate_winding_factor
audit_reduction_sources
analyze_flux_map_reduction
plot_reduction_study
```

Only base MATLAB is required for the new research calculations. The three independent
sets can also run in separate MATLAB processes with `analyze_flux_map_reduction(1)`,
respectively `(2)` and `(3)`; the last completed set collects the combined results.

- Acquisition files: `Resources/DataSets/Set*/raw/`, unchanged.
- Raw and extracted measurements are also available as compressed JSON (`.json.gz`); see the
  [Python loading example](Resources/DataSets/README.md). The MAT files are retained.
- Fresh grouped flux exports: each set's `extracted/evaluation/research_flux_measurements.mat`.
- Numerical output, masks, fitted parameters, source audit and per-group CSVs: `results/`.
- English paper: [main.pdf](winding_factor_study/main.pdf).
- Complete A--H group atlas: [atlas.pdf](winding_factor_study/atlas.pdf).
- Execution, data provenance and review mapping: [RUN.md](winding_factor_study/RUN.md).

Compile `main.tex` and `atlas.tex` twice with `pdflatex` from the
`winding_factor_study` directory after generating the results and figures.

The analysis uses measured currents and identical supported physical points for each
comparison. It distinguishes derivative diagnostics, independent transformed surfaces,
shared-table reconstruction and prediction of entirely omitted excitation groups.
Saved historical GUI tables are a separate provenance control. No productive GUI,
controller software, raw data or magnetic scalar model is changed.
