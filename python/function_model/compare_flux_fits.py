"""Compare saved fits for flux, derivatives, cost, and common support.

MATLAB equivalent: matlab/function_model/compare_flux_fits.m; eval_flux_model.
Planned inputs: Function-model data and saved RBF, B-spline, and baseline models.
Planned outputs: Predictions, Jacobians, masks, and error, size, and runtime tables.
Required preparation: prepare_flux_fit_data, fit_flux_rbf, fit_flux_bspline, and fit_flux_baselines; comparisons do not refit models.
Status: Not implemented yet; to be developed by the user.

Future convention: Control standalone execution through
if __name__ == "__main__":. Importing must not run analyses or open figures.
Call calculations and plotting explicitly. Individual analysis scripts may
show their own figures. The shared entry point will later control execution
order and joint figure display through python/function_model/run_function_model_analysis.py.
This file currently contains only a docstring.
"""
