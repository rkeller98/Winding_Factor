"""Fit the full LUT reference and existing k models A-G on the same training rows.

MATLAB equivalent: matlab/function_model/fit_flux_baselines.m; fit_flux_k_models and fit_reduction_map.
Planned inputs: Function-model data; final training plus validation rows.
Planned outputs: Saved LUT/k models, fit information, and training identities.
Required preparation: prepare_flux_fit_data; preserve existing model families and budgets. No Delaunay port in this task.
Status: Not implemented yet; to be developed by the user.

Future convention: Control standalone execution through
if __name__ == "__main__":. Importing must not run analyses or open figures.
Call calculations and plotting explicitly. Individual analysis scripts may
show their own figures. The shared entry point will later control execution
order and joint figure display through python/function_model/run_function_model_analysis.py.
This file currently contains only a docstring.
"""
