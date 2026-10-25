"""Compare coordinate models A-H for reconstruction and group holdout.

MATLAB equivalent: matlab/k_factor/analyze_flux_map_reduction.m; fit_reduction_map and reduction_lookup.
Planned inputs: Prepared current/flux points, groups, and unchanged MATLAB settings.
Planned outputs: Reduction maps, fit parameters, predictions, common masks, and error tables.
Required preparation: prepare_reduction_data; preserve nested group holdout. Delaunay interpolation is not ported here.
Status: Not implemented yet; to be developed by the user.

Future convention: Control standalone execution through
if __name__ == "__main__":. Importing must not run analyses or open figures.
Call calculations and plotting explicitly. Individual analysis scripts may
show their own figures. The shared entry point will later control execution
order and joint figure display through python/k_factor/run_k_factor_analysis.py.
This file currently contains only a docstring.
"""
