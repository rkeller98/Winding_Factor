"""Prepare common physical data and training/validation/test partitions.

MATLAB equivalent: matlab/function_model/prepare_flux_fit_data.m.
Planned inputs: Set1-3 research_flux_measurements; the temperature artifact if available.
Planned outputs: Datasets with conditions, groups, split masks, and a LUT-supported validation mask.
Required preparation: prepare_reduction_data; respect the optional temperature source and its limited provenance.
Status: Not implemented yet; to be developed by the user.

Future convention: Control standalone execution through
if __name__ == "__main__":. Importing must not run analyses or open figures.
Call calculations and plotting explicitly. Individual analysis scripts may
show their own figures. The shared entry point will later control execution
order and joint figure display through python/function_model/run_function_model_analysis.py.
This file currently contains only a docstring.
"""
