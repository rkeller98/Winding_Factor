"""Compare fresh and historical preprocessing and saved coordinates.

MATLAB equivalent: matlab/k_factor/audit_reduction_sources.m.
Planned inputs: Reduction data and complete historical Set1/2 FluxEval files.
Planned outputs: Tables of flux differences, coordinate conventions, and machine metadata.
Required preparation: prepare_reduction_data; provide historical Total, LUT, and Calc contents separately.
Status: Not implemented yet; to be developed by the user.

Future convention: Control standalone execution through
if __name__ == "__main__":. Importing must not run analyses or open figures.
Call calculations and plotting explicitly. Individual analysis scripts may
show their own figures. The shared entry point will later control execution
order and joint figure display through python/k_factor/run_k_factor_analysis.py.
This file currently contains only a docstring.
"""
