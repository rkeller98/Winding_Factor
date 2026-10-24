%% Rebuild the independent function-model study. No production changes.
root=fileparts(mfilename('fullpath'));addpath(root);
prepare_flux_fit_data;
fit_flux_rbf;
fit_flux_bspline;
fit_flux_baselines;
compare_flux_fits;
compare_fits_with_k_factor;
plot_flux_function_study;
