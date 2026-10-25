function run_function_model_analysis
%% Rebuild the independent function-model study. No production changes.
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'common'));
root=setup_matlab_paths('function_model');
% Shared raw-data preparation is needed only when a compact input is absent.
for s=1:3
    file=fullfile(root,'Resources','DataSets',sprintf('Set%d',s), ...
        'extracted','evaluation','research_flux_measurements.mat');
    if ~isfile(file)
        fprintf('\nPreparing missing grouped measurement inputs: %s\n',file);
        prepare_reduction_data;
        break;
    end
end
fprintf('\nFunction models [1/7]: prepare physical points and data split\n');
prepare_flux_fit_data;
fprintf('\nFunction models [2/7]: fit Gaussian RBF models\n');
fit_flux_rbf;
fprintf('\nFunction models [3/7]: fit cubic B-spline models\n');
fit_flux_bspline;
fprintf('\nFunction models [4/7]: fit full LUT and k models A-G\n');
fit_flux_baselines;
fprintf('\nFunction models [5/7]: compare fluxes, derivatives and evaluation costs\n');
compare_flux_fits;
fprintf('\nFunction models [6/7]: compare with fresh and existing LUT references\n');
compare_fits_with_k_factor;
fprintf('\nFunction models [7/7]: plot all groups and write paper tables\n');
plot_flux_function_study; % Also calls write_flux_function_tables.
drawnow;
end
