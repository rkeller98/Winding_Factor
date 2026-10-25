function run_function_model_analysis(data, export_format)
% Run the function-model study with prepared data, or import default datasets.

matlab_folder = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(matlab_folder, 'common'));
root = setup_matlab_paths('function_model');
if nargin < 2
    export_format = '';
end

if nargin == 0 || isempty(data)
    addpath(root);
    run_matlab_analysis([], [], 'function_model', export_format);
    return;
end

output_folder = fullfile(root, 'results');
if ~isfolder(output_folder)
    mkdir(output_folder);
end
save(fullfile(output_folder, 'reduction_data.mat'), 'data', '-v7');

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

fprintf('\nFunction models [6/7]: compare with this study''s LUT reference\n');
compare_fits_with_k_factor(false);

fprintf('\nFunction models [7/7]: plot all groups and write paper tables\n');
plot_flux_function_study(export_format);

drawnow;
end
