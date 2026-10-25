function run_k_factor_analysis(data, export_format)
% Run the k-factor study with prepared data, or import the default datasets.

matlab_folder = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(matlab_folder, 'common'));
root = setup_matlab_paths('k_factor');
if nargin < 2
    export_format = '';
end

if nargin == 0 || isempty(data)
    addpath(root);
    run_matlab_analysis([], [], 'k_factor', export_format);
    return;
end

output_folder = fullfile(root, 'results');
if ~isfolder(output_folder)
    mkdir(output_folder);
end
save(fullfile(output_folder, 'reduction_data.mat'), 'data', '-v7');

steps = {
    @analyze_winding_factors
    @validate_winding_factor
    @analyze_flux_map_reduction
    @plot_reduction_study
};

descriptions = {
    'Analyze conventional winding factors'
    'Evaluate local derivatives and analytic controls'
    'Fit and compare coordinate reductions A-H'
    'Plot coordinate reduction study and write paper tables'
};

for step_index = 1:numel(steps)
    fprintf('\nk-factor [%d/%d]: %s\n', ...
        step_index, numel(steps), descriptions{step_index});
    if step_index == 1 || step_index == 4
        feval(steps{step_index}, export_format);
    else
        feval(steps{step_index});
    end
    if step_index == 2
        plot_local_winding_factors;
    end
end

drawnow;
end
