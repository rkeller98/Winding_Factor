function data = run_matlab_analysis(source_files, input_type, study, export_format)
% Select files below, or pass a raw path/list or one Eval path as arguments.
% Optional fourth argument: image format, e.g. 'png' or 'pdf'; default: no export.

root = fileparts(mfilename('fullpath'));

if nargin < 1 || isempty(source_files)
    dataset_folder = fullfile(root, 'Resources', 'DataSets');
    source_files = {
        fullfile(dataset_folder, 'Set1', 'raw', ...
            {'Matlab_200.mat', 'Matlab_201.mat', 'Matlab_202.mat'})
        fullfile(dataset_folder, 'Set2', 'raw', ...
            {'flux_eesm_excitation_0A.mat', 'flux_eesm_excitation_2A_to_40A.mat'})
        {fullfile(dataset_folder, 'Set3', 'raw', 'flux_eesm_excitation_3A_to_30A.mat')}
    };
end

if nargin < 2 || isempty(input_type)
    input_type = 'raw';
end

if nargin < 3
    study = 'both';
end
if nargin < 4
    export_format = '';
end

addpath( ...
    fullfile(root, 'matlab', 'common'), ...
    fullfile(root, 'matlab', 'k_factor'), ...
    fullfile(root, 'matlab', 'function_model'));

% A flat path list is one dataset; nested lists describe separate datasets.
if ischar(source_files) || isstring(source_files)
    source_files = cellstr(source_files);
end

if ~iscell(source_files{1})
    source_files = {source_files};
end

data = {};
for dataset_index = 1:numel(source_files)
    measurements = prepare_reduction_data(source_files{dataset_index}, input_type);

    missing_references = ~isfinite(measurements.ie_ref) ...
        | ~isfinite(measurements.rpm_ref);
    if any(missing_references)
        error('The studies require known ie_ref and rpm_ref for grouping.');
    end

    % Each speed is an independent physical map for the existing algorithms.
    speeds = unique(measurements.rpm_ref);
    point_fields = fieldnames(measurements);
    point_fields(strcmp(point_fields, 'meta')) = [];

    for speed_index = 1:numel(speeds)
        speed_rows = measurements.rpm_ref == speeds(speed_index);
        speed_data = measurements;

        for field_index = 1:numel(point_fields)
            field_name = point_fields{field_index};
            speed_data.(field_name) = measurements.(field_name)(speed_rows);
        end

        data{end + 1, 1} = speed_data;
    end
end

output_folder = fullfile(root, 'results');
if ~isfolder(output_folder)
    mkdir(output_folder);
end
save(fullfile(output_folder, 'reduction_data.mat'), 'data', '-v7');

if strcmp(study, 'prepare')
    return;
end

if any(strcmp(study, {'both', 'k_factor'}))
    fprintf('\n[1/2] k-factor and coordinate-reduction study\n');
    run_k_factor_analysis(data, export_format);
end

if any(strcmp(study, {'both', 'function_model'}))
    fprintf('\n[2/2] function-model study\n');
    run_function_model_analysis(data, export_format);
end

drawnow;
fprintf('\nAnalysis finished. Results: %s\nAll analysis figures remain open.\n', ...
    output_folder);
end
