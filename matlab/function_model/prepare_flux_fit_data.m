function prepare_flux_fit_data
% Prepare physical points and the existing training/validation/test split.

matlab_folder = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(matlab_folder, 'common'));
root = setup_matlab_paths('function_model');

output_folder = fullfile(root, 'results', 'flux_fits');
if ~isfolder(output_folder)
    mkdir(output_folder);
end

input_file = fullfile(root, 'results', 'reduction_data.mat');
require_analysis_file(input_file, 'run_matlab_analysis([],[],''prepare'')');
loaded = load(input_file, 'data');
measurements = loaded.data;

data = cell(size(measurements));
inventory = [];

for dataset_index = 1:numel(measurements)
    points = measurements{dataset_index};
    [currents, flux, groups, excitation] = reduction_inputs(points);
    group_count = numel(excitation);

    dataset = struct;
    dataset.name = sprintf('Set%d', dataset_index);
    dataset.x = currents;
    dataset.y = flux;
    dataset.group = groups;
    dataset.levels = excitation;
    dataset.active = {'id', 'iq', 'ie'};
    dataset.input_units = {'A', 'A', 'A'};
    dataset.group_name = 'Nominal excitation';
    dataset.group_unit = 'A';
    dataset.temperature = points.stator_temperature;
    dataset.speed = median(points.rpm_ref);
    dataset.source = points.meta.source_files;
    dataset.fixed_condition = 'Temperature was recorded but not controlled or corrected';

    % Reserve the penultimate group as a training anchor below the upper boundary.
    dataset.validation = ismember(groups, 2:4:(group_count - 2));
    dataset.test = ismember(groups, 4:4:(group_count - 2));
    dataset.boundary = groups == group_count;
    dataset.train = ~(dataset.validation | dataset.test | dataset.boundary);

    data{dataset_index} = dataset;
    inventory = [inventory; ...
        dataset_index, size(currents, 1), group_count, ...
        nnz(dataset.train), nnz(dataset.validation), nnz(dataset.test), ...
        nnz(dataset.boundary), min(dataset.temperature), max(dataset.temperature), ...
        dataset.speed];
end

% Keep the existing optional temperature example as a separate reference case.
fixture_file = fullfile(root, 'Resources', 'GUI_VICE-MeasurementEvalKit', ...
    'SourceCode', 'tests', 'Artifacts', 'json_api_headless_real_run', ...
    'psm_temperature', 'flux_psm_temperature_30C_70C.mat_result.mat');

if isfile(fixture_file)
    loaded = load(fixture_file, 'Total');
    fixture = loaded.Total.pos_2500_rpm.pos_0_A.Meas;

    currents = [fixture.PE_Id_meas_avg(:), fixture.PE_Iq_meas_avg(:), ...
                fixture.CAN_T_Rot_max(:)];
    omega = fixture.PE_w_el_rad_avg(:);
    stator_resistance = fixture.Const_Machine_R_s(:);
    flux = [
        (fixture.PE_Uq_meas_avg(:) - stator_resistance .* currents(:, 2)) ./ omega, ...
        (stator_resistance .* currents(:, 1) - fixture.PE_Ud_meas_avg(:)) ./ omega
    ];

    groups = 1 + (currents(:, 3) > 50);
    normalized_currents = (currents(:, 1:2) - min(currents(:, 1:2))) ...
        ./ (max(currents(:, 1:2)) - min(currents(:, 1:2)));
    spatial_tile = min(floor(4 * normalized_currents(:, 1)), 3) ...
        + 4 * min(floor(4 * normalized_currents(:, 2)), 3);
    split_role = mod(spatial_tile, 5);

    dataset = struct;
    dataset.name = 'Temperature fixture';
    dataset.x = currents;
    dataset.y = flux;
    dataset.group = groups;
    dataset.levels = [30; 70];
    dataset.active = {'id', 'iq', 'theta_rotor'};
    dataset.input_units = {'A', 'A', 'deg C'};
    dataset.group_name = 'Temperature band';
    dataset.group_unit = 'deg C';
    dataset.temperature = currents(:, 3);
    dataset.speed = 2500;
    dataset.source = ['Processed headless test artifact; absent acquisition fixture; ' ...
                      'unverified physical provenance'];
    dataset.fixed_condition = ['PSM-labelled artifact, nominal excitation 0 A; ' ...
                               'no independent temperature holdout'];

    dataset.validation = split_role == 1;
    dataset.test = split_role == 3;
    dataset.boundary = false(size(groups));
    dataset.train = ~(dataset.validation | dataset.test);

    data{end + 1} = dataset;
    inventory = [inventory; ...
        numel(data), size(currents, 1), 2, ...
        nnz(dataset.train), nnz(dataset.validation), nnz(dataset.test), 0, ...
        min(dataset.temperature), max(dataset.temperature), dataset.speed];
else
    fprintf(['Optional temperature fixture absent: %s\n' ...
        'Continuing with the selected measured-current datasets.\n'], fixture_file);
end

for dataset_index = 1:numel(data)
    dataset = data{dataset_index};
    training_rows = dataset.train;
    training_currents = dataset.x(training_rows, :);
    training_flux = dataset.y(training_rows, :);
    current_scale = max(abs(training_currents));

    map = fit_reduction_map( ...
        training_currents, training_flux, 8, [], current_scale, 81, true);
    lut_values = reduction_lookup(map, dataset.x);
    dataset.validation_support = all(isfinite(lut_values), 2) & dataset.validation;
    data{dataset_index} = dataset;

    fprintf(['%s: train %d, validation %d (%d LUT supported), ' ...
        'test %d, upper boundary %d\n'], ...
        dataset.name, nnz(training_rows), nnz(dataset.validation), ...
        nnz(dataset.validation_support), nnz(dataset.test), nnz(dataset.boundary));
end

inventory_columns = {
    'dataset', 'points', 'groups', 'training', 'validation', 'test', ...
    'upper_boundary', 'temperature_min_C', 'temperature_max_C', 'requested_rpm'
};
writetable(array2table(inventory, 'VariableNames', inventory_columns), ...
    fullfile(output_folder, 'data_inventory.csv'));
save(fullfile(output_folder, 'data.mat'), 'data', 'inventory', '-v7');
end
