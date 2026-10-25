function data = prepare_reduction_data(files, input_type)
% Import raw measurements or one Eval MAT into the common measurement struct.
% Currents are in A, mechanical speeds in rpm, and flux linkages in Vs.

if ischar(files) || isstring(files)
    files = cellstr(files);
end

% Common field name             Source channel
channel_map = {
    'id_meas',                  'PE_Id_meas_avg'
    'iq_meas',                  'PE_Iq_meas_avg'
    'ie_meas',                  'I_exc_meas'
    'rpm_meas',                 'CAN_RPM_meas'
    'id_req',                   'PE_Id_req_meas'
    'iq_req',                   'PE_Iq_req_meas'
    'ie_req',                   'I_exc_req'
    'rpm_req',                  'PE_Rpm_req'
    'psi_d',                    'PE_psi_d_meas'
    'psi_q',                    'PE_psi_q_meas'
    'ud_meas',                  'PE_Ud_meas_avg'
    'uq_meas',                  'PE_Uq_meas_avg'
    'omega_el',                 'PE_w_el_rad_avg'
    'stator_resistance',        'Const_Machine_R_s'
    'stator_temperature',       'CAN_T_Stat_max'
    'rotor_temperature',        'CAN_T_Rot_max'
    'machine_id',               'Const_MACHINE_ID_param_none'
    'phase_selector',           'MultiPhase_control_request_sel'
    'current_sensor_delay',     'Const_Current_Sensor_delay_param_sec'
    'angle_sensor_delay',       'Const_Angle_Sensor_delay_param_sec'
    'i_ac_rms',                 'PE_I_AC_meas_Arms'
    'point_id',                 'PE_Index_meas'
};

data = struct;
for channel_index = 1:size(channel_map, 1)
    field_name = channel_map{channel_index, 1};
    data.(field_name) = zeros(0, 1);
end
data.ie_ref = zeros(0, 1);
data.rpm_ref = zeros(0, 1);
data.source_id = zeros(0, 1);
data.repeats = zeros(0, 1);

switch lower(input_type)
    case 'raw'
        data = prepare_raw_files(files, channel_map, data);

    case 'eval'
        if numel(files) ~= 1
            error('Select one Eval file.');
        end
        data = load_eval_file(files{1}, channel_map, data);

    otherwise
        error('Use input_type ''raw'' or ''eval''.');
end

% Required physical quantities determine which points can be analyzed.
% Missing optional channels do not remove an otherwise valid point.
required_values = [data.id_meas, data.iq_meas, data.ie_meas, ...
                   data.psi_d, data.psi_q];
valid_points = all(isfinite(required_values), 2);
point_rows = find(valid_points);

% Sort every channel together; measured excitation remains unchanged.
references = [data.rpm_ref(valid_points), data.ie_ref(valid_points)];
[~, sort_order] = sortrows(references, [1 2]);
point_rows = point_rows(sort_order);

point_fields = fieldnames(data);
for field_index = 1:numel(point_fields)
    field_name = point_fields{field_index};
    data.(field_name) = data.(field_name)(point_rows);
end

data.meta.source_files = files(:);
end


function data = prepare_raw_files(files, channel_map, data)
% Average repetitions within each source file, then calculate stationary flux.

point_fields = fieldnames(data);
channel_count = size(channel_map, 1);

for file_index = 1:numel(files)
    loaded = load(files{file_index});
    variable_names = fieldnames(loaded);
    raw = loaded.(variable_names{1});
    source_channels = erase({raw.Y.Name}, ' ');

    index_channel = strcmp(source_channels, 'PE_Index_meas');
    measurement_index = raw.Y(index_channel).Data(:);
    sample_count = numel(measurement_index);

    samples = struct;
    for channel_index = 1:channel_count
        field_name = channel_map{channel_index, 1};
        source_name = channel_map{channel_index, 2};
        source_index = find(strcmp(source_channels, source_name), 1);

        channel_values = nan(sample_count, 1);
        if ~isempty(source_index)
            channel_values(:) = raw.Y(source_index).Data(:);
        end
        samples.(field_name) = channel_values;
    end

    flux_inputs = [samples.id_meas, samples.iq_meas, samples.ie_meas, ...
                   samples.ud_meas, samples.uq_meas, samples.omega_el, ...
                   samples.stator_resistance];
    valid_samples = measurement_index > 0 & all(isfinite(flux_inputs), 2);

    % The same index can occur at different speed/excitation requests.
    grouping_keys = [samples.rpm_req(valid_samples), ...
                     samples.ie_req(valid_samples), ...
                     measurement_index(valid_samples)];
    [group_keys, ~, group_index] = unique(grouping_keys, 'rows');

    averaged = struct;
    for channel_index = 1:channel_count
        field_name = channel_map{channel_index, 1};
        channel_values = samples.(field_name);
        averaged.(field_name) = accumarray( ...
            group_index, channel_values(valid_samples), [], @mean);
    end

    averaged.psi_d = (averaged.uq_meas ...
        - averaged.stator_resistance .* averaged.iq_meas) ./ averaged.omega_el;
    averaged.psi_q = (averaged.stator_resistance .* averaged.id_meas ...
        - averaged.ud_meas) ./ averaged.omega_el;

    low_speed = abs(averaged.omega_el) <= 1;
    averaged.psi_d(low_speed) = NaN;
    averaged.psi_q(low_speed) = NaN;

    point_count = size(group_keys, 1);
    averaged.ie_ref = group_keys(:, 2);
    averaged.rpm_ref = group_keys(:, 1);
    averaged.source_id = repmat(file_index, point_count, 1);
    averaged.repeats = accumarray(group_index, 1);

    % Concatenate only after averaging, so files cannot share repetitions.
    for field_index = 1:numel(point_fields)
        field_name = point_fields{field_index};
        data.(field_name) = [data.(field_name); averaged.(field_name)];
    end
end
end


function data = load_eval_file(file, channel_map, data)
% Copy saved Meas values and block references without averaging or refitting.

loaded = load(file, 'Total');
speed_fields = fieldnames(loaded.Total);
speed_fields = speed_fields(endsWith(speed_fields, '_rpm'));

for speed_index = 1:numel(speed_fields)
    speed_block = loaded.Total.(speed_fields{speed_index});
    excitation_fields = fieldnames(speed_block);
    excitation_fields = excitation_fields(endsWith(excitation_fields, '_A'));

    for excitation_index = 1:numel(excitation_fields)
        block = speed_block.(excitation_fields{excitation_index});
        measurements = block.Meas;
        point_count = numel(measurements.PE_Id_meas_avg);

        for channel_index = 1:size(channel_map, 1)
            field_name = channel_map{channel_index, 1};
            source_name = channel_map{channel_index, 2};

            channel_values = nan(point_count, 1);
            if isfield(measurements, source_name)
                channel_values(:) = measurements.(source_name)(:);
            end
            data.(field_name) = [data.(field_name); channel_values];
        end

        data.ie_ref = [data.ie_ref; repmat(block.Ref_exc_current, point_count, 1)];
        data.rpm_ref = [data.rpm_ref; repmat(block.Ref_Speed, point_count, 1)];
        data.source_id = [data.source_id; ones(point_count, 1)];
        data.repeats = [data.repeats; nan(point_count, 1)];
    end
end
end
