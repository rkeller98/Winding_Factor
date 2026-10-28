function [eval_set] = load_data(path)

arguments (Input)
    path (1,1) string = "C:\Git\Winding_Factor\Resources\DataSets\Set2\FluxEval\eval_flux_eesm_excitation_0A_to_40A.mat"
end

arguments (Output)
    eval_set (1,1) struct
end

data = load(path);

assert(isfield(data, "Total"), "Dataset is not a FluxEval")
assert(isfield(data, "Inv_Param"), "Missing Inv_Param")

eval_set.ctrl_maps = data.Inv_Param.CtrlMaps;
eval_set.ctrl_maps_4d = data.Inv_Param.Maps4D.CtrlMaps;

channel_map = {
    'id_meas', 'PE_Id_meas_avg';
    'iq_meas', 'PE_Iq_meas_avg';
    'ie_meas', 'I_exc_meas';
    'rpm_meas', 'CAN_RPM_meas'
    'id_req', 'PE_Id_req_meas';
    'iq_req', 'PE_Iq_req_meas';
    'ie_req', 'I_exc_req';
    'rpm_req', 'PE_Rpm_req';
    'psi_d', 'PE_psi_d_meas';
    'psi_q', 'PE_psi_q_meas';
    'ud_meas', 'PE_Ud_meas_avg';
    'uq_meas', 'PE_Uq_meas_avg';
    'omega_el', 'PE_w_el_rad_avg';
    'stator_resistance', 'Const_Machine_R_s';
    'stator_temperature', 'CAN_T_Stat_max';
    'rotor_temperature', 'CAN_T_Rot_max';
    'machine_id', 'Const_MACHINE_ID_param_none';
    'phase_selector', 'MultiPhase_control_request_sel';
    'current_sensor_delay', 'Const_Current_Sensor_delay_param_sec';
    'angle_sensor_delay', 'Const_Angle_Sensor_delay_param_sec';
    'i_ac_rms', 'PE_I_AC_meas_Arms';
    'point_id', 'PE_Index_meas';
    };

measurement = cell2struct( ...
    cell(size(channel_map, 1), 1), ...
    channel_map(:, 1), ...
    1);

f_rpm_field_mask = @(x) ~cellfun(@isempty, regexp(x, "(pos|neg)_(0|[1-9][0-9]*)_rpm"));
f_iexc_field_mask = @(x) ~cellfun(@isempty, regexp(x, "(pos|neg)_(0|[1-9][0-9]*)_A"));

fields = fieldnames(data.Total);
rpm_field_mask = f_rpm_field_mask(fields);
rpm_fields = fields(rpm_field_mask);

for idx_rpm = 1:numel(rpm_fields)
    rpm_data = data.Total.(rpm_fields{idx_rpm});
    fields = fieldnames(rpm_data);
    iexc_field_mask = f_iexc_field_mask(fields);
    iexc_fields = fields(iexc_field_mask);

    iexc_values = cellfun(@(field) rpm_data.(field).Ref_exc_current, iexc_fields);
    [~, sorted_idx] = sort(iexc_values);
    iexc_fields = iexc_fields(sorted_idx);

    for idx_iexc = 1:numel(iexc_fields)
        iexc_data = rpm_data.(iexc_fields{idx_iexc});
        meas = iexc_data.Meas;

        N = numel(meas.(channel_map{1,2}));

        for idx_channel = 1:size(channel_map,1)
            target = channel_map{idx_channel, 1};
            source = channel_map{idx_channel, 2};

            if isfield(meas, source)
                values = meas.(source);
                values = values(:);
                if isscalar(values)
                    values = repmat(values, N, 1);
                end
            else
                values = nan(N,1);
            end

            measurement.(target) = [measurement.(target); values];
        end

    end
end

eval_set.measurement = measurement;
end