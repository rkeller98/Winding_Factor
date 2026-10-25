function prepare_reduction_data
%% Rebuild all three research datasets from acquisition channels.
% Same arithmetic averaging and stationary flux formula for every dataset.
% No new outlier decisions, fitted surfaces or productive GUI changes.
root=setup_matlab_paths;
out = fullfile(root,'results'); if ~isfolder(out), mkdir(out); end
rawfiles = {{'Matlab_200.mat','Matlab_201.mat','Matlab_202.mat'}, ...
    {'flux_eesm_excitation_0A.mat','flux_eesm_excitation_2A_to_40A.mat'}, ...
    {'flux_eesm_excitation_3A_to_30A.mat'}};
channels = {'PE_Id_meas_avg','PE_Iq_meas_avg','I_exc_meas', ...
    'PE_Ud_meas_avg','PE_Uq_meas_avg','PE_w_el_rad_avg', ...
    'PE_Id_req_meas','PE_Iq_req_meas','I_exc_req','PE_Rpm_req', ...
    'CAN_RPM_meas','Const_Machine_R_s','Const_Machine_R_exc', ...
    'Const_Machine_pole_pairs','CAN_T_Stat_max','CAN_T_Rot_max', ...
    'Const_MACHINE_ID_param_none','MultiPhase_control_request_sel', ...
    'Const_Current_Sensor_delay_param_sec','Const_Angle_Sensor_delay_param_sec', ...
    'PE_I_AC_meas_Arms','CAN_U_AC_meas_Vrms'};
data = cell(3,1); inventory = [];
for s=1:3
    values=[]; repeats=[]; file_id=[]; point_id=[]; sample_std=[]; omitted=0;
    for f=1:numel(rawfiles{s})
        file=fullfile(root,'Resources','DataSets',sprintf('Set%d',s),'raw',rawfiles{s}{f});
        require_analysis_file(file,'Restore the original acquisition MAT file under Resources/DataSets/Set*/raw; compact exports cannot replace acquisition channels.');
        loaded=load(file); names=fieldnames(loaded); raw=loaded.(names{1});
        names=erase({raw.Y.Name},' ');
        ind=raw.Y(strcmp(names,'PE_Index_meas')).Data(:);
        a=nan(numel(ind),numel(channels));
        for j=1:numel(channels)
            k=find(strcmp(names,channels{j}),1);
            if ~isempty(k), a(:,j)=raw.Y(k).Data(:); end
        end
        keep=ind>0 & all(isfinite(a(:,1:14)),2);
        omitted=omitted+nnz(~keep);
        [key,~,group]=unique([a(keep,10),a(keep,9),ind(keep)],'rows');
        averaged=nan(size(key,1),size(a,2)); sd=averaged;
        for j=1:size(a,2)
            averaged(:,j)=accumarray(group,a(keep,j),[],@mean);
            sd(:,j)=accumarray(group,a(keep,j),[],@std);
        end
        values=[values;averaged]; repeats=[repeats;accumarray(group,1)];
        sample_std=[sample_std;sd]; file_id=[file_id;f+zeros(size(key,1),1)];
        point_id=[point_id;key(:,3)];
    end
    excitation=unique(values(:,9)); [~,groups]=ismember(values(:,9),excitation);
    x=values(:,1:3); rs=values(:,12); omega=values(:,6);
    psi=[(values(:,5)-rs.*x(:,2))./omega,(rs.*x(:,1)-values(:,4))./omega];
    keep=all(isfinite([x psi]),2) & abs(omega)>1;
    v=struct('x',x(keep,:),'psi',psi(keep,:),'groups',groups(keep), ...
        'excitation',excitation,'values',values(keep,:),'channels',{channels}, ...
        'sample_std',sample_std(keep,:),'repeats',repeats(keep), ...
        'file_id',file_id(keep),'point_id',point_id(keep),'rawfiles',{rawfiles{s}}, ...
        'omitted_samples',omitted,'software',version,'preprocessing', ...
        'Arithmetic channel means per file/speed/excitation/index; index 0 excluded; no outlier rejection');
    data{s}=v;
    % Compact grouped research export, not a saved productive GUI/LUT state.
    Total=struct;
    for g=1:numel(excitation)
        use=v.groups==g; speed=median(v.values(use,10));
        block.Ref_Speed=speed;block.Ref_exc_current=excitation(g);
        block.Meas=struct;
        for j=1:numel(channels),block.Meas.(channels{j})=v.values(use,j)';end
        block.Meas.PE_psi_d_meas=v.psi(use,1)';block.Meas.PE_psi_q_meas=v.psi(use,2)';
        Total.(sprintf('pos_%g_rpm',speed)).(sprintf('pos_%g_A',excitation(g)))=block;
    end
    folder=fullfile(root,'Resources','DataSets',sprintf('Set%d',s),'extracted','evaluation');
    if ~isfolder(folder),mkdir(folder);end
    save(fullfile(folder,'research_flux_measurements.mat'),'Total','-v7');
    for g=1:numel(excitation)
        use=v.groups==g; a=v.values(use,:);
        inventory=[inventory;s,excitation(g),nnz(use),sum(v.repeats(use)), ...
            min(a(:,1:3)),max(a(:,1:3)),mean(a(:,3)),std(a(:,3)), ...
            sqrt(mean((a(:,1:3)-a(:,7:9)).^2)), ...
            min(a(:,11)),max(a(:,11)),median(a(:,12:14)), ...
            min(a(:,15:16)),max(a(:,15:16))];
    end
    fprintf('Set%d: %d operating points, %d groups, %d excluded acquisition records\n', ...
        s,size(v.x,1),numel(excitation),omitted);
    % Standalone fresh flux export, including previously unavailable Set3.
    writetable(array2table([v.groups,v.x,v.psi,v.values(:,7:10),v.repeats,v.file_id,v.point_id], ...
        'VariableNames',{'group','id_A','iq_A','ie_A','psi_d_Vs','psi_q_Vs', ...
        'id_requested_A','iq_requested_A','ie_nominal_A','rpm_requested','repetitions','file','point_index'}), ...
        fullfile(out,sprintf('Set%d_flux_points.csv',s)));
end
audit=array2table(inventory,'VariableNames',{'set','ie_nominal_A','points','samples', ...
    'id_min_A','iq_min_A','ie_min_A','id_max_A','iq_max_A','ie_max_A','ie_mean_A','ie_std_A', ...
    'id_tracking_rms_A','iq_tracking_rms_A','ie_tracking_rms_A','rpm_min','rpm_max', ...
    'Rs_ohm','Re_ohm','pole_pairs','stator_temp_min_C','rotor_temp_min_C', ...
    'stator_temp_max_C','rotor_temp_max_C'});
writetable(audit,fullfile(out,'data_inventory.csv'));
save(fullfile(out,'reduction_data.mat'),'data','audit','channels','rawfiles','-v7');

end
