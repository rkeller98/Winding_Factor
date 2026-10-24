%% Compare historical exports without replacing the new common preprocessing.
root=fileparts(mfilename('fullpath')); out=fullfile(root,'results');
v=load(fullfile(out,'reduction_data.mat'));data=v.data; rows=[]; convention=[];
files={'Matlab_200_201_202_FluxEval.mat','eval_flux_eesm_excitation_0A_to_40A.mat'};
for s=1:2
 source=fullfile(root,'Resources/DataSets',sprintf('Set%d',s),'FluxEval',files{s});
 old=load(source,'Total','LUT'); blocks=old.Total.pos_3000_rpm;names=fieldnames(blocks);
 names=names(contains(names,'_A')); d=data{s};
 for j=1:numel(names)
  block=blocks.(names{j}); m=block.Meas;nominal=block.Ref_exc_current;
  x=[m.PE_Id_meas_avg(:),m.PE_Iq_meas_avg(:),m.I_exc_meas(:)];
  psi=[m.PE_psi_d_meas(:),m.PE_psi_q_meas(:)];
  % Compare at the SAME physical coordinates, not by record index.
  use=d.excitation(d.groups)==nominal;
  F=scatteredInterpolant(d.x(use,1:2),d.psi(use,:),'linear','none');
  e=1000*(F(x(:,1:2))-psi);keep=all(isfinite(e),2);
  rows=[rows;s,nominal,numel(m.I_exc_meas),nnz(use),nnz(keep),sqrt(mean(e(keep,:).^2,1)),max(abs(e(keep,:)),[],'all')];
  if isfield(old,'LUT') && isfield(old.LUT.CtrlMaps,'complete')
   c=old.LUT.CtrlMaps.complete; slice=find(abs(c.I_exc_vec-nominal)<1e-6,1);
   if ~isempty(slice)
    axis=c.Id_breakpoints(slice,:);phys=linspace(min(block.Calc.Grid.Id_fine),max(block.Calc.Grid.Id_fine),numel(axis));
    k=c.EESM_winding_ratio;offset=k*nominal;
    [best,idx]=min([max(abs(axis-phys)),max(abs(axis-phys-offset))]);
    metadata="absent";
    if isfield(c,'Coordinates'),metadata=string(c.Coordinates.d_axis);end
    convention=[convention;{s,nominal,char(metadata),idx,best,offset}];
   end
  end
 end
 fprintf('Historical provenance comparison Set%d complete\n',s);
end
writetable(array2table(rows,'VariableNames',{'set','ie_nominal_A','historical_points','fresh_points', ...
 'matched_support','rmse_d_difference_mVs','rmse_q_difference_mVs','max_difference_mVs'}),fullfile(out,'historical_preprocessing_comparison.csv'));
writetable(cell2table(convention,'VariableNames',{'set','ie_nominal_A','metadata_d_axis', ...
 'inferred_convention_1physical_2reference','axis_alignment_error_A','fixed_reference_offset_A'}),fullfile(out,'historical_coordinates.csv'));
% Metadata and phase-amplitude check from freshly averaged acquisition channels.
meta=[];
for s=1:3
 a=data{s}.values;
 meta=[meta;s,median(a(:,17:20),1),median(vecnorm(a(:,1:2),2,2)./(sqrt(2)*a(:,21)),'omitnan'), ...
 min(a(:,11)),max(a(:,11)),min(a(:,15:16)),max(a(:,15:16))];
end
writetable(array2table(meta,'VariableNames',{'set','machine_id','phase_selector', ...
 'current_sensor_delay_s','angle_sensor_delay_s','median_peak_to_rms_ratio', ...
 'rpm_min','rpm_max','stator_temp_min_C','rotor_temp_min_C','stator_temp_max_C','rotor_temp_max_C'}),fullfile(out,'machine_metadata.csv'));
