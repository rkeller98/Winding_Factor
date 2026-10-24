%% Common physical points and explicit training/validation/test partitions.
root=fileparts(mfilename('fullpath')); addpath(root);
out=fullfile(root,'results','flux_fits'); if ~isfolder(out),mkdir(out);end
data=cell(3,1); inventory=[];
for s=1:3
 a=load(fullfile(root,'Resources','DataSets',sprintf('Set%d',s),'extracted', ...
     'evaluation','research_flux_measurements.mat'),'Total');
 x=[];y=[];nominal=[];temperature=[];speed=[];
 fn=fieldnames(a.Total);
 for j=1:numel(fn)
  groups=fieldnames(a.Total.(fn{j}));
  for g=1:numel(groups)
   b=a.Total.(fn{j}).(groups{g}); m=b.Meas;
   x=[x; m.PE_Id_meas_avg(:),m.PE_Iq_meas_avg(:),m.I_exc_meas(:)];
   y=[y; m.PE_psi_d_meas(:),m.PE_psi_q_meas(:)];
   nominal=[nominal; b.Ref_exc_current+zeros(numel(m.I_exc_meas),1)];
   temperature=[temperature;m.CAN_T_Stat_max(:)];
   speed=[speed;b.Ref_Speed+zeros(numel(m.I_exc_meas),1)];
  end
 end
 [levels,~,group]=unique(nominal); ng=numel(levels);
 d=struct('name',sprintf('Set%d',s),'x',x,'y',y,'group',group,'levels',levels, ...
     'active',{{'id','iq','ie'}},'input_units',{{'A','A','A'}}, ...
     'group_name','Nominal excitation','group_unit','A','temperature',temperature, ...
     'speed',median(speed),'source','Fresh arithmetic channel means; stationary voltage-derived flux', ...
     'fixed_condition','Temperature was recorded but not controlled or corrected');
 % Reserve ng-1 as a training anchor below the separate upper boundary ng.
 d.validation=ismember(group,2:4:ng-2);d.test=ismember(group,4:4:ng-2);
 d.boundary=group==ng; d.train=~(d.validation|d.test|d.boundary);
 data{s}=d;
 inventory=[inventory;s,size(x,1),ng,nnz(d.train),nnz(d.validation),nnz(d.test), ...
     nnz(d.boundary),min(temperature),max(temperature),median(speed)];
end
% A processed test artifact exists, but its acquisition fixture is not available.
% It is a secondary reference exercise, never independent temperature evidence.
file=fullfile(root,'Resources','GUI_VICE-MeasurementEvalKit','SourceCode','tests', ...
 'Artifacts','json_api_headless_real_run','psm_temperature','flux_psm_temperature_30C_70C.mat_result.mat');
if isfile(file)
 a=load(file,'Total');b=a.Total.pos_2500_rpm.pos_0_A;m=b.Meas;
 x=[m.PE_Id_meas_avg(:),m.PE_Iq_meas_avg(:),m.CAN_T_Rot_max(:)];
 omega=m.PE_w_el_rad_avg(:);rs=m.Const_Machine_R_s(:);
 y=[(m.PE_Uq_meas_avg(:)-rs.*x(:,2))./omega, ...
     (rs.*x(:,1)-m.PE_Ud_meas_avg(:))./omega];
 group=1+(x(:,3)>50);z=(x(:,1:2)-min(x(:,1:2)))./(max(x(:,1:2))-min(x(:,1:2)));
 tile=min(floor(4*z(:,1)),3)+4*min(floor(4*z(:,2)),3);role=mod(tile,5);
 d=struct('name','Temperature fixture','x',x,'y',y,'group',group,'levels',[30;70], ...
  'active',{{'id','iq','theta_rotor'}},'input_units',{{'A','A','deg C'}}, ...
  'group_name','Temperature band','group_unit','deg C','temperature',x(:,3),'speed',2500, ...
  'source','Processed headless test artifact; absent acquisition fixture; unverified physical provenance', ...
  'fixed_condition','PSM-labelled artifact, nominal excitation 0 A; no independent temperature holdout');
 d.validation=role==1;d.test=role==3;d.boundary=false(size(group));
 d.train=~(d.validation|d.test); data{4}=d;
 inventory=[inventory;4,size(x,1),2,nnz(d.train),nnz(d.validation),nnz(d.test),0, ...
     min(x(:,3)),max(x(:,3)),2500];
end
for s=1:numel(data)
 d=data{s};tr=d.train;
 map=fit_reduction_map(d.x(tr,:),d.y(tr,:),8,[],max(abs(d.x(tr,:))),81,true);
 d.validation_support=all(isfinite(reduction_lookup(map,d.x)),2)&d.validation;
 data{s}=d;
 fprintf('%s: train %d, validation %d (%d LUT supported), test %d, upper boundary %d\n', ...
 d.name,nnz(tr),nnz(d.validation),nnz(d.validation_support),nnz(d.test),nnz(d.boundary));
end
writetable(array2table(inventory,'VariableNames',{'dataset','points','groups','training', ...
 'validation','test','upper_boundary','temperature_min_C','temperature_max_C','requested_rpm'}), ...
 fullfile(out,'data_inventory.csv'));
save(fullfile(out,'data.mat'),'data','inventory','-v7');
