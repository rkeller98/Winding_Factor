function fit_flux_baselines
%% Fit the original LUT and k methods using the same final training rows.
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'common'));
root=setup_matlab_paths('function_model');
out=fullfile(root,'results','flux_fits');
require_analysis_file(fullfile(out,'data.mat'),'prepare_flux_fit_data');
a=load(fullfile(out,'data.mat'));data=a.data;models=cell(numel(data),8);information=cell(numel(data),1);
for s=1:numel(data)
 d=data{s};use=d.train|d.validation;x=d.x(use,:);y=d.y(use,:);sc=max(abs(x));
 models{s,1}=wrap(fit_reduction_map(x,y,8,[],sc,81,true),d,use);
 if s<=3
  [p,info,sc]=fit_flux_k_models(x,y,d.group(use),81);information{s}=info;
  for m=1:7,models{s,m+1}=wrap(fit_reduction_map(x,y,m,p{m},sc,81,true),d,use);end
 end
 fprintf('%s baselines fitted\n',d.name);
 save(fullfile(out,'baseline_models.mat'),'models','information','-v7');
end

end

function m=wrap(map,d,use)
x=d.x(use,:);
m=struct('kind','lut','map',map,'bounds',[min(x);max(x)],'active',{d.active}, ...
 'input_units',{d.input_units},'conditions',d.fixed_condition,'training_rows',find(use), ...
 'validity','Training physical box and finite interpolation cells; no extrapolation');
end
