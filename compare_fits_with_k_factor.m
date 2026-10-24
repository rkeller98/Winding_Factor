%% Same physical test points: differences from the trained full LUT, not measurement errors.
root=fileparts(mfilename('fullpath'));out=fullfile(root,'results','flux_fits');
a=load(fullfile(out,'data.mat'));data=a.data;
a=load(fullfile(out,'comparison.mat'));study=a.study;differences=[];
for s=1:numel(data)
 d=data{s};r=study{s};use=d.test&r.common;
 for m=1:size(r.prediction,3)
  e=1000*(r.prediction(use,:,m)-r.prediction(use,:,3));
  differences=[differences;s,m,sqrt(mean(e.^2,1)),sqrt(mean(e.^2,'all')), ...
   max(abs(e),[],'all'),size(e,1)];
 end
end
writetable(array2table(differences,'VariableNames',{'dataset','method','rmse_d_mVs', ...
 'rmse_q_mVs','joint_rmse_mVs','max_abs_mVs','points'}),fullfile(out,'differences_from_LUT.csv'));
% Existing all-data LUTs used the held-out groups: comparison only, not validation.
file=fullfile(root,'results','flux_map_reduction.mat');existing=[];
if ~isfile(file)
 unzip(fullfile(root,'winding_factor_study','research_results.zip'),fullfile(root,'results'));
end
old=load(file,'study');
for s=1:3
 d=data{s};r=study{s};p=reduction_lookup(old.study{s}.maps{8},d.x);
 use=d.test&r.common&all(isfinite(p),2);
 for m=1:size(r.prediction,3)
  e=1000*(r.prediction(use,:,m)-p(use,:));
  existing=[existing;s,m,sqrt(mean(e.^2,1)),sqrt(mean(e.^2,'all')), ...
   max(abs(e),[],'all'),size(e,1)];
 end
end
writetable(array2table(existing,'VariableNames',{'dataset','method','rmse_d_mVs', ...
 'rmse_q_mVs','joint_rmse_mVs','max_abs_mVs','points'}),fullfile(out,'existing_LUT_differences.csv'));
save(fullfile(out,'k_comparison.mat'),'differences','existing','-v7');
