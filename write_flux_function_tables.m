function write_flux_function_tables
root=fileparts(mfilename('fullpath'));paper=fullfile(root,'function_model_study');
a=load(fullfile(root,'results','flux_fits','data.mat'));data=a.data;
a=load(fullfile(root,'results','flux_fits','comparison.mat'));study=a.study;labels=a.labels;
out=fullfile(fileparts(paper),'results','flux_fits');M=readtable(fullfile(out,'metrics.csv'));
S=readtable(fullfile(out,'model_sizes.csv'));C=readtable(fullfile(out,'runtime.csv'));
L=readtable(fullfile(out,'derivatives.csv'));D=readtable(fullfile(out,'differences_from_LUT.csv'));
f=fopen(fullfile(paper,'results.tex'),'w');
fprintf(f,'\\newcommand{\\PredictionTable}{\\begin{longtable}{llrrrrr}\\toprule Set & Model & d & q & Joint & Points & Coverage\\\\\\midrule\n');
for s=1:numel(data)
 for m=1:numel(study{s}.models)
  a=M(M.dataset==s&M.mode==2&M.group==0&M.method==m,:);
  fprintf(f,'%d & %s & %.3f & %.3f & %.3f & %d & %.1f\\%%\\\\\n', ...
   s,labels{m},a.rmse_d_mVs,a.rmse_q_mVs,a.joint_rmse_mVs,a.valid_points,100*a.coverage);
 end
end
fprintf(f,'\\bottomrule\\end{longtable}}\n');
fprintf(f,'\\newcommand{\\CostTable}{\\begin{longtable}{llrrrr}\\toprule Set & Model & Terms & Coefficients & KiB & PC median in $\\mu$s\\\\\\midrule\n');
for s=1:numel(data)
 for m=1:numel(study{s}.models)
  a=S(S.dataset==s&S.method==m,:);b=C(C.dataset==s&C.method==m,:);
  fprintf(f,'%d & %s & %d & %d & %.2f & %.1f\\\\\n',s,labels{m},a.terms,a.flux_coefficients,a.online_numeric_bytes/1024,b.median_us);
 end
end
fprintf(f,'\\bottomrule\\end{longtable}}\n');
fprintf(f,'\\newcommand{\\DerivativeTable}{\\begin{longtable}{llrrrrrr}\\toprule Set & Model & $L_{dd}$ & $L_{dq}$ & $L_{de}$ & $L_{qd}$ & $L_{qq}$ & $L_{qe}$\\\\\\midrule\n');
for s=1:3
 for m=[1 2 3 4 9]
  fprintf(f,'%d & %s',s,labels{m});
  for c=1:2
   for r=1:3
    a=L(L.dataset==s&L.method==m&L.output==c&L.input==r,:);fprintf(f,' & %.3f',a.rmse_scaled);
   end
  end
  fprintf(f,'\\\\\n');
 end
end
fprintf(f,'\\bottomrule\\end{longtable}}\n');
fprintf(f,'\\newcommand{\\DifferenceTable}{\\begin{tabular}{lrrrr}\\toprule Set & RBF & Spline & kA & kF\\\\\\midrule\n');
for s=1:3
 fprintf(f,'%d',s);
 for m=[1 2 4 9],a=D(D.dataset==s&D.method==m,:);fprintf(f,' & %.3f',a.joint_rmse_mVs);end
 fprintf(f,'\\\\\n');
end
fprintf(f,'\\bottomrule\\end{tabular}}\n');
fprintf(f,'\\newcommand{\\ReconstructionTable}{\\begin{tabular}{lrrrrr}\\toprule Set & RBF & Spline & LUT & kA & kF\\\\\\midrule\n');
for s=1:3
 fprintf(f,'%d',s);
 for m=[1 2 3 4 9]
  a=M(M.dataset==s&M.mode==1&M.group==0&M.method==m,:);fprintf(f,' & %.3f',a.joint_rmse_mVs);
 end
 fprintf(f,'\\\\\n');
end
fprintf(f,'\\bottomrule\\end{tabular}}\n');
fprintf(f,'\\newcommand{\\ResultSummary}{On the matched excitation tests, joint RMSE for Sets 1/2/3 is ');
for m=1:3
 v=zeros(1,3);
 for s=1:3,a=M(M.dataset==s&M.mode==2&M.group==0&M.method==m,:);v(s)=a.joint_rmse_mVs;end
 fprintf(f,'%.3f/%.3f/%.3f mVs for %s',v,labels{m});
 if m<3,fprintf(f,'; ');else,fprintf(f,'. ');end
end
fprintf(f,'The spline has the smallest joint flux error among these three full-input candidates in all three current datasets; the Gaussian model has the shortest supplied MATLAB evaluation time. ');
fprintf(f,'The spline payload is %.2f/%.2f/%.2f KiB, compared with %.2f/%.2f/%.2f KiB for the full LUT. ', ...
 S.online_numeric_bytes(S.dataset<=3&S.method==2)'/1024,S.online_numeric_bytes(S.dataset<=3&S.method==3)'/1024);
fprintf(f,'This is an accuracy/storage/runtime tradeoff, not an established application pass criterion. ');
if numel(data)>3
 fprintf(f,'The processed artifact has a different ranking and only 44 supported test points out of 142; it cannot establish physical or temperature generalization. ');
end
fprintf(f,'Upper-boundary groups have zero accepted predictions because extrapolation is disabled, rather than an extrapolation accuracy result.}\n');
fclose(f);
end
