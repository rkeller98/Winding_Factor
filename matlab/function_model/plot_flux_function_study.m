function plot_flux_function_study
%% Same surface/error layout as the existing study, in a separate output folder.
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'common'));
root=setup_matlab_paths('function_model');
out=fullfile(root,'results','flux_fits');
paper=fullfile(root,'function_model_study');figs=fullfile(paper,'figures');
if ~isfolder(figs),mkdir(figs);end
atlas=fullfile(figs,'atlas');if ~isfolder(atlas),mkdir(atlas);end
require_analysis_file(fullfile(out,'data.mat'),'prepare_flux_fit_data');
require_analysis_file(fullfile(out,'comparison.mat'),'compare_flux_fits');
require_analysis_file(fullfile(out,'rbf_candidates.csv'),'fit_flux_rbf');
require_analysis_file(fullfile(out,'bspline_candidates.csv'),'fit_flux_bspline');
require_analysis_file(fullfile(out,'differences_from_LUT.csv'),'compare_fits_with_k_factor');
a=load(fullfile(out,'data.mat'));data=a.data;
a=load(fullfile(out,'comparison.mat'));study=a.study;labels=a.labels;metrics=a.metrics;
graphics_properties={'defaultAxesFontSize','defaultTextFontSize'};
graphics_defaults=get(groot,graphics_properties);
restore_graphics=onCleanup(@() set(groot,graphics_properties,graphics_defaults));
set(groot,'defaultAxesFontSize',12,'defaultTextFontSize',12);
blue=[.15 .4 .75];orange=[.95 .4 .1];atlaslist=fopen(fullfile(paper,'atlas_pages.tex'),'w');
for s=1:numel(data)
 d=data{s};r=study{s};nm=numel(r.models);ng=numel(d.levels);res=1000*(r.prediction-d.y);
 keep=repmat(r.common,1,2,nm);limit=max(5,5*ceil(max(abs(res(keep)))/5));
 for g=1:ng
  ix=d.group==g;x=d.x(ix,:);y=d.y(ix,:);
  [D,Q]=ndgrid(linspace(min(x(:,1)),max(x(:,1)),31),linspace(min(x(:,2)),max(x(:,2)),31));
  original=surface_values(x(:,1:2),[y,x(:,3)],[D(:),Q(:)]);
  query=[D(:),Q(:),original(:,3)];reference=original(:,1:2);P=nan(size(query,1),2,nm);
  for m=1:nm,P(:,:,m)=eval_flux_model(r.models{m},query);end
  common=all(isfinite(reshape(P,size(query,1),[])),2)&all(isfinite(reference),2);
  if s==4,role='spatial tiles; no held-out temperature band';
  elseif any(d.boundary(ix)),role='upper boundary: outside training range';
  elseif any(d.test(ix)),role='held-out excitation group';
  elseif any(d.validation(ix)),role='validation group, included in final refit';
  else,role='training group';end
  for m=1:nm
   fig=figure('Visible','on','Position',[20 20 1360 740]);theme(fig,'light');
   layout=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
   title(layout,sprintf('%s | %s | %s %g %s | %s',d.name,labels{m}, ...
    d.group_name,d.levels(g),d.group_unit,role),'FontSize',12);
   for c=1:2
    nexttile;hold on;
    if c==1,X=Q;Y=D;else,X=D;Y=Q;end
    surf(X,Y,reshape(1000*reference(:,c),size(D)),'FaceColor',blue,'FaceAlpha',.7,'EdgeColor','none');
    v=P(:,c,m);v(~common)=nan;
    surf(X,Y,reshape(1000*v,size(D)),'FaceColor','none','EdgeColor',orange,'LineWidth',.5);
    view(42,27);grid on;
    if c==1,xlabel('i_q in A');ylabel('i_d in A');else,xlabel('i_d in A');ylabel('i_q in A');end
    zlabel(sprintf('\\psi_%s in mVs',char('d'+13*(c-1))));title('Reference and reconstruction');
    legend('Reference','Reconstruction','Location','northwest','FontSize',9);
    nexttile;e=1000*(P(:,c,m)-reference(:,c));e(~common)=nan;
    surface(D,Q,zeros(size(D)),reshape(e,size(D)),'EdgeColor','none');view(2);axis tight;
    xlabel('i_d in A');ylabel('i_q in A');title(sprintf('Signed r_%s in mVs',char('d'+13*(c-1))));
    colorbar;clim([-limit limit]);
   end
   colormap(fig,diverging_map);
   file=fullfile(atlas,sprintf('set%d_group%02d_%s.png',s,g,labels{m}));
   exportgraphics(fig,file,'Resolution',100);
   if s==2&&g==4&&m<=2,copyfile(file,fullfile(figs,['example_' lower(labels{m}) '.png']));end
  end
  fprintf('%s group %d/%d figures complete\n',d.name,g,ng);
 end
 methods=strjoin(labels(1:nm),',');levels=strjoin(cellstr(string(d.levels')),',');
 fprintf(atlaslist,'\\AtlasDataset{%d}{%s}{%s}\n',s,levels,methods);
end
fclose(atlaslist);
% Requested summary: error over group, and model size versus validation error.
fig=figure('Visible','on','Position',[20 20 1400 820]);theme(fig,'light');tiledlayout(2,2,'TileSpacing','compact');
for s=1:numel(data)
 nexttile;hold on;d=data{s};r=study{s};
 for m=1:numel(r.models)
  rows=metrics(:,1)==s&metrics(:,2)==2&metrics(:,3)>0&metrics(:,4)==m;
  b=metrics(rows,:);valid=b(:,12)>0;
  plot(d.levels(b(valid,3)),b(valid,11),'-o','DisplayName',labels{m});
 end
 xlabel([d.group_name ' in ' d.group_unit]);ylabel('Joint test RMSE in mVs');title(d.name);grid on;
 legend('Location','best','NumColumns',3,'FontSize',9);
end
exportgraphics(fig,fullfile(figs,'group_errors.png'),'Resolution',130);
b=readtable(fullfile(out,'rbf_candidates.csv'));c=readtable(fullfile(out,'bspline_candidates.csv'));
fig=figure('Visible','on','Position',[20 20 1320 780]);theme(fig,'light');tiledlayout(2,2,'TileSpacing','compact');
for s=1:numel(data)
 nexttile;hold on;bb=b(b.dataset==s,:);cc=c(c.dataset==s,:);
 scatter(bb.numeric_bytes/1024,sqrt(bb.validation_mVs2),38,'o','DisplayName','RBF candidates');
 scatter(cc.numeric_bytes/1024,sqrt(cc.validation_mVs2),38,'s','DisplayName','Spline candidates');
 xlabel('Numeric online storage in KiB');ylabel('Validation joint RMSE in mVs');
 title(data{s}.name);grid on;legend('Location','best');
end
exportgraphics(fig,fullfile(figs,'size_error.png'),'Resolution',130);
write_flux_function_tables;

end

function v=surface_values(x,y,q)
sc=max(x)-min(x);z=x./sc;dt=delaunayTriangulation(z);t=dt.ConnectivityList;e=edges(dt);
lim=9*mean(sum((z(e(:,1),:)-z(e(:,2),:)).^2,2));good=true(size(t,1),1);
for pair=[1 2;1 3;2 3]',good=good&sum((z(t(:,pair(1)),:)-z(t(:,pair(2)),:)).^2,2)<lim;end
[ix,b]=pointLocation(dt,q./sc);keep=isfinite(ix);keep(keep)=good(ix(keep));v=nan(size(q,1),size(y,2));
for c=1:size(y,2),a=y(:,c);v(keep,c)=sum(b(keep,:).*a(t(ix(keep),:)),2);end
end
function map=diverging_map
a=linspace(0,1,128)';map=[[.15+.85*a,.35+.65*a,ones(128,1)]; ...
 [ones(128,1),1-.8*a,1-.85*a]];
end
