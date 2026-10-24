%% Regenerate independent surface comparisons, every group figure and paper tables.
root=fileparts(mfilename('fullpath')); addpath(root); out=fullfile(root,'results');
paper=fullfile(root,'winding_factor_study'); figs=fullfile(paper,'figures'); atlas=fullfile(figs,'atlas');
if ~isfolder(atlas),mkdir(atlas);end
v=load(fullfile(out,'reduction_data.mat')); data=v.data;
v=load(fullfile(out,'flux_map_reduction.mat')); study=v.study; labels=v.labels;
v=load(fullfile(out,'winding_validation.mat')); diagnostics=v.study;
set(groot,'defaultFigureVisible','off','defaultAxesFontSize',14,'defaultTextFontSize',14, ...
    'defaultTextInterpreter','tex','defaultAxesTickLabelInterpreter','tex');
blue=[.05 .32 .72]; orange=[.94 .42 .08];
diverge=[linspace(.1,1,128)',linspace(.3,1,128)',ones(128,1); ...
    ones(128,1),linspace(1,.2,128)',linspace(1,.1,128)'];
atlaslist=fopen(fullfile(paper,'atlas_pages.tex'),'w');
pairrows=[]; contourrows=[]; foldrows=[]; jacrows=[]; boundaryrows=[];
for s=1:3
 d=data{s}; r=study{s}; x=d.x; y=d.psi; ng=numel(d.excitation);
 % Actual measured excitation is interpolated with each original surface.
 originals=cell(ng,1); xygrid=cell(ng,1); originalgrid=cell(ng,1); tri=cell(ng,1);
 for g=1:ng
    ix=d.groups==g; a=x(ix,:); b=y(ix,:);
    [D,Q]=meshgrid(linspace(min(a(:,1)),max(a(:,1)),41),linspace(min(a(:,2)),max(a(:,2)),41));
    original=surface_values(a(:,1:2),[b a(:,3)],[D(:),Q(:)],r.current_scale(1:2));
    xygrid{g}=cat(3,D,Q); originalgrid{g}=original;
    originals{g}=[D(:),Q(:),original(:,3)];
    tri{g}=physical_triangles(a(:,1:2),r.current_scale(1:2));
 end
 % Early view of original and transformed independent surfaces.
 fig=figure('Position',[20 20 1250 780]);theme(fig,'light'); t=tiledlayout(2,2,'TileSpacing','compact'); colors=turbo(ng);
 for c=1:2
  for transformed=0:1
   ax=nexttile; hold(ax,'on');
   for g=1:ng
    ix=d.groups==g; a=x(ix,:); b=y(ix,:); z=a(:,1:2);
    if transformed,z=reduction_coordinates(a,1,r.p{1},r.current_scale);end
    drawx=z(:,1);drawy=z(:,2);if c==1,drawx=z(:,2);drawy=z(:,1);end
    trisurf(tri{g},drawx,drawy,1000*b(:,c),'FaceColor','none', ...
        'EdgeColor',colors(g,:),'LineWidth',.45);
   end
   view(35,25); grid on;
   if c==1
    xlabel('i_q in A');ylabel('i_d in A');if transformed,ylabel('s in A');end
   else
    xlabel('i_d in A');ylabel('i_q in A');if transformed,xlabel('s in A');end
   end
   zlabel(sprintf('\\psi_%s in mVs',char('d'+13*(c-1))));
   title(sprintf('%s, %s flux',choose(transformed,'Original','Transformed'),char('d'+13*(c-1))));
  end
 end
 title(t,sprintf('Set%d: separate surfaces, A k = %.3f A/A; colour denotes excitation group',s,r.p{1}));
 colormap(fig,colors);clim([.5 ng+.5]);cb=colorbar;cb.Layout.Tile='east';
 cb.Ticks=1:ng;cb.TickLabels=string(d.excitation);cb.Label.String='Nominal excitation in A';
 exportgraphics(fig,fullfile(figs,sprintf('geometry_set%d.png',s)),'Resolution',150);close(fig);
 % One physical lower-slice query grid for ALL transformations in each pair.
 pairerror=cell(ng-1,1); pairmask=cell(ng-1,1);
 for g=1:ng-1
   query=originals{g}; orig=originalgrid{g}(:,1:2); error=nan(size(orig,1),2,7);
   hi=d.groups==g+1; reduced=cell(7,1);
   for m=1:7
     z=reduction_coordinates(x(hi,:),m,r.p{m},r.current_scale);
     qz=reduction_coordinates(query,m,r.p{m},r.current_scale); reduced{m}=qz;
     [pred,diagnosis]=mesh_values(z,y(hi,:),qz,tri{g+1},x(hi,1:2));
     low=d.groups==g; zlow=reduction_coordinates(x(low,:),m,r.p{m},r.current_scale);
     [checklow,lowdiagnosis]=mesh_values(zlow,y(low,:),qz,tri{g},x(low,1:2));
     pred(~all(isfinite(checklow),2),:)=nan;
     error(:,:,m)=1000*(pred-orig);
     if g==1
        foldrows=[foldrows;s,m,d.excitation(g),lowdiagnosis];
     end
     foldrows=[foldrows;s,m,d.excitation(g+1),diagnosis];
   end
   common=all(isfinite(error),[2 3]); pairerror{g}=error; pairmask{g}=common;
   for m=1:7
     e=error(common,:,m); n=nnz(common);
     if n>0,a=[sqrt(mean(e.^2,1)),sqrt(mean(e.^2,'all')),max(abs(e),[],'all')];else,a=nan(1,4);end
     pairrows=[pairrows;s,m,d.excitation(g),d.excitation(g+1),a,n,nnz(all(isfinite(orig),2))];
   end
 end
 writetable(array2table(pairrows,'VariableNames',{'set','method','ie_low_A','ie_high_A', ...
    'rmse_d_mVs','rmse_q_mVs','joint_rms_mVs','max_abs_mVs','points','eligible'}), ...
    fullfile(out,'surface_pairs.csv'));
 gp=ceil((ng-1)/2); lim=max(abs(pairerror{gp}),[],'all','omitnan'); lim=max(lim,1);
 fig=figure('Position',[20 20 1180 760]);theme(fig,'light'); t=tiledlayout(2,2,'TileSpacing','compact');
 for c=1:2
  for m=[1 7]
   nexttile; q=originals{gp}; z=reduction_coordinates(q,m,r.p{m},r.current_scale);
   e=pairerror{gp}(:,c,m); e(~pairmask{gp})=nan;
   D=xygrid{gp}(:,:,1); Q=xygrid{gp}(:,:,2);
   surf(reshape(z(:,1),size(D)),reshape(z(:,2),size(Q)),zeros(size(D)),reshape(e,size(D)), ...
       'EdgeColor','none');view(2);axis tight;colorbar;clim([-lim lim]);
   xlabel('s_d in A');ylabel('s_q in A');title(sprintf('%s: \\Delta\\psi_%s in mVs',labels{m}(1),char('d'+13*(c-1))));
  end
 end
 colormap(fig,diverge);title(t,sprintf('Set%d: independently shifted surfaces, %g minus %g A', ...
    s,d.excitation(gp+1),d.excitation(gp)));
 exportgraphics(fig,fullfile(figs,sprintf('surface_set%d.png',s)),'Resolution',150);close(fig);
 % Derivative maps and contour estimator evidence.
 a=diagnostics{s}; l=a.slopes{2}; g=ceil(ng/2); ix=d.groups==g;
 fig=figure('Position',[20 20 1240 660]);theme(fig,'light'); t=tiledlayout(2,3,'TileSpacing','compact');
 values={l.ratio(:,1),l.ratio(:,2),l.klocal,l.eta,1000*vecnorm(l.residual,2,2),1000*l.L(:,2,1)};
 titles={'L_{de}/L_{dd} in A/A','L_{qe}/L_{qd} in A/A','Local LS k in A/A', ...
    '\eta (dimensionless)','Derivative residual norm in mH','L_{qd} in mH'};
 for j=1:6,nexttile;scatter(x(ix,1),x(ix,2),14,values{j}(ix),'filled');title(titles{j});xlabel('i_d in A');ylabel('i_q in A');colorbar;end
 title(t,sprintf('Set%d, nominal excitation %g A; 80-neighbour local affine estimates',s,d.excitation(g)));
 exportgraphics(fig,fullfile(figs,sprintf('derivative_set%d.png',s)),'Resolution',150);close(fig);
 ce=a.contour.centroid_error; drift=a.contour.iq_drift;
 contourrows=[contourrows;s,a.k,1000*sqrt(mean(ce.^2,'all','omitnan')), ...
    median(abs(drift),'all','omitnan'),sum(isfinite(a.contour.ratios),'all')];
 % RMSE over excitation, same physical support across A--H.
 fig=figure('Position',[20 20 1240 780]);theme(fig,'light'); t=tiledlayout(2,2,'TileSpacing','compact');
 for mode=1:2
  for c=1:2
   nexttile; hold on;
   for m=1:8
    rows=r.pergroup(:,2)==mode & r.pergroup(:,3)==m;
    plot(r.pergroup(rows,4),r.pergroup(rows,4+c),'.-','LineWidth',1.1,'DisplayName',labels{m}(1));
   end
   grid on;xlabel('Nominal excitation in A');ylabel(sprintf('RMSE \\psi_%s in mVs',char('d'+13*(c-1))));
   title(choose(mode==2,'Reconstruction','Whole-group prediction'));legend('Location','best','NumColumns',4);
  end
 end
 title(t,sprintf('Set%d: identical supported physical points for all methods',s));
 exportgraphics(fig,fullfile(figs,sprintf('rmse_set%d.png',s)),'Resolution',150);close(fig);
 % Complete atlas: one page for each method and excitation group.
 lim=max(abs(1000*(r.reconstruction-y)),[],'all','omitnan'); lim=max(lim,1);
 for g=1:ng
  D=xygrid{g}(:,:,1); Q=xygrid{g}(:,:,2); orig=originalgrid{g}(:,1:2);
  query=originals{g}; supported=all(isfinite(query),2);
  for m=1:8
   pred=nan(size(orig));pred(supported,:)=reduction_lookup(r.maps{m},query(supported,:));
   e=1000*(pred-orig);
   fig=figure('Position',[20 20 1180 760]);theme(fig,'light'); t=tiledlayout(2,2,'TileSpacing','compact');
   for c=1:2
    nexttile;hold on;
    drawx=D;drawy=Q;if c==1,drawx=Q;drawy=D;end
    surf(drawx,drawy,reshape(1000*orig(:,c),size(D)),'FaceColor',blue,'FaceAlpha',.65,'EdgeColor','none');
    surf(drawx,drawy,reshape(1000*pred(:,c),size(D)),'FaceColor','none','EdgeColor',orange,'LineWidth',.45);
    view(35,25);grid on;
    if c==1,xlabel('i_q in A');ylabel('i_d in A');else,xlabel('i_d in A');ylabel('i_q in A');end
    zlabel(sprintf('\\psi_%s in mVs',char('d'+13*(c-1))));
    title('Original and reconstruction');legend('Original','Reconstruction','Location','best');
    nexttile;surf(D,Q,zeros(size(D)),reshape(e(:,c),size(D)),'EdgeColor','none');view(2);axis tight;
    xlabel('i_d in A');ylabel('i_q in A');title(sprintf('Signed r_%s in mVs',char('d'+13*(c-1))));colorbar;clim([-lim lim]);
   end
   colormap(fig,diverge);title(t,sprintf('Set%d | %s | nominal excitation %g A',s,labels{m},d.excitation(g)));
   name=sprintf('set%d_group%02d_%s.png',s,g,labels{m}(1));
   exportgraphics(fig,fullfile(atlas,name),'Resolution',105);close(fig);
   if g==ceil(ng/2) && m==3
    copyfile(fullfile(atlas,name),fullfile(figs,sprintf('reconstruction_set%d.png',s)));
   end
  end
 end
 fprintf(atlaslist,'\\AtlasDataset{%d}{',s);
 fprintf(atlaslist,'%g,',d.excitation(1:end-1));fprintf(atlaslist,'%g}\n',d.excitation(end));
 % Boundary groups: both extremes separately, with reduced-common and own support.
 for g=[1 ng]
  ix=d.groups==g; common=ix & r.boundary_common;
  for m=1:8
   own=ix & all(isfinite(r.prediction(:,:,m)),2);
   supported_common=common;
   if m==8,supported_common=common & all(isfinite(r.prediction(:,:,m)),2);end
   boundaryrows=[boundaryrows;s,m,d.excitation(g),briefstats(r.prediction(:,:,m)-y,supported_common,nnz(ix)), ...
       briefstats(r.prediction(:,:,m)-y,own,nnz(ix))];
  end
 end
 % G must also be checked where the final held-out physical points are queried.
 for g=1:ng
  tr=d.groups~=g;ix=~tr;sc=max(abs(x(tr,:)),[],1);p=r.fold_p{g,7};
  jac=1+x(ix,3).*(p(4)/sc(1)+2*p(5)*x(ix,1)/sc(1)^2);
  outside=any(x(ix,:)<min(x(tr,:),[],1) | x(ix,:)>max(x(tr,:),[],1),2);
  jacrows=[jacrows;s,d.excitation(g),min(jac),sum(jac<=0),sum(jac<.25),sum(outside),nnz(ix)];
 end
 % Correct H's boundary summary to count only finite predictions.
 inner=d.groups>1 & d.groups<ng;
 r.boundary(8,3:9)=briefstats(r.prediction(:,:,8)-y,r.boundary_common & all(isfinite(r.prediction(:,:,8)),2),nnz(~inner));
 study{s}=r;result=r;
 save(fullfile(out,sprintf('Set%d_reduction.mat',s)),'result','labels','-append');
 fprintf('Set%d figures complete\n',s);
end
fclose(atlaslist);
writetable(array2table(contourrows,'VariableNames',{'set','k','centroid_departure_rms_mVs', ...
    'median_abs_diq_die','finite_secants'}),fullfile(out,'contour_diagnostics.csv'));
writetable(array2table(foldrows,'VariableNames',{'set','method','ie_nominal_A','triangles', ...
    'reversed_triangles','min_area_ratio','physical_reversed_area_fraction','ambiguous_queries'}),fullfile(out,'projected_folding.csv'));
writetable(array2table(jacrows,'VariableNames',{'set','ie_nominal_A','min_G_jacobian','nonpositive', ...
    'below_training_margin','outside_training_box','points'}),fullfile(out,'holdout_jacobians.csv'));
writetable(array2table(boundaryrows,'VariableNames',{'set','method','ie_nominal_A', ...
    'common_rmse_d_mVs','common_rmse_q_mVs','common_joint_mVs','common_max_mVs','common_p95_mVs','common_n','common_coverage', ...
    'own_rmse_d_mVs','own_rmse_q_mVs','own_joint_mVs','own_max_mVs','own_p95_mVs','own_n','own_coverage'}),fullfile(out,'boundary_groups.csv'));
save(fullfile(out,'independent_surface_comparison.mat'),'pairrows','contourrows','foldrows','jacrows','boundaryrows');
save(fullfile(out,'flux_map_reduction.mat'),'study','-append');
% Overall errors and error--coverage curves, no arbitrary acceptance threshold.
fig=figure('Position',[20 20 1320 680]);theme(fig,'light');t=tiledlayout(2,3,'TileSpacing','compact');
for s=1:3
 nexttile; r=study{s}; a=r.metrics(r.metrics(:,2)==2,:);bar(a(:,3),a(:,6));
 xticks(1:8);xticklabels(cellstr(('A':'H')'));ylabel('Joint prediction RMSE in mVs');title(sprintf('Set%d',s));grid on;
end
for s=1:3
 nexttile;hold on;d=data{s};r=study{s};
 for m=1:8
  e=sort(1000*vecnorm(r.prediction(r.cv_common,:,m)-d.psi(r.cv_common,:),2,2));
  plot(e,(1:numel(e))/numel(e),'DisplayName',labels{m}(1));
 end
 xlabel('Vector error in mVs');ylabel('Fraction of common points');grid on;legend('Location','best','NumColumns',4);
end
exportgraphics(fig,fullfile(figs,'overview.png'),'Resolution',150);close(fig);
write_tables(paper,data,study,diagnostics,pairrows,contourrows);

function value=surface_values(x,y,q,sc)
validq=all(isfinite(q),2); value=nan(size(q,1),size(y,2));
z=x./sc; dt=delaunayTriangulation(z);t=dt.ConnectivityList;
edges=unique(sort([t(:,[1 2]);t(:,[1 3]);t(:,[2 3])],2),'rows');
lim=9*mean(sum((z(edges(:,1),:)-z(edges(:,2),:)).^2,2));
good=true(size(t,1),1);
for edge=[1 2;1 3;2 3]'
good=good & sum((z(t(:,edge(1)),:)-z(t(:,edge(2)),:)).^2,2)<lim;
end
[ix,b]=pointLocation(dt,q(validq,:)./sc);keep=isfinite(ix);keep(keep)=good(ix(keep));
rows=find(validq);rows=rows(keep);
for c=1:size(y,2),a=y(:,c);value(rows,c)=sum(b(keep,:).*a(t(ix(keep),:)),2);end
end
function tri=physical_triangles(x,sc)
x=x./sc;dt=delaunayTriangulation(x);tri=dt.ConnectivityList;edges=dt.edges;
lim=9*mean(sum((x(edges(:,1),:)-x(edges(:,2),:)).^2,2));good=true(size(tri,1),1);
for pair=[1 2;1 3;2 3]',good=good & sum((x(tri(:,pair(1)),:)-x(tri(:,pair(2)),:)).^2,2)<lim;end
tri=tri(good,:);
end
function [value,info]=mesh_values(z,v,q,t,physical)
% Preserve source mesh connectivity and reject multiple projected memberships.
u=z(t(:,2),:)-z(t(:,1),:);w=z(t(:,3),:)-z(t(:,1),:);
detz=u(:,1).*w(:,2)-u(:,2).*w(:,1);
a=physical(t(:,2),:)-physical(t(:,1),:);b=physical(t(:,3),:)-physical(t(:,1),:);
detp=a(:,1).*b(:,2)-a(:,2).*b(:,1);ratio=detz./detp;
count=zeros(size(q,1),1);positive=false(size(count));value=nan(size(q,1),size(v,2));
for j=1:size(t,1)
 if abs(detz(j))<1e-16,continue;end
 verts=z(t(j,:),:);lo=min(verts);hi=max(verts);
 ix=find(all(q>=lo-1e-9 & q<=hi+1e-9,2));if isempty(ix),continue;end
 d=q(ix,:)-verts(1,:);
 B=(d(:,1)*w(j,2)-d(:,2)*w(j,1))/detz(j);
 C=(u(j,1)*d(:,2)-u(j,2)*d(:,1))/detz(j);
 keep=B>=-1e-10 & C>=-1e-10 & B+C<=1+1e-10;ix=ix(keep);B=B(keep);C=C(keep);
 if isempty(ix),continue;end
 count(ix)=count(ix)+1;positive(ix)=ratio(j)>0;
 value(ix,:)=(1-B-C).*v(t(j,1),:)+B.*v(t(j,2),:)+C.*v(t(j,3),:);
end
value(count~=1 | ~positive,:)=nan;
info=[size(t,1),nnz(ratio<=0),min(ratio),sum(abs(detp(ratio<=0)))/sum(abs(detp)),nnz(count>1)];
end
function a=briefstats(r,mask,total)
r=1000*r(mask,:);n=size(r,1);if n==0,a=[nan nan nan nan nan 0 0];return;end
v=sort(vecnorm(r,2,2));a=[sqrt(mean(r.^2,1)),sqrt(mean(r.^2,'all')),max(abs(r),[],'all'),v(ceil(.95*n)),n,n/total];
end
function x=choose(condition,a,b)
if condition,x=b;else,x=a;end
end
function write_tables(paper,data,study,diagnostics,pairs,contours)
f=fopen(fullfile(paper,'results.tex'),'w');
fprintf(f,'%% Generated by plot_reduction_study.m. All errors in mVs.\n');
fprintf(f,'\\newcommand{\\ResultTables}{\n');
for s=1:3
 fprintf(f,'\\begin{table}[htbp]\\centering\\small\\caption{Set %d: errors on common support. R: reconstruction; P: interior whole-group prediction. Max: largest absolute component error; Q95: vector-error quantile.}\\begin{tabular}{llrrrrr}\\toprule Method & Mode & RMSE d & RMSE q & Joint & Max & Q95\\\\\\midrule\n',s);
 a=study{s}.metrics;
 for mode=1:2
  for m=1:8
   row=a(a(:,2)==mode & a(:,3)==m,:);typ='R';if mode==2,typ='P';end
   fprintf(f,'%s & %s & %.3f & %.3f & %.3f & %.2f & %.2f\\\\\n',char('A'+m-1),typ,row(4:8));
  end
 end
 fprintf(f,'\\bottomrule\\end{tabular}\\par\\small Common points: reconstruction %d/%d; interior prediction %d/%d.\\end{table}\n', ...
  nnz(study{s}.recon_common),size(data{s}.x,1),nnz(study{s}.cv_common),nnz(data{s}.groups>1 & data{s}.groups<numel(data{s}.excitation)));
end
fprintf(f,'}\n\\newcommand{\\FactorTable}{\\begin{tabular}{lrrrr}\\toprule Set & A & B & C & Derivative LS\\\\\\midrule\n');
for s=1:3,fprintf(f,'%d & %.4f & %.4f & %.4f & %.4f\\\\\n',s,study{s}.p{1:3},diagnostics{s}.slopes{2}.kg);end
fprintf(f,'\\bottomrule\\end{tabular}}\n');
fprintf(f,'\\newcommand{\\ParameterTable}{\\begin{tabular}{llrrrrr}\\toprule Set & Method & $c_1$ & $c_2$ & $c_3$ & $c_4$ & $c_5$\\\\\\midrule\n');
for s=1:3
 for m=4:7
  p=nan(1,5);p(1:numel(study{s}.p{m}))=study{s}.p{m};
  fprintf(f,'%d & %s & %.4f & %.4f & %.4f & %.4f & %.4f\\\\\n',s,char('A'+m-1),p);
 end
end
fprintf(f,'\\bottomrule\\end{tabular}}\n');
fprintf(f,'\\newcommand{\\DerivativeTable}{\\begin{tabular}{lrrrrr}\\toprule Set & $k_{\\rm LS}$ & median $\\eta$ & $\\rho_A$ & $\\rho_{LS}$ & $\\rho_{local}$\\\\\\midrule\n');
for s=1:3
 l=diagnostics{s}.slopes{2};valid=l.valid;a=l.a(valid,:);b=l.b(valid,:);
 fprintf(f,'%d & %.3f & %.3f & %.3f & %.3f & %.3f\\\\\n',s,l.kg,median(l.eta,'omitnan'), ...
 norm(b-diagnostics{s}.k*a,'fro')/norm(b,'fro'),norm(b-l.kg*a,'fro')/norm(b,'fro'),norm(b-l.klocal(valid).*a,'fro')/norm(b,'fro'));
end
fprintf(f,'\\bottomrule\\end{tabular}}\n');
fprintf(f,'\\newcommand{\\SurfaceTable}{\\begin{tabular}{lrrrrrrr}\\toprule Set & A & B & C & D & E & F & G\\\\\\midrule\n');
for s=1:3
 fprintf(f,'%d',s);
 for m=1:7
  rows=pairs(:,1)==s & pairs(:,2)==m; a=pairs(rows,:);
  valid=a(:,9)>0 & isfinite(a(:,7));
  val=sqrt(sum(a(valid,7).^2.*a(valid,9))/sum(a(valid,9)));fprintf(f,' & %.3f',val);
 end
 fprintf(f,'\\\\\n');
 fprintf(f,'\\multicolumn{8}{l}{Common pair anchors: %d/%d.}\\\\\n',sum(a(:,9)),sum(a(:,10)));
end
fprintf(f,'\\bottomrule\\end{tabular}}\n');
fprintf(f,'\\newcommand{\\NumericalTable}{\\begin{tabular}{llrrr}\\toprule Set & Control & RMSE d & RMSE q & Points\\\\\\midrule\n');
for s=1:3
 r=study{s};
 for j=1:4
  label={'A normalized','A raw axes','H normalized','H raw axes'};
  fprintf(f,'%d & %s & %.3f & %.3f & %d\\\\\n',s,label{j},r.scaling_stats(j,[2 3 7]));
 end
 for j=1:3
  fprintf(f,'%d & A raster %d & %.3f & %.3f & %d\\\\\n',s,r.raster_stats(j,[1 2 3 7]));
 end
end
fprintf(f,'\\bottomrule\\end{tabular}}\n');
fprintf(f,'\\newcommand{\\ContourTable}{\\begin{tabular}{lrrr}\\toprule Set & $k_A$ & Centroid RMS in mVs & median $|\\Delta i_q/\\Delta i_e|$\\\\\\midrule\n');
for s=1:3,fprintf(f,'%d & %.4f & %.3f & %.3f\\\\\n',contours(s,1:4));end
fprintf(f,'\\bottomrule\\end{tabular}}\n');
fprintf(f,'\\newcommand{\\SensitivityTable}{\\begin{tabular}{lrrr}\\toprule Set & Neighbours & Derivative LS & $\\rho$\\\\\\midrule\n');
for s=1:3
 a=diagnostics{s}.sensitivity;
 for j=1:3,fprintf(f,'%d & %d & %.4f & %.4f\\\\\n',s,a(j,1:3));end
end
fprintf(f,'\\bottomrule\\end{tabular}}\n');
fprintf(f,'\\newcommand{\\NominalTable}{\\begin{tabular}{lrrr}\\toprule Set & RMS d in mVs & RMS q in mVs & Points\\\\\\midrule\n');
for s=1:3,fprintf(f,'%d & %.3f & %.3f & %d\\\\\n',s,study{s}.reference_stats([1 2 6]));end
fprintf(f,'\\bottomrule\\end{tabular}}\n');
fprintf(f,'\\newcommand{\\OptimizerTable}{\\begin{tabular}{lrrrrrrr}\\toprule Set & B & C & D & E & F & G & Runs/model\\\\\\midrule\n');
for s=1:3
 fprintf(f,'%d',s);r=study{s};ng=numel(data{s}.excitation);
 for m=2:7
  failed=r.fitinfo{m}.exitflag<=0;
  for g=1:ng,failed=failed+(r.fold_info{g}{m}.exitflag<=0);end
  fprintf(f,' & %d',failed);
 end
 fprintf(f,' & %d\\\\\n',ng+1);
end
fprintf(f,'\\bottomrule\\end{tabular}}\n');
fprintf(f,'\\newcommand{\\BoundaryTable}{\\begin{tabular}{lrlrrr}\\toprule Set & Excitation & Method & Joint in mVs & Points & Coverage\\\\\\midrule\n');
for s=1:3
 d=data{s};r=study{s};ng=numel(d.excitation);
 for g=[1 ng]
  for m=[3 7 8]
   use=d.groups==g & r.boundary_common;
   if m==8,use=d.groups==g & all(isfinite(r.prediction(:,:,m)),2);end
   a=briefstats(r.prediction(:,:,m)-d.psi,use,nnz(d.groups==g));
   fprintf(f,'%d & %g & %s & %.3f & %d & %.3f\\\\\n',s,d.excitation(g),char('A'+m-1),a([3 6 7]));
  end
 end
end
fprintf(f,'\\bottomrule\\end{tabular}}\n');
fclose(f);
end
