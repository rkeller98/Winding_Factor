%% Conventional estimator: independently reproduce its secants in all datasets.
root=fileparts(mfilename('fullpath'));addpath(root);out=fullfile(root,'results');
figs=fullfile(root,'winding_factor_study/figures');if ~isfolder(figs),mkdir(figs);end
v=load(fullfile(out,'reduction_data.mat'));data=v.data;factors=cell(3,1);
set(groot,'defaultFigureVisible','off','defaultAxesFontSize',14,'defaultTextFontSize',14);
for s=1:3
 d=data{s};[k,c]=estimate_contour_factor(d.x,d.psi,d.groups);
 factors{s}=struct('k',k,'contour',c);
 fig=figure('Position',[20 20 1250 730]);theme(fig,'light');t=tiledlayout(2,2,'TileSpacing','compact');
 nexttile;plot(1000*c.levels,c.ratios','LineWidth',.6);grid on;
 xlabel('Selected d-flux level in mVs');ylabel('Absolute-numerator secant in A/A');title(sprintf('Conventional factor: %.4f A/A',k));
 nexttile;plot(1000*c.levels,c.signed','LineWidth',.6);grid on;
 xlabel('Selected d-flux level in mVs');ylabel('Signed -\Delta i_d/\Delta i_e in A/A');title('Sign retained for diagnosis');
 nexttile;plot(1000*c.levels,c.iq_drift','LineWidth',.6);grid on;
 xlabel('Selected d-flux level in mVs');ylabel('\Delta i_q/\Delta i_e in A/A');title('Contour-centroid q-current drift');
 nexttile;plot(1000*c.levels,1000*c.centroid_error','LineWidth',.6);grid on;
 xlabel('Selected d-flux level in mVs');ylabel('Centroid flux minus level in mVs');title('A centroid need not lie on its contour');
 title(t,sprintf('Set%d: contour levels and adjacent excitation groups',s));
 exportgraphics(fig,fullfile(figs,sprintf('contours_set%d.png',s)),'Resolution',150);close(fig);
 fprintf('Set%d conventional factor: %.12f A/A\n',s,k);
end
save(fullfile(out,'winding_factors.mat'),'factors');
