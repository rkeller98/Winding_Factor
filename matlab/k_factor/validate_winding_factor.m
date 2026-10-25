function validate_winding_factor
%% Local derivative evidence, estimator diagnostics and analytic numerical controls.
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'common'));
root=setup_matlab_paths('k_factor');
out=fullfile(root,'results');
require_analysis_file(fullfile(out,'reduction_data.mat'),'prepare_reduction_data');
v=load(fullfile(out,'reduction_data.mat'),'data'); data=v.data;
neighbors=[40 80 120]; study=cell(3,1); summary=[];
for s=1:3
    d=data{s}; x=d.x; y=d.psi; [k,contour]=estimate_contour_factor(x,y,d.groups);
    slopes=cell(3,1);
    for h=1:3
        count=neighbors(h); sc=max(x)-min(x); z=x./sc; L=nan(size(x,1),2,3);
        for i=1:size(x,1)
            [distance,ix]=mink(sum((z-z(i,:)).^2,2),count);
            if numel(unique(d.groups(ix)))<2,continue;end
            weight=sqrt(exp(-distance/max(distance)));
            A=[ones(count,1),z(ix,:)-z(i,:)].*weight;
            [~,R]=qr(A,0); if rcond(R)<1e-6,continue;end
            b=A\(y(ix,:).*weight);
            L(i,:,:)=permute(b(2:4,:)./sc',[3 2 1]);
        end
        % Output index before input index: L(:,d,q) = d psi_d / d iq.
        a=L(:,:,1); b=L(:,:,3); an=vecnorm(a,2,2); bn=vecnorm(b,2,2);
        valid=all(isfinite([a b]),2) & an>.1*median(an,'omitnan');
        klocal=sum(a.*b,2)./sum(a.^2,2); klocal(~valid)=nan;
        eta=abs(a(:,1).*b(:,2)-a(:,2).*b(:,1))./(an.*bn);
        eta(~valid | bn<.1*median(bn,'omitnan'))=nan;
        kg=sum(a(valid,:).*b(valid,:),'all')/sum(a(valid,:).^2,'all');
        ratio=b./a;
        for c=1:2,ratio(abs(a(:,c))<.1*median(abs(a(:,c)),'omitnan'),c)=nan;end
        r0=b-k*a; rg=b-kg*a; rl=b-klocal.*a;
        slopes{h}=struct('L',L,'a',a,'b',b,'valid',valid,'eta',eta,'klocal',klocal, ...
            'kg',kg,'ratio',ratio,'residual',r0);
        summary=[summary;s,count,k,kg,median(klocal,'omitnan'),median(eta,'omitnan'), ...
            percentile(eta,.95),norm(r0(valid,:),'fro')/norm(b(valid,:),'fro'), ...
            norm(rg(valid,:),'fro')/norm(b(valid,:),'fro'), ...
            norm(rl(valid,:),'fro')/norm(b(valid,:),'fro'),nnz(valid)];
    end
    common=slopes{1}.valid & slopes{2}.valid & slopes{3}.valid;
    sensitivity=[];
    for h=1:3
        a=slopes{h}.a(common,:); b=slopes{h}.b(common,:); kg=sum(a.*b,'all')/sum(a.^2,'all');
        sensitivity=[sensitivity;neighbors(h),kg,norm(b-kg*a,'fro')/norm(b,'fro'),nnz(common)];
    end
    study{s}=struct('slopes',{slopes},'contour',contour,'k',k,'sensitivity',sensitivity);
    fprintf('Set%d: contour %.8f, derivative global %.8f, supported derivatives %d\n', ...
        s,k,slopes{2}.kg,nnz(slopes{2}.valid));
end
writetable(array2table(summary,'VariableNames',{'set','neighbors','legacy_k','global_derivative_k', ...
    'median_local_k','median_eta','p95_eta','relative_legacy_residual','relative_global_residual', ...
    'relative_local_residual','valid_points'}),fullfile(out,'derivative_metrics.csv'));

% Three small directly evaluated examples, not a test framework.
[D,Q,E]=ndgrid(linspace(-2,2,9),linspace(-1,1,7),linspace(0,3,6));
X=[D(:),Q(:),E(:)]; k=2; S=X(:,1)+k*X(:,3);
Y=[.3*S+.2*X(:,2),-.1*S+.4*X(:,2)];
% Move to a different physical point on the same reduced-coordinate fibre.
Xe=X; Xe(:,3)=Xe(:,3)+.137; Xe(:,1)=Xe(:,1)-k*.137;
R=[.3*(Xe(:,1)+k*Xe(:,3))+.2*Xe(:,2),-.1*(Xe(:,1)+k*Xe(:,3))+.4*Xe(:,2)];
constant_error=max(abs(Y-R),[],'all');
a=2*X(:,3)+.1*X(:,3).^2; kd=2+.2*X(:,3);
Y=[.3*(X(:,1)+a)+.2*X(:,2),-.1*(X(:,1)+a)+.4*X(:,2)];
Xe=X; Xe(:,3)=Xe(:,3)+.137; ae=2*Xe(:,3)+.1*Xe(:,3).^2;
Xe(:,1)=X(:,1)+a-ae;
correct=[.3*(Xe(:,1)+ae)+.2*Xe(:,2),-.1*(Xe(:,1)+ae)+.4*Xe(:,2)];
wrong=[.3*(X(:,1)+X(:,3).*kd)+.2*X(:,2),-.1*(X(:,1)+X(:,3).*kd)+.4*X(:,2)];
integrated_error=max(abs(Y-correct),[],'all'); product_error=max(abs(Y-wrong),[],'all');
% Same s=0,q=0 at (0,0,0) and (-2,0,1), incompatible q fluxes 0 and 2.
pair=[0 0 0;-2 0 1]; pair_y=[pair(:,1)+2*pair(:,3),.5*pair(:,1)+pair(:,2)+3*pair(:,3)];
incompatible_gap=pair_y(2,:)-pair_y(1,:); determinant=1*3-.5*2;
controls=table(constant_error,integrated_error,product_error,incompatible_gap(2),determinant, ...
    'VariableNames',{'constant_max_Vs','integrated_max_Vs','wrong_product_max_Vs','incompatible_gap_Vs','determinant'});
writetable(controls,fullfile(out,'analytic_controls.csv')); disp(controls);
save(fullfile(out,'winding_validation.mat'),'study','summary','neighbors','controls','-v7');

end

function p=percentile(x,q)
x=sort(x(isfinite(x))); p=x(ceil(q*numel(x)));
end
