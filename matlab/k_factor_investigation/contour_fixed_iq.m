function fixed = contour_fixed_iq(groups, iq_queries)
% Unique contour roots at shared physical iq; actual interpolated ie is retained.
if nargin < 2
    iq_queries = [];
end
if isempty(iq_queries)
    q = zeros(0,1);
    if ~any(arrayfun(@(s) isempty(s.points),groups))
        qlo = max(arrayfun(@(s) min(s.points(:,2)),groups));
        qhi = min(arrayfun(@(s) max(s.points(:,2)),groups));
        if qlo <= qhi
            q = linspace(qlo,qhi,201)';
            if qlo <= 0 && qhi >= 0
                q = unique([q;0]);
            end
        end
    end
else
    q = iq_queries(:);
end
roots = cell(numel(q),2);
counts = zeros(numel(q),2);
paired = nan(numel(q),2,3);
for i = 1:numel(q)
    for j = 1:2
        roots{i,j} = fixed_iq_roots(groups(j),q(i));
        counts(i,j) = size(roots{i,j},1);
        if counts(i,j) == 1
            paired(i,j,:) = roots{i,j};
        end
    end
end
fixed_delta = reshape(paired(:,2,:)-paired(:,1,:),numel(q),3);
fixed.q = q;
fixed.roots = roots;
fixed.root_counts = counts;
fixed.unique_pair = all(counts==1,2);
fixed.signed_k = -fixed_delta(:,1)./fixed_delta(:,3);
fixed.absolute_numerator_k = abs(fixed_delta(:,1))./fixed_delta(:,3);
fixed.currents = paired;
end

function roots = fixed_iq_roots(contour,q)
a = contour.a;
b = contour.b;
u = (q-a(:,2))./(b(:,2)-a(:,2));
crossing = u>=0 & u<=1;
roots = a(crossing,:) + u(crossing).*(b(crossing,:)-a(crossing,:));
% A horizontal segment represents a continuum of roots: retain both ends,
% so the slice is explicitly ambiguous and never enters an automatic pair.
horizontal = a(:,2)==q & b(:,2)==q;
roots = [roots;a(horizontal,:);b(horizontal,:)];
if ~isempty(roots)
    roots(:,2) = q;
    roots = uniquetol(roots,1e-10,'ByRows',true,'DataScale',1);
    roots = sortrows(roots,1);
end
end
