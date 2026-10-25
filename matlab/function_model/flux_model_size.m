function [terms, coefficients, bytes] = flux_model_size(m)
% Mandatory double-precision evaluation arrays; exclude provenance and work arrays.
switch m.kind
    case 'rbf'
        terms = size(m.centers, 1);
        coefficients = numel(m.coefficients);
        count = numel(m.centers)+1+coefficients+numel(m.origin)+numel(m.scale)+numel(m.bounds);
    case 'bspline'
        terms = size(m.coefficients, 1);
        coefficients = numel(m.coefficients);
        count = coefficients+sum(cellfun(@numel, m.knots))+1+numel(m.origin)+numel(m.scale)+numel(m.bounds);
    otherwise
        a = m.map;
        terms = a.nodes;
        coefficients = 2*a.nodes;
        count = coefficients+sum(cellfun(@numel, a.axes))+numel(a.p)+numel(a.current_scale)+ ...
            numel(a.origin)+numel(a.scale)+numel(m.bounds);
end
bytes = 8*count;
end
