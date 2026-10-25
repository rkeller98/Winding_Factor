function [study, figures] = run_k_factor_investigation(set_ids, export_format)
% One tabbed figure per measured set, one ideal reference and one overview.
% No export unless explicitly requested, e.g. run_k_factor_investigation(1,'png').
if nargin < 1 || isempty(set_ids)
    set_ids = 1:3;
end
if nargin < 2
    export_format = '';
end
[~,ideal_figure] = investigate_ideal_machine(export_format);
local = investigate_local_factors(set_ids);
study = cell(size(set_ids));
figures = gobjects(numel(set_ids)+2,1);
figures(1) = ideal_figure;
for j = 1:numel(set_ids)
    s = set_ids(j);
    [contour,figures(j+1)] = investigate_contour_averaging(s,[],50,export_format);
    symmetry = investigate_flux_symmetry(s,contour.group_indices,export_format);
    resistance = investigate_resistance_sensitivity(s,[],contour.group_indices,export_format);
    for factor = {'psi_d','psi_q','ls','kd','kq'}
        plot_local_factor(local{j},factor{1},[],[],export_format);
    end
    plot_local_factor_comparison(local{j},[],export_format);
    study{j} = struct('set',s,'contour',contour,'symmetry',symmetry, ...
        'resistance',resistance,'local',local{j});
end
figures(end) = plot_investigation_comparison(study,export_format);
end
