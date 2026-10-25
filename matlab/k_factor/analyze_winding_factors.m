function analyze_winding_factors(export_format)
% Calculate the contour factor and display its field reconstruction.
% Image export is optional: analyze_winding_factors('png') or ('pdf').
if nargin < 1
    export_format = '';
end

matlab_folder = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(matlab_folder, 'common'));
root = setup_matlab_paths('k_factor');

output_folder = fullfile(root, 'results');
figure_folder = fullfile(output_folder, 'contour_factor');

input_file = fullfile(output_folder, 'reduction_data.mat');
require_analysis_file(input_file, 'run_matlab_analysis([],[],''prepare'')');
loaded = load(input_file, 'data');
data = loaded.data;
factors = cell(size(data));

graphics_properties = {'defaultAxesFontSize', 'defaultTextFontSize'};
graphics_defaults = get(groot, graphics_properties);
restore_graphics = onCleanup(@() set(groot, graphics_properties, graphics_defaults));
set(groot, 'defaultAxesFontSize', 14, 'defaultTextFontSize', 14);

for dataset_index = 1:numel(data)
    [currents, flux, groups, excitation] = reduction_inputs(data{dataset_index});
    [k, contour, fig] = estimate_contour_factor(currents, flux, groups, true, excitation);
    factors{dataset_index} = struct('k', k, 'contour', contour);
    fig.Name = sprintf('Set%d: contour-factor analysis', dataset_index);

    tab_group = findobj(fig, 'Type', 'uitabgroup', 'Tag', 'contour_tabs');
    selected_tab = tab_group.SelectedTab;
    layouts = gobjects(0);
    if ~isempty(export_format)
        layouts = findall(fig, 'Type', 'tiledlayout');
    end

    for layout_index = 1:numel(layouts)
        layout = layouts(layout_index);
        plot_tab = ancestor(layout, 'uitab');

        % Activate the plot tab and its parent tab before exporting its layout.
        container = plot_tab;
        while ~isequal(container, fig)
            if isa(container, 'matlab.ui.container.Tab')
                container.Parent.SelectedTab = container;
            end
            container = container.Parent;
        end

        drawnow;
        file_name = sprintf('%s_set%d.png', plot_tab.Tag, dataset_index);
        figure_file = fullfile(figure_folder, file_name);
        export_analysis_plot(layout, figure_file, export_format, 150);
    end

    tab_group.SelectedTab = selected_tab;

    % Keep the manuscript's existing contour-figure path available.
    if strcmpi(export_format, 'png')
        copyfile(fullfile(figure_folder, sprintf('contours_set%d.png', dataset_index)), ...
            fullfile(root, 'winding_factor_study', 'figures', ...
            sprintf('contours_set%d.png', dataset_index)));
    end

    fprintf('Set%d conventional factor: %.12f A/A\n', dataset_index, k);
end

save(fullfile(output_folder, 'winding_factors.mat'), 'factors');
if ~isempty(export_format)
    fprintf('Contour-factor figures: %s\n', figure_folder);
end
end
