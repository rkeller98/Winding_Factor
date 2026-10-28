function winding_factor_app(input_file)
% File picker for the original contour-factor analysis and all its plot tabs.
app_folder = fileparts(mfilename('fullpath'));
if ~isdeployed
    root = fileparts(app_folder);
    addpath(fullfile(root, 'matlab', 'common'), fullfile(root, 'matlab', 'k_factor'));
end

plots = gobjects(0);
factors = {};
data = {};
source_file = '';
busy = false;
last_folder = pwd;

app = figure('Name', 'Winding Factor Analysis', 'NumberTitle', 'off', ...
    'Tag', 'winding_factor_app', 'MenuBar', 'none', 'ToolBar', 'none', ...
    'Color', [.2 .2 .2], 'Position', [200 160 800 590], ...
    'CloseRequestFcn', @close_app);
movegui(app, 'center');
banner = axes(app, 'Units', 'normalized', 'Position', [0 .34 1 .66]);
image(banner, imread(fullfile(app_folder, 'assets', 'splash.png')));
axis(banner, 'image', 'off');
uicontrol(app, 'Style', 'text', 'String', 'WINDING FACTOR ANALYSIS', ...
    'Units', 'normalized', 'Position', [.04 .285 .92 .05], ...
    'BackgroundColor', [.2 .2 .2], 'ForegroundColor', [1 .7 .15], ...
    'FontSize', 19, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
select_button = uicontrol(app, 'Style', 'pushbutton', 'String', 'Select Flux Eval', ...
    'Tag', 'select_flux_eval', 'Units', 'normalized', 'Position', [.04 .19 .29 .07], ...
    'BackgroundColor', [1 .7 .15], 'FontSize', 13, 'FontWeight', 'bold', ...
    'Callback', @select_file);
save_button = uicontrol(app, 'Style', 'pushbutton', 'String', 'Save results...', ...
    'Units', 'normalized', 'Position', [.36 .19 .25 .07], ...
    'Enable', 'off', 'FontSize', 11, 'Callback', @save_results);
status = uicontrol(app, 'Style', 'text', ...
    'String', 'Select a FluxEval MAT file. Analysis and visualization start automatically.', ...
    'Tag', 'analysis_status', 'Units', 'normalized', 'Position', [.04 .025 .92 .14], ...
    'BackgroundColor', [.2 .2 .2], 'ForegroundColor', [.94 .94 .94], ...
    'FontSize', 11, 'HorizontalAlignment', 'left');
drawnow;

if nargin > 0 && ~isempty(input_file)
    load_file(char(input_file));
end
if isdeployed && isgraphics(app)
    uiwait(app);
end

    function select_file(~, ~)
        [name, folder] = uigetfile({'*.mat', 'FluxEval MAT file (*.mat)'}, ...
            'Select Flux Eval', last_folder);
        if isequal(name, 0), return; end
        load_file(fullfile(folder, name));
    end

    function load_file(file)
        if busy, return; end
        busy = true;
        select_button.Enable = 'off';
        save_button.Enable = 'off';
        app.Pointer = 'watch';
        cleanup = onCleanup(@finish_loading);
        try
            status.String = sprintf('Reading FluxEval...\n%s', file);
            drawnow;
            imported = prepare_reduction_data(file, 'eval');
            if isempty(imported.id_meas)
                error('WindingFactor:EmptyInput', 'The FluxEval contains no finite measured currents and fluxes.');
            end
            speeds = unique(imported.rpm_ref);
            if any(~isfinite(speeds))
                error('WindingFactor:MissingSpeed', 'The FluxEval contains no valid speed references.');
            end
            new_factors = cell(numel(speeds), 1);
            new_data = cell(numel(speeds), 1);
            summaries = cell(numel(speeds), 1);
            [last_folder, name] = fileparts(file);
            for index = 1:numel(speeds)
                status.String = sprintf('Calculating contours and plots: %g rpm (%d/%d)...', ...
                    speeds(index), index, numel(speeds));
                drawnow;
                rows = imported.rpm_ref == speeds(index);
                dataset = imported;
                fields = fieldnames(dataset);
                for field_index = 1:numel(fields)
                    field = fields{field_index};
                    if isnumeric(dataset.(field))
                        dataset.(field) = dataset.(field)(rows);
                    end
                end
                [currents, flux, groups, excitation] = reduction_inputs(dataset);
                if numel(excitation) < 2
                    error('WindingFactor:TooFewGroups', ...
                        '%g rpm: at least two excitation-current groups are required.', speeds(index));
                end
                if any(~isfinite(excitation))
                    error('WindingFactor:MissingExcitation', 'The FluxEval contains invalid excitation references.');
                end
                [k, contour, fig] = estimate_contour_factor(currents, flux, groups, true, excitation);
                plots(end+1) = fig;
                if ~isfinite(k)
                    delete(fig);
                    error('WindingFactor:NoContours', ...
                        '%g rpm: no shared contours in the 25-75%% flux range; k cannot be determined.', speeds(index));
                end
                fig.Name = sprintf('%s | %g rpm | k = %.6f A/A', name, speeds(index), k);
                add_pointwise_factor_views(fig, currents, flux, groups, excitation, contour.reduction, index);
                new_factors{index} = struct('k', k, 'contour', contour, 'rpm', speeds(index));
                new_data{index} = dataset;
                summaries{index} = sprintf('%g rpm: k = %.8f A/A | %d points | %d excitation groups', ...
                    speeds(index), k, size(currents, 1), numel(excitation));
                fig.UserData = new_factors{index};
                drawnow;
            end
            factors = new_factors;
            data = new_data;
            source_file = file;
            status.String = sprintf('%s\n%s', [name '.mat'], strjoin(summaries, newline));
            figure(plots(end));
        catch exception
            status.String = sprintf('Analysis failed:\n%s', exception.message);
            errordlg(exception.message, 'FluxEval analysis', 'non-modal');
        end
    end

    function finish_loading()
        busy = false;
        if ~isgraphics(app), return; end
        app.Pointer = 'arrow';
        select_button.Enable = 'on';
        if ~isempty(factors), save_button.Enable = 'on'; end
    end

    function save_results(~, ~)
        [~, name] = fileparts(source_file);
        [file, folder] = uiputfile('*.mat', 'Save winding-factor results', ...
            fullfile(last_folder, [name '_winding_factors.mat']));
        if isequal(file, 0), return; end
        try
            save(fullfile(folder, file), 'factors', 'data', 'source_file', '-v7.3');
        catch exception
            errordlg(exception.message, 'Save results', 'non-modal');
        end
    end

    function close_app(~, ~)
        if busy, return; end
        delete(plots(isgraphics(plots)));
        delete(app);
    end
end
