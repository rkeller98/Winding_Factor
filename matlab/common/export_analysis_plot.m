function export_analysis_plot(plot_handle, file_name, export_format, resolution)
% Save a plot only when a format was explicitly requested.
if isempty(export_format)
    return;
end

[folder, name] = fileparts(file_name);
if ~isfolder(folder)
    mkdir(folder);
end
file_name = fullfile(folder, [name '.' char(export_format)]);

exportgraphics(plot_handle, file_name, 'Resolution', resolution);
end
