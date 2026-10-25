function [outputFile, data] = extract_meas_from_fluxeval(fluxEvalPath, outputFile)
% All Meas fields as MAT: grouped Total and cells {excitation current,speed}.
% Example: f = extract_meas_from_fluxeval('C:\path\FluxEval.mat'); data = load(f);
input = load(fluxEvalPath, 'Total');
rpm_fields = fieldnames(input.Total);
rpm_fields = rpm_fields(endsWith(rpm_fields, '_rpm'));
speeds = zeros(size(rpm_fields));
for r = 1:numel(rpm_fields)
    group = input.Total.(rpm_fields{r});
    names = fieldnames(group);
    names = names(endsWith(names, '_A'));
    speeds(r) = double(group.(names{1}).Ref_Speed);
end
[~, order] = sort(speeds);
rpm_fields = rpm_fields(order);
Total = struct;
Meas = struct;
rpm_ref = {};
exc_ref = {};
exc_fields = {};
for r = 1:numel(rpm_fields)
    group = input.Total.(rpm_fields{r});
    names = fieldnames(group);
    names = names(endsWith(names, '_A'));
    currents = cellfun(@(name) double(group.(name).Ref_exc_current), names);
    [~, order] = sort(currents);
    names = names(order);
    for c = 1:numel(names)
        slice = group.(names{c});
        fields = fieldnames(slice);
        slice = rmfield(slice, fields(~strcmp(fields, 'Meas') & ~startsWith(fields, 'Ref_')));
        Total.(rpm_fields{r}).(names{c}) = slice;
        rpm_ref{c,r} = slice.Ref_Speed;
        exc_ref{c,r} = slice.Ref_exc_current;
        exc_fields{c,r} = names{c};
        fields = fieldnames(slice.Meas);
        for f = 1:numel(fields)
            Meas.(fields{f}){c,r} = slice.Meas.(fields{f});
        end
    end
end
data.Total = Total;
data.Meas = Meas;
data.rpm_fields = rpm_fields;
data.exc_fields = exc_fields;
data.rpm_ref = rpm_ref;
data.exc_ref = exc_ref;
if nargin < 2 || isempty(outputFile)
    [folder, name] = fileparts(fluxEvalPath);
    [~, parent] = fileparts(folder);
    if strcmpi(parent, 'FluxEval')
        folder = fullfile(fileparts(folder), 'extracted', 'evaluation');
    end
    outputFile = fullfile(folder, [name '_measurements.mat']);
end
folder = fileparts(outputFile);
if ~isempty(folder) && ~exist(folder, 'dir')
    mkdir(folder);
end
save(outputFile, '-struct', 'data', '-v7');
end
