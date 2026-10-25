function extract_measurements(datasetRoot)
% Extract all local sets again. Originals are in raw/ and FluxEval/.
if nargin == 0
    datasetRoot = fileparts(mfilename('fullpath'));
end
addpath(fileparts(fileparts(mfilename('fullpath'))));
sets = dir(fullfile(datasetRoot, 'Set*'));
sets = sets([sets.isdir]);
for k = 1:numel(sets)
    setRoot = fullfile(datasetRoot, sets(k).name);
    files = dir(fullfile(setRoot, 'raw', '*.mat'));
    for j = 1:numel(files)
        source = fullfile(files(j).folder, files(j).name);
        folder = fullfile(setRoot, 'extracted', 'raw');
        if ~exist(folder, 'dir'), mkdir(folder); end
        [~, name] = fileparts(source);
        extract_raw(source, fullfile(folder, [name '_measurements.mat']));
        fprintf('Raw data: %s\n', source);
    end
    files = dir(fullfile(setRoot, 'FluxEval', '*.mat'));
    for j = 1:numel(files)
        extract_meas_from_fluxeval(fullfile(files(j).folder, files(j).name));
        fprintf('Evaluation: %s\n', files(j).name);
    end
end
end

function extract_raw(source, destination)
input = load(source);
names = fieldnames(input);
raw = input.(names{1});
names = strrep({raw.Y.Name}, ' ', '');
rpm = raw.Y(strcmp(names, 'PE_Rpm_req')).Data;
exc = raw.Y(strcmp(names, 'I_exc_req')).Data;
index = raw.Y(strcmp(names, 'PE_Index_meas')).Data;
Total = struct;
for speed = reshape(unique(rpm), 1, [])
    speedName = axis_name(speed, 'rpm');
    for current = reshape(unique(exc(rpm == speed)), 1, [])
        currentName = axis_name(current, 'A');
        rows = find(rpm == speed & exc == current);
        [~, order] = sort(index(rows));
        rows = rows(order);
        Meas = struct;
        for f = 1:numel(names)
            values = raw.Y(f).Data;
            Meas.(names{f}) = reshape(values(rows), 1, []);
        end
        slice.Meas = Meas;
        slice.Ref_Speed = speed;
        slice.Ref_exc_current = current;
        slice.SampleIndex = reshape(rows, 1, []);
        slice.Time = reshape(raw.X.Data(rows), 1, []);
        Total.(speedName).(currentName) = slice;
    end
end
Signals = rmfield(raw.Y, 'Data');
save(destination, 'Total', 'Signals', '-v7');
end

function name = axis_name(value, unit)
name = sprintf('pos_%g_%s', value, unit);
if value < 0, name = sprintf('neg_%g_%s', -value, unit); end
name = strrep(name, '.', 'p');
end
