function export_measurements_json(datasetRoot)
% Add JSON copies of the raw and extracted MAT files; keep the MAT files.
if nargin == 0
    datasetRoot = fileparts(mfilename('fullpath'));
end
patterns = {'Set*/raw/*.mat', 'Set*/extracted/raw/*.mat', ...
    'Set*/extracted/evaluation/*.mat'};
for p = 1:numel(patterns)
    files = dir(fullfile(datasetRoot, patterns{p}));
    for k = 1:numel(files)
        source = fullfile(files(k).folder, files(k).name);
        [~, name] = fileparts(source);
        output = fullfile(files(k).folder, [name '.json']);
        bytes = unicode2native(jsonencode(load(source)), 'UTF-8');
        fid = fopen(output, 'w');
        fwrite(fid, bytes, 'uint8');
        fclose(fid);
        gzip(output);
        delete(output);
        fprintf('%s.gz: %d JSON bytes\n', output, numel(bytes));
    end
end
end
