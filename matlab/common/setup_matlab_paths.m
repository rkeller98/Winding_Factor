function root = setup_matlab_paths(study)
% Add common helpers, the study's dependencies, and the MeasEval package root.

common_folder = fileparts(mfilename('fullpath'));
matlab_folder = fileparts(common_folder);
root = fileparts(matlab_folder);
geometry_folder = fullfile(root, 'Resources', 'MeasEval', 'matlab');

if ~isfolder(geometry_folder)
    error('WindingFactor:MissingDependency', ...
        ['MeasEval MATLAB resources missing: %s. ' ...
         'Restore Resources/MeasEval from this repository.'], geometry_folder);
end

addpath(common_folder, geometry_folder);

if nargin > 0
    addpath(fullfile(matlab_folder, study));

    if strcmp(study, 'function_model')
        addpath(fullfile(matlab_folder, 'k_factor'));
    end
end
end
