function root=setup_matlab_paths(study)
% Add only common helpers, the requested study, and the MeasEval package root.
folder=fileparts(mfilename('fullpath'));
matlab_folder=fileparts(folder);
root=fileparts(matlab_folder);
geometry=fullfile(root,'Resources','MeasEval','matlab');
if ~isfolder(geometry)
    error('WindingFactor:MissingDependency', ...
        'MeasEval MATLAB resources missing: %s. Restore Resources/MeasEval from this repository.',geometry);
end
addpath(folder,geometry);
if nargin>0
    addpath(fullfile(matlab_folder,study));
end
end
