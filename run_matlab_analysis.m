function run_matlab_analysis
% Run both studies in order; each study resolves its own common dependencies.
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'matlab','k_factor'),fullfile(root,'matlab','function_model'));
fprintf('\n[1/2] k-factor and coordinate-reduction study\n');
run_k_factor_analysis;
fprintf('\n[2/2] function-model study\n');
run_function_model_analysis;
drawnow;
fprintf('\nBoth studies finished. Results: %s\nAll analysis figures remain open.\n',fullfile(root,'results'));
end
