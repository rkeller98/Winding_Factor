function write_flux_function_tables
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'common'));
root = setup_matlab_paths('function_model');
paper = fullfile(root,'function_model_study');
out = fullfile(root,'results','flux_fits');
require_analysis_file(fullfile(out,'data.mat'),'prepare_flux_fit_data');
require_analysis_file(fullfile(out,'comparison.mat'),'compare_flux_fits');
require_analysis_file(fullfile(out,'metrics.csv'),'compare_flux_fits');
require_analysis_file(fullfile(out,'model_sizes.csv'),'compare_flux_fits');
require_analysis_file(fullfile(out,'runtime.csv'),'compare_flux_fits');
require_analysis_file(fullfile(out,'derivatives.csv'),'compare_flux_fits');
require_analysis_file(fullfile(out,'differences_from_LUT.csv'),'compare_fits_with_k_factor');
a = load(fullfile(root,'results','flux_fits','data.mat'));
data = a.data;
current_sets = find(cellfun(@(d) strcmp(d.active{3},'ie'), data));
a = load(fullfile(root,'results','flux_fits','comparison.mat'));
study = a.study;
labels = a.labels;
M = readtable(fullfile(out,'metrics.csv'));
S = readtable(fullfile(out,'model_sizes.csv'));
C = readtable(fullfile(out,'runtime.csv'));
L = readtable(fullfile(out,'derivatives.csv'));
D = readtable(fullfile(out,'differences_from_LUT.csv'));
f = fopen(fullfile(paper,'results.tex'),'w');
fprintf(f,'\\newcommand{\\PredictionTable}{\\begin{longtable}{llrrrrr}\\toprule Set & Model & d & q & Joint & Points & Coverage\\\\\\midrule\n');

for s = 1:numel(data)
    for m = 1:numel(study{s}.models)
        a = M(M.dataset == s & M.mode == 2 & M.group == 0 & M.method == m, :);
        fprintf(f,'%d & %s & %.3f & %.3f & %.3f & %d & %.1f\\%%\\\\\n', ...
            s, labels{m}, a.rmse_d_mVs, a.rmse_q_mVs, a.joint_rmse_mVs, a.valid_points, 100*a.coverage);
    end
end
fprintf(f,'\\bottomrule\\end{longtable}}\n');
fprintf(f,'\\newcommand{\\CostTable}{\\begin{longtable}{llrrrr}\\toprule Set & Model & Terms & Coefficients & KiB & PC median in $\\mu$s\\\\\\midrule\n');

for s = 1:numel(data)
    for m = 1:numel(study{s}.models)
        a = S(S.dataset == s & S.method == m, :);
        b = C(C.dataset == s & C.method == m, :);
        fprintf(f,'%d & %s & %d & %d & %.2f & %.1f\\\\\n', s, labels{m}, a.terms, a.flux_coefficients, a.online_numeric_bytes/1024, b.median_us);
    end
end
fprintf(f,'\\bottomrule\\end{longtable}}\n');
fprintf(f,'\\newcommand{\\DerivativeTable}{\\begin{longtable}{llrrrrrr}\\toprule Set & Model & $L_{dd}$ & $L_{dq}$ & $L_{de}$ & $L_{qd}$ & $L_{qq}$ & $L_{qe}$\\\\\\midrule\n');

for s = current_sets(:)'
    for m = [1 2 3 4 9]
        fprintf(f,'%d & %s', s, labels{m});

        for c = 1:2
            for r = 1:3
                a = L(L.dataset == s & L.method == m & L.output == c & L.input == r, :);
                fprintf(f,' & %.3f', a.rmse_scaled);
            end
        end
        fprintf(f,'\\\\\n');
    end
end
fprintf(f,'\\bottomrule\\end{longtable}}\n');
fprintf(f,'\\newcommand{\\DifferenceTable}{\\begin{tabular}{lrrrr}\\toprule Set & RBF & Spline & kA & kF\\\\\\midrule\n');

for s = current_sets(:)'
    fprintf(f,'%d', s);

    for m = [1 2 4 9]
        a = D(D.dataset == s & D.method == m, :);
        fprintf(f,' & %.3f', a.joint_rmse_mVs);
    end
    fprintf(f,'\\\\\n');
end
fprintf(f,'\\bottomrule\\end{tabular}}\n');
fprintf(f,'\\newcommand{\\ReconstructionTable}{\\begin{tabular}{lrrrrr}\\toprule Set & RBF & Spline & LUT & kA & kF\\\\\\midrule\n');

for s = current_sets(:)'
    fprintf(f,'%d', s);

    for m = [1 2 3 4 9]
        a = M(M.dataset == s & M.mode == 1 & M.group == 0 & M.method == m, :);
        fprintf(f,' & %.3f', a.joint_rmse_mVs);
    end
    fprintf(f,'\\\\\n');
end
fprintf(f,'\\bottomrule\\end{tabular}}\n');
fprintf(f,'\\newcommand{\\ResultSummary}{Joint prediction RMSE (mVs): ');

for s = current_sets(:)'
    fprintf(f,'Set %d: ', s);

    for m = 1:3
        a = M(M.dataset == s & M.mode == 2 & M.group == 0 & M.method == m, :);
        fprintf(f,'%s %.3f; ', labels{m}, a.joint_rmse_mVs);
    end
end
fprintf(f,'Upper-boundary groups are evaluated without extrapolation.}\n');
fclose(f);
end
