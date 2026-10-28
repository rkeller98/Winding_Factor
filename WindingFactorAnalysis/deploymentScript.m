function deploymentScript(include_runtime)
% Build the Windows EXE and an offline installer, following MeasEval's build.
if nargin < 1, include_runtime = true; end
app_folder = fileparts(mfilename('fullpath'));
root = fileparts(app_folder);
original_path = path;
cleanup = onCleanup(@() path(original_path));
restoredefaultpath;
addpath(app_folder, fullfile(root, 'matlab', 'common'), fullfile(root, 'matlab', 'k_factor'));

assets = fullfile(app_folder, 'assets');
options = compiler.build.StandaloneApplicationOptions(fullfile(app_folder, 'winding_factor_app.m'));
options.ExecutableName = 'WindingFactorAnalysis';
options.ExecutableVersion = '1.0.0.0';
options.OutputDir = fullfile(app_folder, 'application');
options.EmbedArchive = true;
options.AutoDetectDataFiles = false;
options.AdditionalFiles = {assets};
options.ExecutableIcon = fullfile(assets, 'icon.png');
options.ExecutableSplashScreen = fullfile(assets, 'icon.png');
options.Verbose = true;
result = compiler.build.standaloneWindowsApplication(options);

if include_runtime
    if ~isfile(mcrinstaller)
        compiler.runtime.download;
    end
    package = compiler.package.InstallerOptions(result);
    package.ApplicationName = 'Winding Factor Analysis';
    package.InstallerName = 'WindingFactorAnalysis_Setup';
    package.OutputDir = fullfile(app_folder, 'installer');
    package.RuntimeDelivery = 'installer';
    package.InstallerIcon = fullfile(assets, 'icon.png');
    package.InstallerSplash = fullfile(assets, 'icon.png');
    package.InstallerLogo = fullfile(assets, 'installation.png');
    package.Summary = 'Select a FluxEval MAT file to analyze and visualize winding factors.';
    package.Version = '1.0.0.0';
    package.Shortcut = fullfile(options.OutputDir, 'WindingFactorAnalysis.exe');
    package.Verbose = true;
    compiler.package.installer(result, 'Options', package);
end
end
