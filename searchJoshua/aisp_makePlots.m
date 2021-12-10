function DSet = aisp_makePlots(dataDir, parsDir, figDir, firstOrAll, ...
    configFile, varargin)

% INPUT
% dataDir: Directory containing the original unfitted dataset
% parsDir: Directory containing the fitted parameter information
% figDir: Where to save the resulting figures
% firstOrAll: 'first' or 'all'. Collect the parameters associated with the first
% round of fitting only (up to Config.Nreps), or collect all parameters
% including any that were later scheduled using cluster_fcn_fancy
% configFile: string. File path to matlab file to be loaded.
% varargin: Boolean. Override error if there appear to be duplicate fits? Use with
% caution, probably indicates a bug.

addpath('./plotFuns')
addpath('../lautils-mat/stats')
addpath('./circstat-matlab')
addpath('./visualSearch')
addpath('./visualSearch/plotFuns')
addpath('./visualSearch/analysisFuns')
addpath('./visualSearch/modellingTools')

if ~isempty(varargin)
    overrideDupError = varargin{1};
else
    overrideDupError = false;
end

Config = load(configFile);

accumulate_pars(dataDir, parsDir, firstOrAll, configFile, overrideDupError)

%% Comparison between models
Figures = plot_likelihoods(dataDir, parsDir, configFile);
figure(Figures.Likelihoods);
mT_exportNicePdf(4, 6.5, figDir, 'modelComparison')


%% Comparison between models and data

[bigPlot, indiviudalPlots] = plotAllModelFits(dataDir, parsDir, configFile);

% Models together
figure(bigPlot)
scale = 1.1;
mT_exportNicePdf(10.5*scale, 10.5*scale, figDir, 'allModelFits', true);

% Models individually
for iM = 1 : length(Config.ModelList)
    figure(indiviudalPlots{iM});
    mT_exportNicePdf(15.9, 7, figDir, ['hitAndFa_model' num2str(iM)], true);
end


%% Look at how close different runs of the same fit ended

[DSet, ~] = getData(dataDir);
DSet = convertToDSetFormat(DSet, parsDir, configFile);

mT_plotFitEndPoints(DSet, false, 2)


