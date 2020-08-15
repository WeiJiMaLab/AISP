function aisp_makePlots(dataDir, parsDir, figDir, firstOrAll)

% INPUT
% dataDir: Directory containing the original unfitted dataset
% parsDir: Directory containing the fitted parameter information
% figDir: Where to save the resulting figures
% firstOrAll: 'first' or 'all'. Collect the parameters associated with the first
% round of fitting only (up to Config.Nreps), or collect all parameters
% including any that were later scheduled using cluster_fcn_fancy

addpath('./plotFuns')
addpath('../lautils-mat/stats')
addpath('./circstat-matlab')
addpath('./visualSearch')
addpath('./visualSearch/plotFuns')
addpath('./visualSearch/analysisFuns')
addpath('./visualSearch/modellingTools')

Config = load('Config.mat');

accumulate_pars(dataDir, parsDir, firstOrAll)

%% Comparison between models
Figures = plot_likelihoods(dataDir, parsDir);
figure(Figures.Likelihoods);
mT_exportNicePdf(5, 5, figDir, 'modelComparison')


%% Comparison between models and data

[bigPlot, indiviudalPlots] = plotAllModelFits(dataDir, parsDir);

% Models together
figure(bigPlot)
mT_exportNicePdf(9, 10.5, figDir, 'allModelFits', true);

% Models individually
for iM = 1 : length(Config.ModelList)
    figure(indiviudalPlots{iM});
    mT_exportNicePdf(15.9, 7, figDir, ['hitAndFa_model' num2str(iM)], true);
end


%% Look at how close different runs of the same fit ended

[DSet, ~] = getData(dataDir);
DSet = convertToDSetFormat(DSet, parsDir);

mT_plotFitEndPoints(DSet, true, 2)


