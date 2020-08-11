function aisp_makePlots(dataDir, parsDir, figDir)

% INPUT
% dataDir: Directory containing the original unfitted dataset
% parsDir: Directory containing the fitted parameter information
% figDir: Where to save the resulting figures

addpath('../lautils-mat/stats')
addpath('./circstat-matlab')
addpath('./visualSearch')
addpath('./visualSearch/plots')
addpath('./visualSearch/modellingTools')

Config = load('Config.mat');

accumulate_pars(dataDir, parsDir)

%% Comparison between models
Figures = plot_likelihoods(dataDir, parsDir);
figure(Figures.Likelihoods);
mT_exportNicePdf(14, 15.9, figDir, 'modelComparison')


%% Comparison between models and data
for iM = 1 : length(Config.ModelList)
    Figures = aisp_compareModelAndData(dataDir, parsDir, Config.ModelList{iM});
    
    figure(Figures.AllSumStatsAllData)
    mT_exportNicePdf(14, 15.9, figDir, ['allDataAllStats_', Config.ModelList{iM}])
end