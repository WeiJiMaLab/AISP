function aisp_makePlots(dataDir)

addpath('../lautils-mat/stats')
addpath('./visualSearch')
addpath('./visualSearch/plots')
addpath('./visualSearch/modellingTools')

plotsDir = './figures';

Config = load('Config.mat');

%% Comparison between models
Figures = plot_likelihoods(dataDir);
figure(Figures.Likelihoods);
mT_exportNicePdf(14, 15.9, plotsDir, 'modelComparison')


%% Comparison between models and data
for iM = 1 : length(Config.ModelList)
    Figures = aisp_compareModelAndData(dataDir, Config.ModelList{iM});
    
    figure(Figures.AllSumStatsAllData)
    mT_exportNicePdf(14, 15.9, plotsDir, ['allDataAllStats_', Config.ModelList{iM}])
end