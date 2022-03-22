function standalonePerformancePlotter()
% Make the simulated performance plots. Directories to add to path, and
% directories from which to load data are specified relative to the 
% directory AISP/searchJoshua. Therefore this function must be called 
% while AISP/searchJoshua is the current working directory.

addReqPaths()
addpath('./plotFuns')
addpath('./mat-comp-model-tools')
addpath('./simFuns')

dataDir = fullfile('Data','StandardFormat_participantExcluded.mat');
parsDir = fullfile('Data','Fits','10_nyuRun_6thModel');
figDir = 'figures';

Config = produceConfig(true);

aisp_plotSimPerformance(dataDir, parsDir, Config.ModelList, ...
    Config.ModelLabel);
%mT_exportNicePdf(8*0.8, 10*0.8, figDir, 'simPerformance')