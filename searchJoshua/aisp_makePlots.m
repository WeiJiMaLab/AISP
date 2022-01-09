function DSet = aisp_makePlots(dataDir, parsDir, figDir, firstOrAll, ...
    Config, varargin)

% INPUT
% dataDir: Directory containing the original unfitted dataset
% parsDir: Directory containing the fitted parameter information
% figDir: Where to save the resulting figures
% firstOrAll: 'first' or 'all'. Collect the parameters associated with the first
% round of fitting only (up to Config.Nreps), or collect all parameters
% including any that were later scheduled using cluster_fcn_fancy
% Config: struct. Has the following fields...
%   ModelLabel: Cell array of labels to use for each model
%   ModelList: Cell array of names of the models, as they were named during
%       fitting
%   Nreps: The minimim number of repitions for each model that was 
%       used during fitting.
% varargin: Boolean. Override error if there appear to be duplicate fits? Use with
% caution, probably indicates a bug.

addReqPaths()
addpath('./plotFuns')
addpath('./mat-comp-model-tools')

if ~isempty(varargin)
    overrideDupError = varargin{1};
else
    overrideDupError = false;
end
% 
% accumulate_pars(dataDir, parsDir, firstOrAll, Config, overrideDupError)
% 
% %% Comparison between models
% Figures = plot_likelihoods(dataDir, parsDir, Config);
% figure(Figures.Likelihoods);
% mT_exportNicePdf(8, 6.5, figDir, 'modelComparison')
% 
% figure(Figures.Stds);
% mT_exportNicePdf(8, 6.5, figDir, 'modelStandardDeviations')
% 
% 
% %% Look at how close different runs of the same fit ended
% 
% [DSet, ~] = getData(dataDir);
% DSet = convertToDSetFormat(DSet, parsDir, Config);
% 
% mT_plotFitEndPoints(DSet, false, 2)
% 

%% Comparison between models and data

[bigPlot, indiviudalPlots] = plotAllModelFits(dataDir, parsDir, Config);

% Models together
figure(bigPlot)
scale = 1.1;
mT_exportNicePdf(20.5*scale, 10.5*scale, figDir, 'allModelFits', true);

% Models individually
for iM = 1 : length(Config.ModelList)
    figure(indiviudalPlots{iM});
    mT_exportNicePdf(15.9, 7, figDir, ['hitAndFa_model' num2str(iM)], true);
end



