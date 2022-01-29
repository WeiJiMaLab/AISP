function fig = aisp_plotSimPerformance(dataDir, parsDir, ModelList, ...
    ModelLabels)
% Simulate performance of observers using the various models, under 
% handpicked param values, while the average noise level varies. Note, 
% simulates a new dataset in which, for all trials, the stimulus has the
% same number of items in it.

% INPUT
% dataDir: Directory containing the original unfitted dataset
% parsDir: Directory containing the fitted parameter information. The 
%   paramters used for the simulation are based off the fitted parameters.
%   See code below.
% ModelList: Cell array of names of the models, as they were named during
%       fitting
% ModelLabels: Cell array of labels to use for each model


%% Choose params to simulate with

templateModel = 'Bayes';
bestParams = aisp_loadBestFits(dataDir, parsDir, templateModel, 'array');
[~, Nptpnts] = getData(dataDir);
assert(size(bestParams, 1) == Nptpnts)
meanParams = mean(bestParams, 1);
ParamStruct = paramVec2Struct(meanParams, templateModel, 'to struct');

% For this plot we fix some params to specific values
ParamStruct.Beta0 = 0;
ParamStruct.LapseRate = 0;
numSamples = 10;
lnKappaRange = -4 : 0.2 : 4;
kappaRange = exp(lnKappaRange);


%% Run simulation

accAcrossLnKappa = nan(length(lnKappaRange), length(ModelList));

for iM = 1 : length(ModelList)
    modelName = ModelList{iM};
    
    if any(strcmp(modelName, {'impSamp'  'jointPostSamp'}))
        ParamStruct.NumSamples = numSamples;
    else
        assert(sum(strcmp(modelName, ...
            {'Bayes', 'PE', 'PE2', 'PE_imagineL'}))==1)
    end
    
    for iK = 1 : length(lnKappaRange)
        ParamStruct.LnKappa_x(:) = lnKappaRange(iK);
        
        % Stimulate stimulus
        nItems = [2, 3, 4, 6];
        setSizeCond = [1, 2, 3, 4];
        kappa_s = [0, 1.5];
        blockType = [1, 2];
        distStats.mu_s = 0;
        
        iSetSize = 2;
        iBlockType = 1;
        thisNItems = nItems(iSetSize);
        thisSetSizeCond = setSizeCond(iSetSize);
        distStats.kappa_s = kappa_s(iBlockType);
        thisBlockType = blockType(iBlockType);
        
        Data = aisp_simSingleCondStimulus(100000, thisNItems, ...
            thisSetSizeCond, distStats, thisBlockType);
        Data.Response = nan(size(Data.Target));
        SimDSet = struct();
        SimDSet.P(1).Data = Data;
        
        % Simulate responses
        SimDSet = aisp_simRespAndAcc(SimDSet, modelName, {ParamStruct});
        
        assert(all(ismember(SimDSet.P(1).Data.Accuracy, [0, 1])))
        accAcrossLnKappa(iK, iM) = mean(SimDSet.P(1).Data.Accuracy);
    end
end

fig = figure;
lines = semilogx(kappaRange, accAcrossLnKappa);
set(gca, 'xdir', 'reverse')
ylabel('Proportion correct')
xlabel('Sensory precision (log concentration parameter)')
yline(0.5, '--')
legendLabels = ModelLabels;
legend(lines, legendLabels)
legend box off

disp('here')
