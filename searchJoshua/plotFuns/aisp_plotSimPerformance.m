function fig = aisp_plotSimPerformance(dataDir, parsDir, ModelList, ...
    ModelLabels, varargin)
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
% varargin{1}: scalar. Sets number of items to use in the simulated 
%   stimuli. Determines which one of four options to use, does not directly
%   determine the number of stimuli.
% varargin{2}: scalar. Sets the concentration parameter of the distractors 
%   in the simulated stimuli. Determines which one of two options to use, 
%   does not directly determine the concentration parameter.
% varargin{3}: scalar. Determines how many trials to use for simulating
%   each indiviudal point in the resulting plots.


if (length(varargin) > 0) && (~isempty(varargin{1}))
    iSetSize = varargin{1};
else
    iSetSize = 2;
end

if (length(varargin) > 1) && (~isempty(varargin{2}))
    iBlockType = varargin{2};
else
    iBlockType = 1;
end

if (length(varargin) > 2) && (~isempty(varargin{3}))
    nTrialPerSim = varargin{3};
else
    nTrialPerSim = 100000;
end

numSamples = 10;
lnKappaRange = -4 : 0.2 : 7;


%% Choose params to simulate with

templateModel = 'Bayes';
bestParams = aisp_loadBestFits(dataDir, parsDir, templateModel, 'array');
[DSet, Nptpnts] = getData(dataDir);
assert(size(bestParams, 1) == Nptpnts)
meanParams = mean(bestParams, 1);


%% Run simulation

accAcrossLnKappa = nan(length(lnKappaRange), length(ModelList));

for iM = 1 : length(ModelList)
    modelName = ModelList{iM};
    
    % For this plot we fix some params to specific values
    ParamStruct = paramVec2Struct(meanParams, templateModel, 'to struct');
    OrigParamStruct = ParamStruct;
    ParamStruct.Beta0 = 0;
    ParamStruct.LapseRate = 0;
    
    if any(strcmp(modelName, {'impSamp'  'jointPostSamp'}))
        ParamStruct.NumSamples = numSamples;
    else
        assert(sum(strcmp(modelName, ...
            {'Bayes', 'PE', 'PE2', 'PE_imagineL'}))==1)
    end
    
    disp('Using the following param struct...')
    disp('Apart from LnKappa_x which will be varied')
    disp(ParamStruct)
    
    for iK = 1 : length(lnKappaRange)
        ParamStruct.LnKappa_x(:) = lnKappaRange(iK);
        
        % Stimulate stimulus
        nItems = [2, 3, 4, 6];
        setSizeCond = [1, 2, 3, 4];
        kappa_s = [0, 1.5];
        blockType = [1, 2];
        distStats.mu_s = 0;
        
        % Checks
        SomeData = DSet.P(1).Data;
        realItemCondCombos = unique([SomeData.SetSize , ...
            SomeData.SetSizeCond], 'rows');
        assert(size(realItemCondCombos, 2) == 2)
        assert(isequal([nItems', setSizeCond'], ...
            sortrows(realItemCondCombos)))
        
        thisNItems = nItems(iSetSize);
        thisSetSizeCond = setSizeCond(iSetSize);
        distStats.kappa_s = kappa_s(iBlockType);
        thisBlockType = blockType(iBlockType);
        
        Data = aisp_simSingleCondStimulus(nTrialPerSim, thisNItems, ...
            thisSetSizeCond, distStats, thisBlockType);
        Data.Response = nan(size(Data.Target));
        SimDSet = struct();
        SimDSet.P(1).Data = Data;
        
        % Simulate responses
        SimDSet = aisp_simRespAndAcc(SimDSet, modelName, {ParamStruct});
        
        assert(all(ismember(SimDSet.P(1).Data.Accuracy, [0, 1])))
        accAcrossLnKappa(iK, iM) = mean(SimDSet.P(1).Data.Accuracy);
    end
    
    disp(['Model ', num2str(iM), ' simulation complete.'])
end

% Find the fitted value for ln kappa_x for this experimental
% setting
fittedLnKappa_x = OrigParamStruct.LnKappa_x(thisSetSizeCond);

fig = figure;
kappaRange = exp(lnKappaRange);
lines = semilogx(kappaRange, accAcrossLnKappa);
set(gca, 'xdir', 'reverse')
ylabel('Proportion correct')
xlabel('Sensory precision (concentration parameter, \kappa)')
yline(0.5, '--')
hold on
semilogx(exp([fittedLnKappa_x, fittedLnKappa_x]), [0.5, 0.45], 'k')
legendLabels = ModelLabels;
legend(lines, legendLabels)
legend box off
set(gca, 'TickDir', 'out');
set(gca, 'box', 'off');
