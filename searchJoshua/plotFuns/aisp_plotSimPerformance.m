function aisp_plotSimPerformance(DSet, Config)
% Simulate performance of observers using the various models, under 
% handpicked param values, while the average noise level varies

% INPUT
% DSet: Stimuli is this dataset will be used as the basis of the
% stimulation
% Config: struct. Has the following fields...
%   ModelLabel: Cell array of labels to use for each model
%   ModelList: Cell array of names of the models, as they were named during
%       fitting


% Choose params to simulate with
params = {'LnKappa_x', 'LnKappa_x', 'LnKappa_x', 'LnKappa_x', ...
    'Beta0', 'Beta1', 'LapseRate', 'NumSamples'};
paramVec = [2.7, 2.6, 2.2, 1.7, 0, 1, 0, 10];
aisp_checkParamOrder(params, true)

lnKappaRange = -1.6 : 0.2 : 3;

for iM = 1 : length(Config.ModelList)
    modelName = Config.ModelList{iM};
    accAcrossKappa = nan(length(lnKappaRange), 1);
    
    for iK = 1 : length(lnKappaRange)
        % Get params into required format
        if any(strcmp(modelName, {'Bayes', 'PE', 'PE2', 'PE_imagineL'}))
            thisParamVec = paramVec(1:end-1);
        else
            thisParamVec = paramVec;
        end
        ParamStruct = paramVec2Struct(thisParamVec, modelName, 'to struct');
        ParamStruct.LnKappa_x = ParamStruct.LnKappa_x + lnKappaRange(iK);
        paramStructs = {ParamStruct};
        paramStructs = repmat(paramStructs, 1, length(DSet.P));
        
        nReps = 10;
        allAcc = nan(nReps, 1);
        for iRep = 1 : nReps
            SimDSet = aisp_simRespAndAcc(DSet, modelName, paramStructs);
            acc = mT_stackData(SimDSet.P, @(str)mean(str.Data.Accuracy));
            allAcc(iRep) = mean(acc);
        end
        
        accAcrossKappa(iK) = mean(allAcc);
        disp(['Finished kappa iter ', num2str(iK)])
    end
end

% THIS FUN IS WORK IN PROGRESS