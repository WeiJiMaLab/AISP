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

for iM = 1 : length(Config.ModelList)
    modelName = Config.ModelList{iM};

    % Get params into required format
    ParamStruct = paramVec2Struct(paramVec, modelName, 'toStruct');
    paramStructs = {ParamStruct};
    paramStructs = repmapt(paramStructs, 1, length(DSet.P));
    
    SimDSet = aisp_simRespAndAcc(DSet, modelName, paramStructs);
end

% THIS FUN IS WORK IN PROGRESS