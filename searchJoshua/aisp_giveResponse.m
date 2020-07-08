function resp = aisp_giveResponse(Model, ParamStruct, Data, percepts)
% Make a response assuming that the trial is not the result of a lapse

% INPUT
% stim   [num trials X num locations] array describing the perceived orientation 
%        of the Gabor patches


% What response is given in each case?
relKappaX = ParamStruct.Kappa_x(Data.SetSizeCond);

if strcmp(Model, 'bayes')
    d = makeBaysianDecision(percepts, Data.SetSize, relKappaX, ...
        Data.KappaS, 0, ParamStruct);
    % TODO add a check to ensure mu_s = 0 as have passed to this function

elseif strcmp(Model, 'pointEst')
    d = aisp_computePointEstDV(percepts, Data.SetSize, relKappaX, ...
        Data.KappaS, 0, ParamStruct);
    % TODO add a check to ensure mu_s = 0 as have passed to this function
        
elseif strcmp(Model, 'optimalPointEst')
    
end


% Compute the probability of a target present response 
pPresent = (ParamStruct.LapseRate / 2) + ((1 - ParamStruct.LapseRate) * ...
    (1 ./ (1 + exp(-ParamStrcut.Beta0 - (ParamStrcut.Beta1 * d)))));

% Simulate a response
if size(pPresent, 2) ~=1; error('Bug'); end
resp = rand(size(pPresent)) < pPresent;