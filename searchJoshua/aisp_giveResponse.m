function resp = aisp_giveResponse(type, ParamStruct, Data, percepts)
% Make a response assuming that the trial is not the result of a lapse

% INPUT
% type: Which model to simulate with?

% TODO add a check to ensure mu_s = 0 as have passed to all the functions below

% What response is given in each case?
kappaX = exp(ParamStruct.LnKappa_x);
relKappaX = kappaX(Data.SetSizeCond);

if strcmp(type, 'bayes')
    d = aisp_computeBaysianDV(percepts, Data.SetSize, relKappaX, ...
        Data.KappaS, 0);

elseif strcmp(type, 'PE')
    d = aisp_computePointEstDV(percepts, Data.SetSize, relKappaX, ...
        Data.KappaS, 0);
        
elseif strcmp(type, 'PE2')
    d = aisp_computeOptimalPointEstDV(percept, Data.SetSize, relKappaX, ...
        Data.KappaS, 0);
    
end


% Compute the probability of a target present response 
pPresent = (ParamStruct.LapseRate / 2) + ((1 - ParamStruct.LapseRate) * ...
    (1 ./ (1 + exp(-ParamStruct.Beta0 - (ParamStruct.Beta1 * d)))));

% Simulate a response
if size(pPresent, 2) ~=1; error('Bug'); end
resp = rand(size(pPresent)) < pPresent;