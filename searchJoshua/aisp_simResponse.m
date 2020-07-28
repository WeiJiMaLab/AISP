function resp = aisp_simResponse(type, ParamStruct, Data)
% Simulate a response assuming that the trial is not the result of a lapse

% INPUT
% type: Which model to simulate with?

%% Simulate percepts

% For all models we need to add noise to the real orientations to generate
% simulated percepts
kappaX = exp(ParamStruct.LnKappa_x);
relKappaX = kappaX(Data.SetSizeCond);
percepts = aisp_addNoiseToStim(relKappaX, Data.Orientation);


%% Simulate responses based on these percepts

% What response is given in each case?
kappaX = exp(ParamStruct.LnKappa_x);
relKappaX = kappaX(Data.SetSizeCond);

if strcmp(type, 'Bayes')
    d = aisp_computeBaysianDV(percepts, Data.SetSize, relKappaX, ...
        Data.KappaS, 0);

elseif strcmp(type, 'PE')
    d = aisp_computePointEstDV(percepts, Data.SetSize, relKappaX, ...
        Data.KappaS, 0);
        
elseif strcmp(type, 'PE2')
    d = aisp_computeOptimalPointEstDV(percepts, Data.SetSize, relKappaX, ...
        Data.KappaS, 0);
    
end


% Compute the probability of a target present response 
pPresent = (ParamStruct.LapseRate / 2) + ((1 - ParamStruct.LapseRate) * ...
    (1 ./ (1 + exp(-ParamStruct.Beta0 - (ParamStruct.Beta1 * d)))));

% Simulate a response
if size(pPresent, 2) ~=1; error('Bug'); end
resp = rand(size(pPresent)) < pPresent;


if any(isnan(resp(:))); error('Bug'); end