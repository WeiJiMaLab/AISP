function resp = aisp_simResponse(modelName, ParamStruct, Data)
% Simulate a response

% INPUT
% type: Which model to simulate with?

if any(strcmp(modelName, {'impSamp', 'jointPostSamp'}))
    assert(isequal(size(ParamStruct.NumSamples), [1, 1]))
    ParamStruct.NumSamples = randRound(ParamStruct.NumSamples);
end


%% Simulate percepts

% For all models we need to add noise to the real orientations to generate
% simulated percepts
kappaX = exp(ParamStruct.LnKappa_x);
relKappaX = kappaX(Data.SetSizeCond);
percepts = aisp_addNoiseToStim(relKappaX, Data.Orientation, ...
                                'efficientSamp');


%% Simulate responses based on these percepts

% What response is given in each case?
if strcmp(modelName, 'Bayes')
    d = aisp_computeBaysianDV(percepts, Data.SetSize, relKappaX, ...
        Data.KappaS, 0);

elseif strcmp(modelName, 'PE')
    d = aisp_computePointEstDV(percepts, Data.SetSize, relKappaX, ...
        Data.KappaS, 0, 'stimAndTarg');
    
elseif strcmp(modelName, 'PE_imagineL')
    d = aisp_computePointEstDV(percepts, Data.SetSize, relKappaX, ...
        Data.KappaS, 0, 'stimAndTarg_incImagine');
        
elseif strcmp(modelName, 'PE2')
    d = aisp_computeOptimalPointEstDV(percepts, Data.SetSize, relKappaX, ...
        Data.KappaS, 0);
    
elseif strcmp(modelName, 'impSamp')
    d = aisp_computeImpSampDV(percepts, Data.SetSize, relKappaX, ...
        Data.KappaS, 0, ParamStruct.NumSamples);
    
elseif strcmp(modelName, 'jointPostSamp')
    d = aisp_computeJointPostSampDV(percepts, Data.SetSize, relKappaX, ...
        Data.KappaS, 0, ParamStruct.NumSamples, true);
else
    error('Bug')
end


% Compute the probability of a target present response 
pPresent = (ParamStruct.LapseRate / 2) + ((1 - ParamStruct.LapseRate) * ...
    (1 ./ (1 + exp(-ParamStruct.Beta0 - (ParamStruct.Beta1 * d)))));

% Simulate a response
if size(pPresent, 2) ~=1; error('Bug'); end
resp = rand(size(pPresent)) < pPresent;


if any(isnan(resp(:))); error('Bug'); end