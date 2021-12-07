function resp = aisp_simResponse(modelName, ParamStruct, Data, varargin)
% Simulate a response

% INPUT
% modelName: str. Which model to simulate with?
% varargin{1}: bool. Run checks? If false, does not run some checks, 
%   therefore potentially saving time, but with less chance to detect 
%   bugs. Default is true.

if length(varargin) >= 1
    runChecks = varargin{1};
else
    runChecks = true;
end

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
    
elseif any(strcmp(modelName, {'impSamp', 'jointPostSamp'}))
    uniqSizes = unique(Data.SetSize);
    assert(~any(isnan(uniqSizes)))
    assert(length(size(percepts)) == 2)
    d = nan(length(Data.SetSize), 1);
    
    for iU = 1 : length(uniqSizes)
        relTs = Data.SetSize == uniqSizes(iU);
        
        if strcmp(modelName, 'impSamp')
            d(relTs) = aisp_computeImpSampDV( ...
                percepts(relTs, :), ...
                uniqSizes(iU), ...
                relKappaX(relTs), ...
                Data.KappaS(relTs), ...
                0, ParamStruct.NumSamples, runChecks);
            
        elseif strcmp(modelName, 'jointPostSamp')
            d(relTs) = aisp_computeJointPostSampDV( ...
                percepts(relTs, :), ...
                uniqSizes(iU), ...
                relKappaX(relTs), ...
                Data.KappaS(relTs), ...
                0, ParamStruct.NumSamples, runChecks);
        end
    end
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