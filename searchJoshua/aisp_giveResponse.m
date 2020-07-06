function resp = aisp_giveResponse(Model, ParamStruct, Data, percepts)
% Make a response assuming that the trial is not the result of a lapse

% INPUT
% stim   [num trials X num locations] array describing the perceived orientation 
%        of the Gabor patches


% What response is given in each case?
if strcmp(Model, 'bayes')
    relKappaX = ParamStruct.Kappa_x(Data.SetSizeCond);
    d = makeBaysianDecision(percepts, Data.SetSize, relKappaX, ...
        Data.KappaS, ParamStruct);

elseif strcmp(Model, 'pointEst')
    
%     % Are we ignoring the set size when applying thresholds?
%     if strcmp(Model.SetSizeThresh, 'variable')
%        relSetSizes = Data.SetSizeCond; 
%         
%     elseif strcmp(Model.SetSizeThresh, 'fixed')
%         relSetSizes = ones(size(Data.SetSizeCond));
% 
%     end
%     
%     % Are we ignoring the block type when applying thresholds?
%     if strcmp(Model.BlockTypes, 'use')
%         relBlockTypes = Data.BlockType; 
%         
%     elseif strcmp(Model.BlockTypes, 'ignore')
%         relBlockTypes = ones(size(Data.SetSizeCond));
%         
%     end
%     
%     resp = makeMinRuleDecision(ParamStruct.Thresh, percepts, relSetSizes, ...
%             relBlockTypes, 0);
        
elseif strcmp(Model, 'optimalPointEst')
    
end


% Compute the probability of a target present response 
pPresent = (ParamStruct.LapseRate / 2) + ((1 - ParamStruct.LapseRate) * ...
    (1 ./ (1 + exp(-ParamStrcut.Beta0 - (ParamStrcut.Beta1 * d)))));

% Simulate a response
if size(pPresent, 2) ~=1; error('Bug'); end
resp = rand(size(pPresent)) < pPresent;