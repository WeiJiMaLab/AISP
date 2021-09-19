function d = aisp_computeSamplingDV(percept, nItems, kappa_x, kappa_s, ...
    mu_s)
% Compute the decision variable of an observer who uses sampling to
% evaluate tricky integrals 

% INPUT
% percept   [numTrials x setSize x numSamples] array of stimulus percepts
% nItems
% kappa_x   Observer's belief about the concetration parameter of the 
%           measurement noise
% kappa_s   Overserver's belief about the concentration parameter of the 
%           distractor distribution
% mu_s      Center of the distractor distribution

% JCT, 2021

% Check input
assert(mu_s == 0)

if size(percept, 2) > 8; error('Bug'); end
assert(length(size(percept)) == 2)

inputVectors = {nItems, kappa_s, kappa_x};

for iInputVec = 1 : length(inputVectors)
    vecSize = size(inputVectors{iInputVec});
    if (length(vecSize) ~= 2) || (vecSize(2) ~= 1)
        error('Bug')
    end
end

assert(length(nItems) == size(percept, 1))


% Process each possible number of items in turn
targetSamples = zeros(size(percept));
targetSamples(isnan(percept)) = nan;
targetSamples = addNoise(targetSamples, kappa_s);

%%% WOKRING HERE
% Need to also samples the absent samples as above, then samples target
% and set to zero

% Sample target location, and set these samples to 0
absentSamples = nan(size(percept));
for iNumItems = 1 : size(percept, 2)
    relTrials = nItems == iNumItems;
    
    if sum(relTrials) == 0
        continue
    end
    
    