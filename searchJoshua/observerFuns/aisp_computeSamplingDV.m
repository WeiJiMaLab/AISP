function d = aisp_computeSamplingDV(percept, nItems, kappa_x, kappa_s, ...
    mu_s, nSamples)
% Compute the decision variable of an observer who uses sampling to
% evaluate tricky integrals 

% INPUT
% percept   [numTrials x setSize] array of stimulus percepts
% nItems
% kappa_x   Observer's belief about the concetration parameter of the 
%           measurement noise
% kappa_s   Overserver's belief about the concentration parameter of the 
%           distractor distribution
% mu_s      Center of the distractor distribution
% nSamples  How many samples to take per trial

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

% Draw samples assuming target is absent. Will use a 
% [numTrials x setSize x numSamples] array
absentSamples = drawAbsentSamples(percept, kappa_s, nSamples);

% Draw samples assuming target is present. Can just draw more samples
% as if absent, and then add targets
targetSamples = drawAbsentSamples(percept, kappa_s, nSamples);

% Simulating target locations is a little tricky as not all positions in
% percept are occupied. Some are nan's representing that no stimulus was
% presented there at all. Create an array representing the
% probability that location at index i or lower is selected as the target. 
activeLocs = ~isnan(percept);
cumulatProb = cumsum(activeLocs, 2);
cumulatProb = cumulatProb ./ cumulatProb(:, end);

% Draw a value between 0 and 1, and then use to select a location
% acording to the above cumulative probabilities
uDraw = rand([size(targetSamples, 1), 1, size(targetSamples, 3)]);

lowerCut = [zeros(size(cumulatProb, 1), 1), cumulatProb(:, 2:end)];
upperCut = cumulatProb;

targMask = (uDraw > lowerCut) & (uDraw < upperCut);
assert(isequal(size(targMask), size(targetSamples)))
numTargs = sum(targMask, 2);
assert(all(numTargs(:) == 1))

targetSamples(targMask) = 0;
assert(isequal(isnan(targetSamples), isnan(percept)))

% Evaluate the samples
% WORKING HERE AND WORKING IN evalOneTermInD. Just finished but not checked

assert(all(size(d) == [size(percept, 1), 1]))

end

function absentSamples = drawAbsentSamples(percept, kappas, nSamples)
% Draw samples assuming target is absent. 

% INPUT
% percept: [numTrials x setSize] array of stimulus percepts. Only used 
%   to determine how many samples to draw, and which positions should 
%   be set to nan
% kappas: vector as long as the first dimension of percept, and 
%   determines the kappa for drawing orientations
% nSamples: scalar. How many samples to draw?

% OUTPUT
% samples: [numTrials x setSize x numSamples] array

assert(length(size(percept)) == 2)

absentSamples = zeros(size(percept));
absentSamples(isnan(percept)) = nan;

currentNDims = length(size(absentSamples));
assert(currentNDims == 2);
absentSamples = repmat(absentSamples, [ones(1, currentNDims), nSamples]);

absentSamples = addNoise(absentSamples, kappa_s);

end


function term = evalOneTermInD(percept, stimSamples, kappa_x)
% Evaluates one of the log-sum-exp terms in the expression derived for d
% in the derivations

summed = sum(cos(percept - stimSamples), 2, 'omitnan');
exponent = kappa_x .* summed;
term = logsumexp(exponent, 3);

assert(isequal(shape(term), [size(percept, 1), 1]))

end







    
    