function samples = drawSamples(percept, kappa_s, nSamples, targPresent)
% Draw samples for the stimulus from p(s | C)

% INPUT
% percept: [numTrials x setSize] array of stimulus percepts. Only used 
%   to determine how many trials to draw samples for, and which positions 
%   should be set to nan
% kappa_s: vector as long as the first dimension of percept, and 
%   determines the kappa for drawing orientations
% nSamples: scalar. How many samples to draw?
% targPresent: str. If 'targPres' draw from p(s | C ='target present'), 
%   else if 'targAbs' draw from p(s | C ='target absent')

% OUTPUT
% samples: [numTrials x setSize x numSamples] array

% JCT, 2021

if nSamples == 0
    error('No samples requested')
end

samples = drawAbsentSamples(percept, kappa_s, nSamples);

if strcmp(targPresent, 'targPres')
    % Just need to add target locations to samples drawn as if no 
    % target present

    % Simulating target locations is a little tricky as not all positions in
    % percept are occupied. Some are nan's representing that no stimulus was
    % presented there at all. Create an array representing the
    % probability that location at index i or lower is selected as the target. 
    activeLocs = ~isnan(percept);
    cumulatProb = cumsum(activeLocs, 2);
    cumulatProb = cumulatProb ./ cumulatProb(:, end);

    % Draw a value between 0 and 1, and then use to select a location
    % acording to the above cumulative probabilities
    uDraw = rand([size(samples, 1), 1, size(samples, 3)]);

    lowerCut = [zeros(size(cumulatProb, 1), 1), cumulatProb(:, 1:end-1)];
    upperCut = cumulatProb;

    targMask = (uDraw > lowerCut) & (uDraw < upperCut);
    assert(isequal(size(targMask), size(samples)))
    numTargs = sum(targMask, 2);
    if ~all(numTargs(:) == 1)
       error('Bug') 
    end

    samples(targMask) = 0;
    
    samplesSubsec = samples(:, :, randi(size(samples, 3), [1, 1]));
    if ~isequal(isnan(samplesSubsec), isnan(percept))
        error('Bug') 
    end
    
elseif ~strcmp(targPresent, 'targAbs')
    error('Incorrect use of inputs')
end

assert(isequal(find3Dsize(samples), [size(percept), nSamples]))

end

function absentSamples = drawAbsentSamples(percept, kappa_s, nSamples)
% Draw samples assuming target is absent. 

assert(length(size(percept)) == 2)

absentSamples = zeros(size(percept));
absentSamples(isnan(percept)) = nan;

currentNDims = length(size(absentSamples));
assert(currentNDims == 2);
absentSamples = repmat(absentSamples, [ones(1, currentNDims), nSamples]);

absentSamples = addNoise(absentSamples, kappa_s, 'efficientSamp');
assert(isequal(find3Dsize(absentSamples), [size(percept), nSamples]))

end


