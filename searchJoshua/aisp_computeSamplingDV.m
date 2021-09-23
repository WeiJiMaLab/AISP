function d = aisp_computeSamplingDV(percept, nItems, kappa_x, kappa_s, ...
    mu_s)
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

% Draw samples assuming target is absent
absentSamples = zeros(size(percept));
absentSamples(isnan(percept)) = nan;
absentSamples = addNoise(absentSamples, kappa_s);

% Draw samples assuming target is present
targetSamples = zeros(size(percept));
targetSamples(isnan(percept)) = nan;
targetSamples = addNoise(targetSamples, kappa_s);

% Simulating target locations is a little tricky as not all positions in
% percept are occupied. Some are nan's representing that no stimulus was
% presented there at all
warning('Check code looks good in operation') % WORKING HERE -- do this check
targetLoc = false(size(percept));
toDo = true(size(targetLoc, 1), 1);
while sum(toDo) > 0
    proposeTargLoc = randi(size(targetLoc, 2), sum(toDo), 1);
    linIdx = sub2ind(size(targetLoc), find(toDo), proposeTargLoc);
    targetLoc(linIdx) = true;
    
    targetLoc(isnan(targetSamples)) = false;
    
    numTargets = sum(targetLoc, 2);
    toDo = numTargets ~= 1;
    assert(all(numTargets(toDo) == 0))
end

targetSamples(targetLoc) = 0;
assert(isequal(isnan(targetSamples), isnan(percept)))

% Evaluate the samples
vals = % WORKING HERE -- need to impliment a check that a row isn't just nans because
% in this case prod() returns 1
absIntegral = prod(vmpdf(percept, absentSamples, kappa_x), 2, 'omitnan');
targetIntegral = nanprod(vmpdf(percept, targetSamples, kappa_x), 2);

d = log(targetIntegral / absIntegral);
assert(all(size(d) == [size(percept, 1), 1]))










    
    