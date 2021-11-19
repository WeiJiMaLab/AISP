function samples = drawSamples(kappa_s, nTrials, setSize, nSamples, ...
    targPresent, runChecks)
% Draw samples for the stimulus from p(s | C)

% INPUT
% kappa_s: vector as long as nTrials, and determines the kappa for drawing 
%   orientations
% nTrials: scalar.
% setSize: scalar. Will draw nSamples for nTrials trials in which this
%   many stimuli were presented each trial.
% nSamples: scalar. How many samples to draw?
% targPresent: str. If 'targPres' draw from p(s | C ='target present'), 
%   else if 'targAbs' draw from p(s | C ='target absent')
% runChecks bool. If true check potentially costly assertions

% OUTPUT
% samples: [numTrials x setSize x numSamples] array

% JCT, 2021

if nSamples == 0
    error('No samples requested')
end
assert(size(kappa_s, 2)==1)
assert(size(kappa_s, 1)== nTrials)

samples = drawAbsentSamples(kappa_s, nTrials, setSize, nSamples, runChecks);

if strcmp(targPresent, 'targPres')
    % Just need to add target locations to samples drawn as if no 
    % target present
    targLocs = randi(setSize, [nTrials, 1, nSamples]);
    correspondTrial = nan(nTrials, 1);
    correspondTrial(:) = 1:nTrials;
    correspondTrial = repmat(correspondTrial, [1, 1, nSamples]);
    correspondSample = nan(1, 1, nSamples);
    correspondSample(:) = 1:nSamples;
    correspondSample = repmat(correspondSample, [nTrials, 1, 1]);
    targIdx = sub2ind([nTrials, setSize, nSamples], ...
        correspondTrial(:), targLocs(:), correspondSample(:));
    
    samples(targIdx) = 0;
    
    if runChecks
        numTargs = sum(samples == 0, 2);
        assert(isequal(find3Dsize(numTargs), [nTrials, 1, nSamples]))
        assert(all(numTargs(:) == 1))
    end
    
elseif ~strcmp(targPresent, 'targAbs')
    error('Incorrect use of inputs')
end

if runChecks
    assert(isequal(find3Dsize(samples), [nTrials, setSize, nSamples]))
end

end

function absentSamples = drawAbsentSamples(kappa_s, nTrials, setSize, ...
    nSamples, runChecks)
% Draw samples assuming target is absent. 

absentSamples = nan([nTrials, setSize, nSamples]);

uniKappas = unique(kappa_s);
for iKappa = 1 : length(uniKappas)
    thisKappa = uniKappas(iKappa);
    match = thisKappa == kappa_s;
    numReq = sum(match)*setSize*nSamples;
    
    unshaped = sampVm(0, thisKappa, numReq);
    absentSamples(match, :, :) = reshape(unshaped, ...
        [sum(match), setSize, nSamples]);
end

if runChecks
    assert(isequal(find3Dsize(absentSamples), ...
        [nTrials, setSize, nSamples]))
    assert(~any(isnan(absentSamples(:))))
end

end


