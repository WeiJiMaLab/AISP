function [proposeCat, allPropsl] = drawProposals(kappa_s, nTrials, ...
    nItems, nPropsls, runChecks)
% Draw independent proposals for the MCMC smapler

% INPUT
% kappa_s: vector as long as number of trials. Overserver's belief about 
%   the concentration parameter of the distractor distribution
% nItems: scalar. All trials passed must have the name number of 
%           stimuli items, the number given here.
% nPropsls: How many proposals to generate for each trial

% OUTPUT
% proposeCat: [nTrials, 1, nSamples] array, describing the cateogry of 
%  each proposal
% allPropsl: [nTrials, nItems, nSamples] array of stimulus orientations
%   for each proposal

% JCT, 2021

assert(length(kappa_s) == nTrials)
assert(length(nItems) == 1)

proposeCat = nan(nTrials, 1, nPropsls);
allPropsl = nan(nTrials, nItems, nPropsls);

uniKappas = unique(kappa_s);
for iKappa = 1 : length(uniKappas)
    thisKappa = uniKappas(iKappa);
    match = thisKappa == kappa_s;
    nMatch = sum(match);
    
    [theseProposeCat, thesePropsl] = drawThisKappaPropsls(nMatch, ...
        nItems, nPropsls, thisKappa, runChecks);
    
    proposeCat(match, 1, :) = theseProposeCat;
    allPropsl(match, :, :) = thesePropsl;
end

if runChecks
    assert(~any(isnan(proposeCat(:))))
    assert(~any(isnan(allPropsl(:))))
    
    numTargets = sum(allPropsl == 0, 2);
    assert(isequal(size(numTargets), size(proposeCat)))
    
    numTargetsUnderPres = numTargets(proposeCat == 1);
    assert(all(numTargetsUnderPres(:) == 1))
    
    numTargetsUnderAbs = numTargets(proposeCat == 0);
    assert(all(numTargetsUnderAbs(:) == 0))
end

end

function [proposeCat, allPropsl] = drawThisKappaPropsls(nTrials, ...
    nItems, nPropsls, kappa_s, runChecks)
% Like draw proposals but just works for a single value of kappa_s

assert(length(kappa_s) == 1)

proposeCat = randi(2, [nTrials*nPropsls, 1]) -1;
unshapedPropsls = nan(nTrials*nPropsls, nItems);

numAbsent = sum(proposeCat == 0);
thisKappa_s = repmat(kappa_s, [numAbsent, 1]); % TODO this could be made 
% more efficient
absentSamples = drawSamples(thisKappa_s, numAbsent, nItems, ...
    1, 'targAbs', runChecks);
unshapedPropsls(proposeCat == 0, :) = absentSamples;

numPresent = sum(proposeCat == 1);
thisKappa_s = repmat(kappa_s, [numPresent, 1]); % TODO this could be made 
% more efficient
targetSamples = drawSamples(thisKappa_s, numPresent, nItems, ...
    1, 'targPres', runChecks);
unshapedPropsls(proposeCat == 1, :) = targetSamples;

if runChecks
    assert(~any(isnan(unshapedPropsls(:))))
end

% Reshape to [nTrials x 1 x nPropsls] or [nTrials x nItems x nPropsls]
proposeCat = reshape(proposeCat, [nTrials, 1, nPropsls]);
allPropsl = nan(nTrials, nItems, nPropsls);
for iPropsl = 1 : nPropsls
    startRow = 1 + ((iPropsl -1)*nTrials);
    endRow = iPropsl * nTrials;
    
    if runChecks
        nanInPropsl = isnan(allPropsl(:, :, iPropsl));
        assert(all(nanInPropsl(:)))
    end
    
    allPropsl(:, :, iPropsl) = unshapedPropsls(startRow:endRow, :);
end
assert(endRow == size(unshapedPropsls, 1))

end