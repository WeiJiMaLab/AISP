function d = aisp_computeImpSampDV(percept, nItems, kappa_x, ...
    kappa_s, mu_s, nSamples, runChecks)
% Compute the decision variable of an observer who uses what we are calling
% importance sampling to evaluate tricky integrals 

% INPUT
% percept   [numTrials x setSize] array of stimulus percepts. In each 
%           row all nans should come after all non-nans. In other words
%           percept(:, 1:nItems) should all be not nans, and the remaining
%           columns should only contain nans.
% nItems    scalar. All trials passed must have the name number of 
%           stimuli items, the number given here.
% kappa_x   Observer's belief about the concetration parameter of the 
%           measurement noise
% kappa_s   Overserver's belief about the concentration parameter of the 
%           distractor distribution
% mu_s      Center of the distractor distribution
% nSamples  How many samples to take per trial
% runChecks bool. If true check potentially costly assertions

% JCT, 2021

dvFunsInputChecks(percept, mu_s, nItems, kappa_s, kappa_x)
assert(length(nItems) == 1)

% As specified in the description of the input, several columns in percept
% just contain nans. We can get rid of them to improve efficiency.
if runChecks
    toDiscard = percept(:, nItems+1:end);
end
percept = percept(:, 1:nItems);

if runChecks
    assert(~any(isnan(percept(:))))
    assert(all(isnan(toDiscard(:))))
end
% WORKING HERE
% Create [numTrials x setSize x numSamples] array
absentSamples = drawSamples(percept, kappa_s, nSamples, 'targAbs');
targetSamples = drawSamples(percept, kappa_s, nSamples, 'targPres');

d = evalOneTermInD(percept, targetSamples, kappa_x) - ...
        evalOneTermInD(percept, absentSamples, kappa_x);

assert(all(size(d) == [size(percept, 1), 1]))

end

function term = evalOneTermInD(percept, stimSamples, kappa_x)
% Evaluates one of the log-sum-exp terms in the expression derived for d
% in the derivations

exponent = kappaTimesSumCosTerms(percept, stimSamples, kappa_x);
term = logsumexp(exponent, 3);

assert(isequal(size(term), [size(percept, 1), 1]))

end







    
    