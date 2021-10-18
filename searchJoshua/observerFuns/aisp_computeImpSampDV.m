function d = aisp_computeImpSampDV(percept, nItems, kappa_x, ...
    kappa_s, mu_s, nSamples)
% Compute the decision variable of an observer who uses what we are calling
% importance sampling to evaluate tricky integrals 

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

dvFunsInputChecks(percept, mu_s, nItems, kappa_s, kappa_x)

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







    
    