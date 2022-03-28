function d = aisp_computeJointPostSampDV(percept, nItems, kappa_x, ...
    kappa_s, mu_s, nSamples, runChecks)
% Compute the decision variable of an observer who uses Metropolis Hastings
% sampling with imdependent draws from the prior for the proposal 
% distribution

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
percept = trimPercept(percept, nItems, runChecks);

% Because proposal samples are independent we can draw them all now
nTrials = size(percept, 1);
[proposeCat, allPropsl] = drawProposals(kappa_s, nTrials, ...
    nItems, nSamples, runChecks);

phi = kappaTimesSumCosTerms(percept, allPropsl, kappa_x);

if runChecks
    assert(isequal(find3Dsize(phi), [nTrials, 1, nSamples])) 
end

catSamples = nan(nTrials, 1, nSamples);
currentCat = nan(nTrials, 1);
currentPhi = nan(nTrials, 1);

initalCat = proposeCat(:, 1, 1);
catSamples(:, 1, 1) = initalCat;
currentCat(:) = initalCat;
currentPhi(:) = phi(:, 1, 1);

for iPropsl = 2 : nSamples
    pAccept = exp(phi(:, 1, iPropsl) - currentPhi);
    
    if runChecks
      assert(isequal(size(pAccept), [nTrials, 1]))  
    end
        
    u = rand([nTrials, 1]);
    accepted = u < pAccept;
    
    currentCat(accepted) = proposeCat(accepted, 1, iPropsl);
    catSamples(:, 1, iPropsl) = currentCat;
    currentPhi(accepted) = phi(accepted, 1, iPropsl);
end

d = log( (1 + sum(catSamples==1, 3)) ./ (1 + sum(catSamples==0, 3)) );

if runChecks
    assert(isequal(size(d), [nTrials, 1]))  
    sumBothCats = sum(catSamples==1, 3) + sum(catSamples==0, 3);
    assert(all(sumBothCats == nSamples))
end

end
