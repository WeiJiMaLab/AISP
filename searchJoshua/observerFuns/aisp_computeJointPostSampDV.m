function d = aisp_computeJointPostSampDV(percept, nItems, kappa_x, ...
    kappa_s, mu_s, nSamples, runChecks)
% Compute the decision variable of an observer who uses Metropolis Hastings
% sampling with imdependent draws from the prior for the proposal 
% distribution

% INPUT
% percept   [numTrials x setSize] array of stimulus percepts
% nItems
% kappa_x   Observer's belief about the concetration parameter of the 
%           measurement noise
% kappa_s   Overserver's belief about the concentration parameter of the 
%           distractor distribution
% mu_s      Center of the distractor distribution
% nSamples  How many samples to take per trial
% runChecks bool. If true check potentially costly assertions

% JCT, 2021

dvFunsInputChecks(percept, mu_s, nItems, kappa_s, kappa_x)

% Because proposal samples are independent can draw them all now
nTrials = size(percept, 1);
proposeCat = randi(2, [nTrials, 1, nSamples]) -1;

absentSamples = drawSamples(percept, kappa_s, nSamples, 'targAbs');
targetSamples = drawSamples(percept, kappa_s, nSamples, 'targPres');

allPropsl = nan(nTrials, size(percept, 2), nSamples);

for iPropsl = 1 : nSamples
    presTrials = proposeCat(:, 1, iPropsl) == 1;
    allPropsl(presTrials, :, iPropsl) = ...
        targetSamples(presTrials, :, iPropsl);
    
    absTrials = proposeCat(:, 1, iPropsl) == 0;
    allPropsl(absTrials, :, iPropsl) = ...
        absentSamples(absTrials, :, iPropsl);
end

if runChecks
   allNan = all(isnan(allPropsl), 3);
   noneNan = ~any(isnan(allPropsl), 3);
   allOrNone = allNan | noneNan;
   assert(all(allOrNone(:)))
   
   assert(isequal(allNan, isnan(percept)))
end

phi = kappaTimesSumCosTerms(percept, allPropsl, kappa_x);
assert(isequal(find3Dsize(phi), [nTrials, 1, nSamples])) 

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
assert(isequal(size(d), [nTrials, 1]))  


end
