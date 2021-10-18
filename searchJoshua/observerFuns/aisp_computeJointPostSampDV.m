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

nAbsent = sum(proposeCat == 0, 3);
nPresent = sum(proposeCat == 1, 3);
allPropsl = nan(nTrials, size(percept, 2), nSamples);

for iTrial = 1 : nTrials 
    absentPropsl = drawSamples(percept(iTrial, :), kappa_s(iTrial), ...
                                nAbsent(iTrial, 1), 'targAbs');
    presentPropsl = drawSamples(percept(iTrial, :), kappa_s(iTrial), ...
                                nPresent(iTrial, 1), 'targPres');

    if runChecks
        assert(isequal(size(absentPropsl), [1, size(percept, 2), ...
            nAbsent(iTrial)]))
        assert(isequal(size(presentPropsl), [1, size(percept, 2), ...
            nPresent(iTrial)]))
    end
                            
    allPropsl(iTrial, :, proposeCat == 0) = absentPropsl;
    allPropsl(iTrial, :, proposeCat == 1) = presentPropsl;
end

if runChecks
   allNan = all(isnan(allPropsl), 3);
   noneNan = ~any(isnan(allPropsl), 3);
   allOrNone = allNan | noneNan;
   assert(all(allOrNone(:)))
end

phi = kappaTimesSumCosTerms(percept, allPropsl, kappa_x);
assert(isequal(shape(phi), [nTrials, 1, nSamples])) 

catSamples = nan(nTrials, 1, nSamples);
currentCat = nan(nTrials, 1);
currentPhi = nan(nTrial, 1);

initalCat = proposeCat(:, 1, 1);
catSamples(:, 1, 1) = initalCat;
currentCat(:) = initalCat;
currentPhi(:) = phi(:, 1, 1);

for iPropsl = 2 : nSamples
    pAccept = exp(phi(:, 1, iPropsl) - currentPhi);
    
    if runChecks
      assert(isequal(shape(pAccept), [nTrials, 1]))  
    end
        
    u = rand([nTrials, 1]);
    accepted = u < pAccept;
    
    currentCat(accepted) = proposeCat(accepted, 1, iPropsl);
    catSamples(:, 1, iPropsl) = currentCat;
    currentPhi(accepted) = phi(accepted, 1, iPropsl);
end

d = log( (1 + sum(catSamples==1, 3)) ./ (1 + sum(catSamples==0, 3)) );
assert(isequal(shape(d), [nTrials, 1]))  


