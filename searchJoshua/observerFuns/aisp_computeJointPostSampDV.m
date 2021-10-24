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
    allPropsl = addProposals(allPropsl, proposeCat, percept, kappa_s, ...
        'absent', nAbsent, nPresent, iTrial, runChecks);
    allPropsl = addProposals(allPropsl, proposeCat, percept, kappa_s, ...
        'present', nAbsent, nPresent, iTrial, runChecks);
end

if runChecks
   allNan = all(isnan(allPropsl), 3);
   noneNan = ~any(isnan(allPropsl), 3);
   allOrNone = allNan | noneNan;
   assert(all(allOrNone(:)))
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

function allPropsl = addProposals(allPropsl, proposeCat, percept, ...
    kappa_s, targCase, nAbsent, nPresent, iTrial, runChecks)
% Add proposals to the matrix containing all proposals, allPropsl. Either
% add proposals for the case where the target is present, or proposals 
% for the case where the target is absent, depending on targCase.

% INPUT
% targCase: str. 'present' or 'absent'. Which case to add proposals for.
% nAbsent, nPresent: vectors. Long as the number of trials. Describing 
% how many target absent, and target present proposals to draw for each
% trial
% iTrial: Which trial to add proposals for.

assert(size(proposeCat, 2) == 1)

if strcmp(targCase, 'present')
    nSamples = nPresent(iTrial, 1);
    sampleType = 'targPres';
    relProposals = proposeCat(iTrial, 1, :) == 1;
elseif strcmp(targCase, 'absent')
    nSamples = nAbsent(iTrial, 1);
    sampleType = 'targAbs';
    relProposals = proposeCat(iTrial, 1, :) == 0;
else
    error('Bug')
end

if nSamples == 0
    % Take a break Matlab
else
    thisPropsl = drawSamples(percept(iTrial, :), kappa_s(iTrial), ...
        nSamples, sampleType);
    
    if runChecks
        assert(isequal(find3Dsize(thisPropsl), ...
            [1, size(percept, 2), nSamples]))
    end
    
    allPropsl(iTrial, :, relProposals) = thisPropsl;
end

end
