function d = aisp_computeBaysianDV(percept, nItems, kappa_x, kappa_s, ...
    mu_s, runChecks)
% Compute the decision variable of a Bayesian observer

% NOTE
% I have only focussed on making the code efficient for when mu_s==0. Lots
% of the same tricks could be used for when this is not the case but they
% have not been implimented.

% INPUT
% percept   [numTrials x setSize] array of stimulus percepts
% nItems    [numTrials x 1] vectors. Different numbers of items in 
%           different trials is permitted.
% kappa_x   Observer's belief about the concetration parameter of the 
%           measurement noise
% kappa_s   Overserver's belief about the concentration parameter of the 
%           distractor distribution
% mu_s      Center of the distractor distribution
% runChecks bool. If true check potentially costly assertions.

% JCT, 2020

% Check input
if size(percept, 2) > 8; error('Bug'); end
assert(length(size(percept)) == 2)

% We use implicit expansion below so it is very important that all input
% vectors are the expected shape.
inputVectors = {nItems, kappa_s, kappa_x};

for iInputVec = 1 : length(inputVectors)
    vecSize = size(inputVectors{iInputVec});
    if (length(vecSize) ~= 2) || (vecSize(2) ~= 1)
        error('Bug')
    end
end

% If kappa_s is passed as a single value, expand it into a vector so that
% the calculations for 'termC' below, work.
if size(kappa_s, 1) == 1
    kappa_s = repmat(kappa_s, size(percept, 1), 1);
end

% Same for kappa_x
if size(kappa_x, 1) == 1
    kappa_x = repmat(kappa_x, size(percept, 1), 1); 
end

% Do some computations to speed up the loglikelihood calculation later

% We need to compute log(besseli(0, kappa_s) but there are only a few values of
% kappa_s so we can speed the computation up by not doing the same thing over
% and over again. Lets call it termA.
termA = aisp_computeLogBesseliForDuplicatedValues(kappa_s);

if runChecks
    assert(isequal(size(termA), size(kappa_s)))
end

% If mu_s aways equals zero, then the formula for the local loglikelihood 
% ratio requires us to compute cos(percept) twice. Lets just do it once
if length(mu_s) == 1 && mu_s == 0
    termB = cos(percept);
    termBShortcut = true; 
else
    error('Case not coded up.')
end

% Compute local loglikelkihood ratio. The shortcuts we use (described
% above) depend on the situation. Note instead of using repmat to make
% everything the same size, we rely on MATLABs explicit expansion. (Will
% only work on MATLAB R2016b and later.)
if termBShortcut
    
    % If kappa_s == 0, then there is a whole term we do not need to calculate.
    % Only calculate it when kappa_s ~= 0
    termC = NaN(size(percept));

    % Then change this value to the correct one wherever kappa_s~=0 
    calcTrials = kappa_s ~= 0;
    
    % Set the termC for all activeLocs to the value at kappa_s==0
    trialTermC = logBesseliWithNoOverflow(kappa_x(~calcTrials));
    
    termC(~calcTrials, :, :) = repmat(trialTermC, ...
        1, size(termC, 2), size(termC, 3));
    
    if sum(calcTrials) > 0
        arg = ( (kappa_x(calcTrials).^2) + ...
            (kappa_s(calcTrials).^2) + ...
            (2*kappa_x(calcTrials).*kappa_s(calcTrials).* ...
            termB(calcTrials, :, :)) ).^0.5;
        termC(calcTrials, :, :) = logBesseliWithNoOverflow(arg);    
    end
    
    d_loc = (kappa_x .* termB) + termA - termC;
    
    if runChecks
        assert(isequal(size(d_loc), size(percept)))
    end
end

% Compute overal loglikelihood ratio that target is present vs. absent, d
d = log( (1./nItems) .* nansum(exp(d_loc), 2) );

if runChecks
    assert(isequal(size(d), [size(percept, 1), 1]))
end

end


function result = logBesseliWithNoOverflow(z)
% Compute log(besseli(0, z)) even for very large values of z using the 
% scaling option on besseli provided by matlab

result = log(besseli(0, z, 1)) + abs(real(z));


end




