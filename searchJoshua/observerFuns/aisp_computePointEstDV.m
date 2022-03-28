function d = aisp_computePointEstDV(percept, nItems, kappa_x, kappa_s, mu_s, ...
    maxOver, runChecks)
% Compute the decision variable for the point estimate observer

% INPUT
% percept   [numTrials x setSize x numSamples] array of stimulus percepts
% nItems    [numTrials x 1] vectors. Different numbers of items in 
%           different trials is permitted.
% kappa_x   Observer's belief about the concetration parameter of the 
%           measurement noise
% kappa_s   Overserver's belief about the concentration parameter of the 
%           distractor distribution
% mu_s      Center of the distractor distribution
% maxOver   What kind of point estimate observer should the function behave
%           like? An observer who maximises over both stimuli and the
%           location of the target ('stimAndTarg'), or and obeserver who
%           maximises over both stimulus and location even when the target
%           is not present and hence 'location' doesn't really exist 
%           ('stimAndTarg_incImagine')?
% runChecks bool. If true check potentially costly assertions.

%% Setup
% Check input
if size(percept, 2) > 8; error('Bug'); end
assert(length(size(percept)) == 2)
if ~ismember(maxOver, {'stimAndTarg', 'stimAndTarg_incImagine'})
    error('Incorrect use of inputs.')
end

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


%% Main calculations

if ~(length(mu_s) == 1 && mu_s == 0)
    error('Not coded up yet')
end

cosPercept = cos(percept);
kappa_d = sqrt((kappa_x.^2) + (kappa_s.^2) + (2*kappa_x.*kappa_s.*cosPercept));

% Find the term that needs to be maximised over
toMax = (kappa_x.*cosPercept) - kappa_d;

if runChecks
    assert(isequal(size(kappa_d), size(percept)))
    assert(isequal(size(toMax), size(percept)))
end

% If the observer uses an "imaginary" variable for L even when the
% target is absent then we must also include an additional term
% featuring N
if strcmp(maxOver, 'stimAndTarg_incImagine')
    additionalTerm = log(nItems);
else
    additionalTerm = 0;
end

d = log((2*pi) ./ nItems) ...
    + additionalTerm ...
    + aisp_computeLogBesseliForDuplicatedValues(kappa_s) ...
    + max(toMax, [], 2);

if runChecks
    assert(isequal(size(d), [size(percept, 1), 1]))
end


