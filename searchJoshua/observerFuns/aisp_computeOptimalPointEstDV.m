function d = aisp_computeOptimalPointEstDV(percept, nItems, kappa_x, ...
    kappa_s, mu_s, runChecks)
% Compute the decision variable for the optimal point estimate observer

% INPUT
% nItems    [numTrials x 1] vectors. Different numbers of items in 
%           different trials is permitted.
% runChecks bool. If true check potentially costly assertions.

assert(length(size(percept)) == 2)

% First compute the standard point estiamte
d = aisp_computePointEstDV(percept, nItems, kappa_x, kappa_s, mu_s, ...
    'stimAndTarg', runChecks);

% Then calculate the optimal offset
uniqueCases = unique([nItems, kappa_x], 'rows');
if size(uniqueCases, 1) > 4
    error('Unepected number of cases given experiment code was written for.')
end
uniqueKappa_s = unique(kappa_s);

[ofset, assocNItems, assocKappa_x, assocKappa_s] = ...
    aisp_computeOptimalPointEstOfset(uniqueCases(:, 1), uniqueCases(:, 2), ...
    uniqueKappa_s, mu_s, runChecks);

% Apply the offset
allTrialsOfset = nan(size(nItems));

for iCase = 1 : length(ofset(:))
    thisNItems = assocNItems(iCase);
    thisKappa_x = assocKappa_x(iCase);
    thisKappa_s = assocKappa_s(iCase);
    
    relTrials = all( ...
        [thisNItems, thisKappa_x, thisKappa_s] == [nItems, kappa_x, kappa_s], 2);
    allTrialsOfset(relTrials) = ofset(iCase);
end

if any(isnan(allTrialsOfset(:))); error('Bug'); end

d = d - allTrialsOfset;