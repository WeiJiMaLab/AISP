function percept = trimPercept(percept, nItems, runChecks)
% Removes columns of nans from percept (assuming percept is provided in
% the fornat described below)

% INPUT
% percept   [numTrials x setSize] array of stimulus percepts. In each 
%           row all nans should come after all non-nans. In other words
%           percept(:, 1:nItems) should all be not nans, and the remaining
%           columns should only contain nans.
% nItems    scalar. All trials passed must have the same number of 
%           stimuli items, the number given here.

% JCT, 2021

if runChecks
    toDiscard = percept(:, nItems+1:end);
end

percept = percept(:, 1:nItems);

if runChecks
    assert(~any(isnan(percept(:))))
    assert(all(isnan(toDiscard(:))))
end