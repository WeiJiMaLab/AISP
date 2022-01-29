function Data = aisp_simSingleCondStimulus(nTrials, nItems, setSizeCond, ...
    distStats, blockType)
% Simulate stimuli for a single condition

% INPUT
% nTrials   Number to simulate
% nItems    Number of items in the display in this condition 
% setSizeCond
%           All the different possible number of nItems in the experiment, should
%           be assigned a number running from 1 to the number of set size 
%           conditions. What is the number for the nItems in this simulation?
%           May be used, for example, to index into ParamStruct.Kappa_x
% distStats Struct of distractor statistics. Requred fields...
%   mu_s      Distractor mean
%   kappa_s   Distoractor von Mises concentration parameter
% blockType Each different block type in the experiment should be assigned a
%           number. What is the number for the current block type.

assert(isequal(size(nItems), [1, 1]))
assert(isequal(size(setSizeCond), [1, 1]))
assert(isequal(size(blockType), [1, 1]))
assert(isequal(size(distStats.mu_s), [1, 1]))
assert(isequal(size(distStats.kappa_s), [1, 1]))


% Set the probability of the target being present to 0.5
Data.Target = logical(randi([0 1], nTrials, 1));

% Randomise target locations
Data.TargetLoc = randi([1 nItems], nTrials, 1);
Data.TargetLoc(~Data.Target) = NaN;

% Randomise distractor orientations
Data.Orientation = circ_vmrnd_fixed(distStats.mu_s, distStats.kappa_s, ...
    [nTrials, nItems]);

% Set the orientations of the targets to zero.
% First need to find the indcies in Data.Orientation which correspond to targets
targetIndex = sub2ind([size(Data.Orientation)], ...
    find(Data.Target), Data.TargetLoc(Data.Target));
Data.Orientation(targetIndex) = 0;

% Pad the orientation data with NaNs until it is 6 columns wide
orientationDataSize = size(Data.Orientation);
Data.Orientation = ...
    [Data.Orientation, ...
    NaN(orientationDataSize(1), 6-orientationDataSize(2))];

numTargs = sum(Data.Orientation == 0, 2);
assert(all(numTargs(Data.Target)==1))
assert(all(numTargs(~Data.Target)==0))

% Add some properties of the stimulus to the data struct
Data.SetSize = repmat(nItems, nTrials, 1);
Data.BlockType = repmat(blockType, nTrials, 1);
Data.SetSizeCond = repmat(setSizeCond, nTrials, 1);
Data.KappaS = repmat(distStats.kappa_s, nTrials, 1);




