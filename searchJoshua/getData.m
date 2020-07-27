function [DSet, Nptpnts] = getData(dataDir)

if strcmp(dataDir, 'local')
    dataDir = 'D:\Research data _ BACKED UP\Visual search\Main Study\In use\StandardFormat_participantExcluded.mat';
end
Loaded = load(dataDir);
DSet = Loaded.DSet;

% Throughout we assume that the target is at the orientation of 0. Check this.
% TODO in some places it is not made clear that this is built into the code and
% not simply an option -- change
for iP = 1 : length(DSet.P)
    targetTrials = logical(DSet.P(iP).Data.Target);
    targetLoc = DSet.P(iP).Data.TargetLoc(targetTrials);
    allOrientations = DSet.P(iP).Data.Orientation(targetTrials, :);
    targetIdx = sub2ind(size(allOrientations), [1:sum(targetTrials)]', targetLoc);
    targetVals = allOrientations(targetIdx);
    
    if any(targetVals(:) ~= 0); error('Bug'); end
end

Nptpnts = length(DSet.P);