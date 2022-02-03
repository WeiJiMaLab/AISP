function [DSet, SimDSet] = aisp_simDataset(dataDir, parsDir, model)
% Takes the dataset given by DSet and simulates new responses using "model", and
% the best fitting parameters. Note that nothing except responses and accuracy 
% is simulated and SimDSet is in all other respects identical to DSet.

[DSet, ~] = getData(dataDir);
allParamStructs = aisp_loadBestFits(dataDir, parsDir, model);

%% Simulate data

SimDSet = aisp_simRespAndAcc(DSet, model, allParamStructs);

% TODO remove the duplication here
for iPtpnt = 1 : length(SimDSet.P)  
    SimDSet.P(iPtpnt).Data ...
        = computeStimStats(SimDSet.P(iPtpnt).Data, 'circ', false);
end

for iPtpnt = 1 : length(DSet.P)  
    DSet.P(iPtpnt).Data ...
        = computeStimStats(DSet.P(iPtpnt).Data, 'circ', false);
end






