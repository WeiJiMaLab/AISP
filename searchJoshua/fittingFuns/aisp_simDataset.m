function [DSet, SimDSet] = aisp_simDataset(dataDir, parsDir, model)
% Takes the dataset given by DSet and simulates new responses using "model", and
% the best fitting parameters. Note that nothing except responses and accuracy 
% is simulated and SimDSet is in all other respects identical to DSet.

[DSet, Nptpnts] = getData(dataDir);
SimDSet = DSet;

%% Find the best fitting parameters

fname = [parsDir '/pars_' model '.mat'];
f = load(fname);

bestPars = aisp_collectBestFittingParams(f.nLogLs, f.pars);


%% Simulate data
for iP = 1 : Nptpnts
    PtpntData = SimDSet.P(iP).Data;
    error('Inputs to paramVec2Struct changed')
    ParamStruct = paramVec2Struct(bestPars(iP, :), 'to struct');
    
    resp = aisp_simResponse(model, ParamStruct, PtpntData);
    SimDSet.P(iP).Data.Response = resp;
    SimDSet.P(iP).Data.Accuracy = SimDSet.P(iP).Data.Response == SimDSet.P(iP).Data.Target;
end

%% Add summary statitics 

% TODO remove the duplication here
for iPtpnt = 1 : length(SimDSet.P)  
    SimDSet.P(iPtpnt).Data ...
        = computeStimStats(SimDSet.P(iPtpnt).Data, 'circ', false);
end

for iPtpnt = 1 : length(DSet.P)  
    DSet.P(iPtpnt).Data ...
        = computeStimStats(DSet.P(iPtpnt).Data, 'circ', false);
end






