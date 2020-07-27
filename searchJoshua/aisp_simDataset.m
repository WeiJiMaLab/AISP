function [SimDSet, DSet] = aisp_simDataset(dataDir, model)
% Takes the dataset given by DSet and simulates new responses using "model", and
% the best fitting parameters. Note that nothing except responses and accuracy 
% is simulated and SimDSet is in all other respects identical to DSet.

[DSet, Nptpnts] = getData(dataDir);
SimDSet = DSet;

%% Find the best fitting parameters

fname = ['pars/pars_' model '.mat'];
f = load(fname);

bestPars = aisp_collectBestFittingParams(f.nLogLs, f.pars);


%% Simulate data
for iP = 1 : Nptpnts
    PtpntData = SimDSet.P(iP).Data;
    ParamStruct = paramVec2Struct(bestPars(iP, :), 'to struct');
    
    resp = aisp_simResponse(model, ParamStruct, PtpntData);
    SimDSet.P(iP).Data.Response = resp;
    SimDSet.P(iP).Data.Accuracy = SimDSet.P(iP).Data.Response == SimDSet.P(iP).Data.Target;
end







