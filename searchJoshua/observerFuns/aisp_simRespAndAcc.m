function SimDSet = aisp_simRespAndAcc(DSet, model, paramStructs)
% Takes the dataset given by DSet and simulates new responses using 
% "model", and the best fitting parameters as described by the parameter
% structures. Note that nothing except responses and accuracy 
% is simulated and SimDSet is in all other respects identical to DSet.

% INPUT
% ParamStructs: cell array. As long as the number of participants. Each 
%   element is a paramter stucture describing the parameters to use for 
%   the correspoding participant.

SimDSet = DSet;

for iP = 1 : length(DSet.P)
    PtpntData = SimDSet.P(iP).Data;
    ParamStruct = paramStructs{iP};
    
    % Send all the nan's in PtpntData.Orientation to the end of each row
    % as expected by the sim response functions
    designMat = struct2DesignMat(PtpntData, 'to matrix', true);
    PtpntData = struct2DesignMat(designMat, 'to struct');
    
    resp = aisp_simResponse(model, ParamStruct, PtpntData);
    SimDSet.P(iP).Data.Response = resp;
    SimDSet.P(iP).Data.Accuracy = SimDSet.P(iP).Data.Response == ...
        SimDSet.P(iP).Data.Target;
end
