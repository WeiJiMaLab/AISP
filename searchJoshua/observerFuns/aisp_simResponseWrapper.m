function response = aisp_simResponseWrapper(data, pars, modelName, ...
    varargin)
% Simulates responses for ibs-sampling, but accepts inputs in the form required
% for the IBS code.

% INPUT
% varargin{1}: bool. Run checks? If false, does not run some checks, 
%   therefore potentially saving time, but with less chance to detect 
%   bugs. Default is true.

if length(varargin) >= 1
    runChecks = varargin{1};
else
    runChecks = true;
end

DesignStruct = struct2DesignMat(data, 'to struct');
ParamStruct = paramVec2Struct(pars, modelName, 'to struct');

response = aisp_simResponse(modelName, ParamStruct, DesignStruct, ...
    runChecks);


