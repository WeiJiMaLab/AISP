function response = aisp_simResponseWrapper(data, pars, type)
% Simulates responses for ibs-sampling, but accepts inputs in the form required
% for the IBS code.

%% Change format of data
DesignStruct = struct2DesignMat(data, 'to struct');
ParamStruct = paramVec2Struct(pars, 'to struct');


%% Simulate responses
response = aisp_simResponse(type, ParamStruct, DesignStruct);


