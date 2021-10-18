function response = aisp_simResponseWrapper(data, pars, type)
% Simulates responses for ibs-sampling, but accepts inputs in the form required
% for the IBS code.

DesignStruct = struct2DesignMat(data, 'to struct');
ParamStruct = paramVec2Struct(pars, 'to struct');

response = aisp_simResponse(type, ParamStruct, DesignStruct);


