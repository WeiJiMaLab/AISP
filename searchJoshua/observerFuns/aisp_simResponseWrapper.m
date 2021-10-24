function response = aisp_simResponseWrapper(data, pars, modelName)
% Simulates responses for ibs-sampling, but accepts inputs in the form required
% for the IBS code.

DesignStruct = struct2DesignMat(data, 'to struct');
ParamStruct = paramVec2Struct(pars, modelName, 'to struct');

response = aisp_simResponse(modelName, ParamStruct, DesignStruct);


