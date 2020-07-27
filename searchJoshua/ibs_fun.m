function response = ibs_fun(data, pars, type)
% simulates responses for ibs-sampling

%% Change format of data
DesignStruct = struct2DesignMat(data, 'to struct');
ParamStruct = paramVec2Struct(pars, 'to struct');


%% Simulate responses
response = aisp_simResponse(type, ParamStruct, DesignStruct);


