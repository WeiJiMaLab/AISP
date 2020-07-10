function response = ibs_fun(data, pars, type)
% simulates responses for ibs-sampling

%% Change format of data
DesignStruct = struct2DesignMat(data, 'to struct');
ParamStruct = paramVec2Struct(pars, 'to struct');


%% Simulate responses

% For all models we need to add noise to the real orientations to generate
% simulated percepts
relKappaX = ParamStruct.Kappa_x(DesignStruct.SetSizeCond);
percepts = aisp_addNoiseToStim(relKappaX, DesignStruct.Orientation);

% Simulate responses based on these percepts
response = aisp_giveResponse(type, ParamStruct, DesignStruct, percepts);