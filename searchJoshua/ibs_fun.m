function response = ibs_fun(data, pars, type)
% simulates responses for ibs-sampling

try

%% Change format of data
DesignStruct = struct2DesignMat(data, 'to struct');
ParamStruct = paramVec2Struct(pars, 'to struct');


%% Simulate responses

% For all models we need to add noise to the real orientations to generate
% simulated percepts
kappaX = exp(ParamStruct.LnKappa_x);
relKappaX = kappaX(DesignStruct.SetSizeCond);
percepts = aisp_addNoiseToStim(relKappaX, DesignStruct.Orientation);

% Simulate responses based on these percepts
response = aisp_giveResponse(type, ParamStruct, DesignStruct, percepts);


catch
    disp('here')
end