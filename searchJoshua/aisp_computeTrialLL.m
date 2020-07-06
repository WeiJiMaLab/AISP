function trialLL = aisp_computeTrialLL(Model, ParamStruct, Data, ~, varargin)
% Computes the loglikelihood for the trials passed to it. Complies with the
% specification required for the 'findMaximimLikelihood' function in the
% modellingTools repo. See README there, and description of 
% 'Settings.ComputTrialLL' for more info.

% INPUTS
% Model     This must be specified when set up the function handle in 
%           'Settings.ComputTrialLL'. See 'findDefaultModelSettings' for a
%           description of all model varients avaliable.
% ParamStruct
%           Structure containing a field for each fitted parameter.
% Data      Strcuture containing field for each feature of the
%           stimulus/behaviour. Fields contain arrays which are num trials long,
%           along the first dimention.
% varargin  Leave empty, or specify 'unitTest', to run giveResponseImp2 function
%           instead of giveResponse, so that outputs can be compared.

if isempty(varargin)
    unitTest = false;
elseif strcmp(varargin{1}, 'unitTest')
    unitTest = true;
else
    error('Incorrect use of inputs.')
end


% Simulation size. How many draws to simulate each trial?
nDraws = 1000; % Was 1000 or 5000

% Bootstrap von Mises samples to speed things up?
sampleShortcut = false;


% Duplicate data accordingly, along the third dimention
duplicatedOrientations = repmat(Data.Orientation, 1, 1, nDraws);
 
    
% Add noise to the stimuli to generate percepts
percepts = addNoiseToStim(Model, ParamStruct, Data, duplicatedOrientations, ...
    sampleShortcut);


% Make a response assuming that the trail is not the result of a lapse
if ~unitTest
    predictedResp = aisp_giveResponse(Model, ParamStruct, Data, percepts);
    
elseif unitTest
    error('Not implimented yet')
%     disp('********** UNIT TEST ************')
%     predictedResp = giveResponseImp2(Model, ParamStruct, Data, percepts);
    
end


% What are the likelihoods of the responses?
likelihoodYes = sum(predictedResp, 3)/size(predictedResp, 3);


% Due to using a finite number of draws, the estimated likelihood can be 0 or 1.
% Avoid taking the log of zero by constraining the likelihood to fall within 
% [(1/nDraws), 1 - (1/nDraws)]
likelihoodYes( likelihoodYes > (1 - (1/nDraws)) ) = 1 - (1/nDraws);
    
    
likelihoodYes( likelihoodYes < (1/nDraws) ) = 1/nDraws;
 

likelihoods = [1-likelihoodYes, likelihoodYes];

% TODO Remove lapses here
% Do we also need to account for lapses?
if strcmp(Model.Lapses, 'yes')
    
    error('Need to remove lapses here')
    
    likelihoods = ((1 - ParamStruct.LapseRate) * likelihoods) + ...
        (ParamStruct.LapseRate * repmat([0.5 0.5], size(likelihoods, 1), 1));
    
    
end


% What is the loglikelihood of the actual response
relevantLikelihood = sub2ind([size(likelihoods, 1), 2], ...
    [1 : size(likelihoods, 1)]', ...
    (Data.Response +1));


trialLL = log(likelihoods(relevantLikelihood));

    
end