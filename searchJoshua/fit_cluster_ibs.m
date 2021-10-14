function fit_cluster_ibs(iRep, iPtpnt, type, DSet, idx, varargin)

% INPUT
% idx: Index of the job. Used when saving the results
% varargin{1}: Boolean. Default false. Whether to display information helpful
% for debuging.

if length(varargin) >= 1
    debugMode = varargin{1};
else
    debugMode = false;
end

% In debug model ibslike_var uses persistant variables. Clear these for a fresh
% start
if debugMode
    clear ibslike_var.m
end

addpath(genpath('../bads/'))
addpath(genpath('../ibs/'))
addpath('./visualSearch')
addpath('./visualSearch/analysisFuns')
addpath('../lautils-mat/stats')
addpath('../lautils-mat/math')

% Create the directory we will use for saving the results
if ~exist('pars', 'dir')
    mkdir('pars')
end

% Has this fit been run previously?
saveFile = sprintf('./pars/pars_%d_%s_%d_%d.mat',idx,type,iPtpnt,iRep);
if exist(saveFile, 'file')
    warning('Skipping this fit as it has already been completed.')
    return
end

DatSubj = DSet.P(iPtpnt).Data;
designMat = struct2DesignMat(DatSubj, 'to matrix');

% Fitting settings
opt_varLimit = 4; %4; can be higher if needed
opt_bads = bads;
opt_bads.NoiseFinalSamples = 100;
opt_bads.NoiseSize = sqrt(opt_varLimit);
opt_ibs = ibslike;
opt_ibs.Nreps = 1;
opt_ibs.MaxIter = 2 * (10^4);
% opt_ibs.Vectorized = 'on';

RespFun = @(pars, data) aisp_simResponseWrapper(data, pars, type);
fun_handle = @(pars) ibslike_var(RespFun, pars, DatSubj.Response, designMat, ...
    opt_ibs, opt_varLimit, debugMode);

[X0,LB,UB,PLB,PUB] = get_bads_bounds();

[pars,nLogL] = bads(fun_handle,X0,LB,UB,PLB,PUB,opt_bads);
save(saveFile,'pars','nLogL')


