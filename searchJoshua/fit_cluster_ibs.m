function fit_cluster_ibs(iRep, iPtpnt, type, DSet)

addpath(genpath('../bads/'))
addpath(genpath('../ibs/'))
addpath('./visualSearch')
addpath('../lautils-mat/stats')

% Create the directory we will use for saving the results
if ~exist('pars', 'dir')
    mkdir('pars')
end

% Has this fit been run previously?
saveFile = sprintf('./pars/pars_%s_%d_%d.mat',type,iPtpnt,iRep);
if exist(saveFile, 'file')
    warning('Skipping this fit as it has already been completed.')
    return
end

DatSubj = DSet.P(iPtpnt).Data;
designMat = struct2DesignMat(DatSubj, 'to matrix');

opt_bads = bads;
opt_bads.NoiseFinalSamples = 100;
opt_bads.NoiseSize = 5;

opt_ibs = ibslike;
opt_ibs.Nreps = 50; %500;
opt_ibs.MaxIter = 2 * (10^4);
opt_ibs.Vectorized = 'on';

opt_varLimit = Inf; %4

RespFun = @(pars, data) aisp_simResponseWrapper(data, pars, type);
fun_handle = @(pars) ibslike_var(RespFun, pars, DatSubj.Response, designMat, ...
    opt_ibs, opt_varLimit);

[X0,LB,UB,PLB,PUB] = get_bads_bounds();

[pars,nLogL] = bads(fun_handle,X0,LB,UB,PLB,PUB,opt_bads);
save(saveFile,'pars','nLogL')


