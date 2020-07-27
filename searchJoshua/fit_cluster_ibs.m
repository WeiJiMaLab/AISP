function fit_cluster_ibs(iRep, iSubj, type, DSet)

addpath(genpath('../bads/'))
addpath(genpath('../ibs/'))
addpath('./visualSearch')
addpath('../lautils-mat/stats')

DatSubj = DSet.P(iSubj).Data;
designMat = struct2DesignMat(DatSubj, 'to matrix');

opt_bads = bads;
opt_bads.NoiseFinalSamples = 100;
opt_bads.NoiseSize = 5;

opt_ibs = ibslike;
opt_ibs.Nreps = 25;
opt_ibs.MaxIter = 2 * (10^4);
opt_ibs.Vectorized = 'on';

opt_varLimit = 4;

FUN = @(pars, data) ibs_fun(data, pars, type);
fun_handle = @(pars) ibslike_var(FUN, pars, DatSubj.Response, designMat, ...
    opt_ibs, opt_varLimit);

[X0,LB,UB,PLB,PUB] = get_bads_bounds();

switch type
    case 'bayes'
        [pars,nLogL] = bads(fun_handle,X0,LB,UB,PLB,PUB,opt_bads);
        save(sprintf('./pars/pars_Bayes_%d_%d.mat',iSubj,iRep),'pars','nLogL')
    case 'PE'
        [pars,nLogL] = bads(fun_handle,X0,LB,UB,PLB,PUB,opt_bads);
        save(sprintf('./pars/pars_Freq_%d_%d.mat',iSubj,iRep),'pars','nLogL')
    case 'PE2'
        [pars,nLogL] = bads(fun_handle,X0,LB,UB,PLB,PUB,opt_bads);
        save(sprintf('./pars/pars_Freq2_%d_%d.mat',iSubj,iRep),'pars','nLogL')
end