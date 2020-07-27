function fit_cluster_ibs(iRep, iSubj, type, DSet)

addpath(genpath('../bads/'))
addpath(genpath('../ibs/'))
addpath('./visualSearch')
addpath('../lautils-mat/stats')

DatSubj = DSet.P(iSubj).Data;
designMat = struct2DesignMat(DatSubj, 'to matrix');

options = bads;
options.NoiseFinalSamples = 100;
options.NoiseSize = 5;

opt_ibs = ibslike;
opt_ibs.Nreps = 25;
opt_ibs.MaxIter = 10000;
opt_ibs.Vectorized = 'on';

FUN = @(pars, data) ibs_fun(data, pars, type);
fun_handle = @(pars) ibslike(FUN, pars, DatSubj.Response, designMat, opt_ibs);

[X0,LB,UB,PLB,PUB] = get_bads_bounds();

switch type
    case 'bayes'
        [pars,nLogL] = bads(fun_handle,X0,LB,UB,PLB,PUB,options);
        save(sprintf('./pars/pars_Bayes_%d_%d.mat',iSubj,iRep),'pars','nLogL')
    case 'PE'
        [pars,nLogL] = bads(fun_handle,X0,LB,UB,PLB,PUB,options);
        save(sprintf('./pars/pars_Freq_%d_%d.mat',iSubj,iRep),'pars','nLogL')
    case 'PE2'
        [pars,nLogL] = bads(fun_handle,X0,LB,UB,PLB,PUB,options);
        save(sprintf('./pars/pars_Freq2_%d_%d.mat',iSubj,iRep),'pars','nLogL')
end