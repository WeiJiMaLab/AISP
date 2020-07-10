function fit_cluster_ibs(iRep, iSubj, type)

DSet = getData();
DatSubj = DSet.P(iSubj).Data;
designMat = struct2DesignMat(DatSubj, 'to matrix');

addpath(genpath('../bads/'))
addpath(genpath('../ibs/'))

options = bads;
options.NoiseFinalSamples = 100;
options.NoiseSize = 5;

opt_ibs = ibslike;
opt_ibs.Nreps = 50;
opt_ibs.NegLogLikeThreshold = size(DatSubj.Response, 1) * log(2) + 100;

FUN = @(pars, data) ibs_fun(data, pars, type);
fun_handle = @(pars) ibslike(FUN, pars, DatSubj.Response, designMat, opt_ibs);

[X0,LB,UB,PLB,PUB] = get_bads_bounds();

switch type
    case 'bayes'
        [pars,likelihood] = bads(fun_handle,X0,LB,UB,PLB,PUB,options);
        save(sprintf('~/AISP/searchJoshua/pars/pars_Bayes_%d_%d.mat',iSubj,iRep),'pars','likelihood')
    case 'PE'
        [pars,likelihood] = bads(fun_handle,X0,LB,UB,PLB,PUB,options);
        save(sprintf('~/AISP/searchJoshua/pars/pars_Freq_%d_%d.mat',iSubj,iRep),'pars','likelihood')
    case 'PE2'
        [pars,likelihood] = bads(fun_handle,X0,LB,UB,PLB,PUB,options);
        save(sprintf('~/AISP/searchJoshua/pars/pars_Freq2_%d_%d.mat',iSubj,iRep),'pars','likelihood')
end