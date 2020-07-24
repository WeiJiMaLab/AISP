function fit_cluster_ibs(iRep, iSubj, type)

[stimulus, response, performance] = readdata(iSubj);

data = [stimulus, response];

addpath(genpath('../bads/'))
addpath(genpath('../ibs/'))

var_limit = 4;

options = bads;
options.NoiseFinalSamples = 100;
options.NoiseSize = sqrt(var_limit);
opt_ibs = ibslike;
opt_ibs.Nreps = 1;
opt_ibs.MaxIter = 20000;
opt_ibs.Vectorized = true;
FUN = @(pars, stimulus) ibs_fun(stimulus,pars,type);
fun_handle = @(pars) ibslike_var(FUN,pars,response,stimulus,opt_ibs, var_limit);

[X0, LB, UB, PLB, PUB] = get_bads_bounds();

[pars,likelihood] = bads(fun_handle,X0,LB,UB,PLB,PUB,options);

switch type
    case 'bayes'
        save(sprintf('~/AISP/searchShan/pars/pars_ibs_Bayes_%d_%d.mat',iSubj,iRep),'pars','likelihood')
    case 'freq'
        save(sprintf('~/AISP/searchShan/pars/pars_ibs_Freq_%d_%d.mat',iSubj,iRep),'pars','likelihood')
end