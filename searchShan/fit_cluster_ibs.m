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


switch type
    case 'bayes'
        [pars,likelihood] = bads(fun_handle,X0,LB,UB,PLB,PUB,options);
        save(sprintf('~/AISP/searchShan/pars/pars_ibs_Bayes_%d_%d.mat',iSubj,iRep),'pars','likelihood')
    case 'freq'
        [pars,likelihood] = bads(fun_handle,X0,LB,UB,PLB,PUB,options);
        save(sprintf('~/AISP/searchShan/pars/pars_ibs_Freq_%d_%d.mat',iSubj,iRep),'pars','likelihood')
    case 'sample1'
        [pars,likelihood] = bads(fun_handle,X0,LB,UB,PLB,PUB,options);
        save(sprintf('~/AISP/searchShan/pars/pars_ibs_Sample1_%d_%d.mat',iSubj,iRep),'pars','likelihood')
    case 'sample'
        [pars,likelihood] = bads(fun_handle,[X0,10],[LB,1],[UB,10000],[PLB,1],[PUB,1000],options);
        save(sprintf('~/AISP/searchShan/pars/pars_ibs_Sample_%d_%d.mat',iSubj,iRep),'pars','likelihood')
    case 'cssample'
        [pars,likelihood] = bads(fun_handle,[X0,10],[LB,1],[UB,10000],[PLB,1],[PUB,1000],options);
        save(sprintf('~/AISP/searchShan/pars/pars_ibs_cssample_%d_%d.mat',iSubj,iRep),'pars','likelihood')
end