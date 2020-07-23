function fit_cluster_ibs(iRep, iSubj, type)

data = getData();

addpath(genpath('../bads/'))
addpath(genpath('../ibs/'))

var_limit = 4;

options = bads;
options.NoiseFinalSamples = 100;
options.NoiseSize = sqrt(var_limit);
datSubj = data(data(:,1)==iSubj,:);
opt_ibs = ibslike;
opt_ibs.Nreps = 1;
opt_ibs.MaxIter = 100000;
FUN = @(pars,data) ibs_fun(data,pars,type);
fun_handle = @(pars) ibslike_var(FUN,pars,datSubj(:,4),datSubj,opt_ibs, var_limit);

[X0,LB,UB,PLB,PUB] = get_bads_bounds();

switch type
    case 'bayes'
        [pars,likelihood] = bads(fun_handle,X0,LB,UB,PLB,PUB,options);
        save(sprintf('~/AISP/collinearityYanli/pars/pars_Bayes_%d_%d.mat',iSubj,iRep),'pars','likelihood')
    case 'PE'
        [pars,likelihood] = bads(fun_handle,X0,LB,UB,PLB,PUB,options);
        save(sprintf('~/AISP/collinearityYanli/pars/pars_Freq_%d_%d.mat',iSubj,iRep),'pars','likelihood')
    case 'PE2'
        [pars,likelihood] = bads(fun_handle,X0,LB,UB,PLB,PUB,options);
        save(sprintf('~/AISP/collinearityYanli/pars/pars_Freq2_%d_%d.mat',iSubj,iRep),'pars','likelihood')
end