function fit_cluster_ibs(iRep, iSubj, type)

[stimulus, response, performance] = readdata(iSubj);

data = [stimulus, response];

addpath(genpath('../bads/'))
addpath(genpath('../ibs/'))

options = bads;
options.NoiseFinalSamples = 100;
options.NoiseSize = 5;
opt_ibs = ibslike;
opt_ibs.Nreps = 25;
opt_ibs.NegLogLikeThreshold = size(data,1)*log(2)+200;
%fun_handle = @(pars) likelihood_optim(datSubj,pars,type);
FUN = @(pars, data) ibs_fun(data,pars,type);
fun_handle = @(pars) ibslike(FUN,pars,response,stimulus,opt_ibs);
switch type
    case 'bayes'
        [pars,likelihood] = bads(fun_handle,[5,1,0,0.01],[0,0,-5,0],[20,50,5,0.5],[2,0,-5,0],[10,20,5,0.5],options);
        save(sprintf('~/AISP/searchShan/pars/pars_ibs_Bayes_%d_%d.mat',iSubj,iRep),'pars','likelihood')
    case 'freq'
        [pars,likelihood] = bads(fun_handle,[1,1,0,0.01],[0,0,-5,0],[20,50,5,0.5],[2,0,-5,0],[10,20,5,0.5],options);
        save(sprintf('~/AISP/searchShan/pars/pars_ibs_Freq_%d_%d.mat',iSubj,iRep),'pars','likelihood')
end