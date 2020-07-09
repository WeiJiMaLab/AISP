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

lowerBound = [0,0,-5,0.0001];
upperBound = [20,50,5,0.5];
lowerPlausible = [2,0,-5,0.0001];
upperPlausible = [10,20,5,0.5];

x0 = lowerPlausible + rand(size(lowerBound)) .* (upperPlausible-lowerPlausible);

switch type
    case 'bayes'
        [pars,likelihood] = bads(fun_handle,x0,lowerBound,upperBound,lowerPlausible,upperPlausible,options);
        save(sprintf('~/AISP/searchShan/pars/pars_ibs_Bayes_%d_%d.mat',iSubj,iRep),'pars','likelihood')
    case 'freq'
        [pars,likelihood] = bads(fun_handle,x0,lowerBound,upperBound,lowerPlausible,upperPlausible,options);
        save(sprintf('~/AISP/searchShan/pars/pars_ibs_Freq_%d_%d.mat',iSubj,iRep),'pars','likelihood')
end