function fit_cluster_ibs(iRep,iSubj,type)

data = getData();

addpath(genpath('../bads/'))
addpath(genpath('../ibs/'))

options = bads;
options.NoiseFinalSamples = 100;
options.NoiseSize = 5;
datSubj = data(data(:,1)==iSubj,:);
opt_ibs = ibslike;
opt_ibs.Nreps = 25;
opt_ibs.NegLogLikeThreshold = size(datSubj,1) * log(2) + 100;
%fun_handle = @(pars) likelihood_optim(datSubj,pars,type);
FUN = @(pars,data) ibs_fun(data,pars,type);
fun_handle = @(pars) ibslike(FUN,pars,datSubj(:,4),datSubj,opt_ibs);
switch type
    case 'bayes'
        [pars,likelihood] = bads(fun_handle,[2,1,0,-1,-1,-1,-4,0,0.01],[-3,-3,-3,-3,-3,-3,-10,-2,eps],[6,6,6,6,6,6,10,5,0.25],[-1,-1,-1,-1,-1,-1,-3,-2,0.01],[5,4,3,3,3,2,10,5,0.1],options);
        save(sprintf('~/AISP/categorizationWill/pars/pars_ibs_Bayes_%d_%d.mat',iSubj,iRep),'pars','likelihood')
    case 'freq'
        [pars,likelihood] = bads(fun_handle,[2,1,0,-1,-1,-1,-2,0,0.01],[-5,-5,-5,-5,-5,-5,-10,-2,eps],[5,5,5,5,5,5,10,5,0.75],[-3,-3,-3,-3,-3,-3,-4,-2,0.01],[5,4,3,3,3,2,10,5,0.1],options);
        save(sprintf('~/AISP/categorizationWill/pars/pars_ibs_Freq_%d_%d.mat',iSubj,iRep),'pars','likelihood')
    case 'freq2'
        [pars,likelihood] = bads(fun_handle,[2,1,0,-1,-1,-1,-2,0,0.01],[-5,-5,-5,-5,-5,-5,-10,-2,eps],[5,5,5,5,5,5,10,5,0.75],[-3,-3,-3,-3,-3,-3,-4,-2,0.01],[5,4,3,3,3,2,10,5,0.1],options);
        save(sprintf('~/AISP/categorizationWill/pars/pars_ibs_Freq2_%d_%d.mat',iSubj,iRep),'pars','likelihood')
    case 'freq3'
        [pars,likelihood] = bads(fun_handle,[2,1,0,-1,-1,-1,-2,0,0.01],[-5,-5,-5,-5,-5,-5,-10,-2,eps],[5,5,5,5,5,5,10,5,0.75],[-3,-3,-3,-3,-3,-3,-4,-2,0.01],[5,4,3,3,3,2,10,5,0.1],options);
        save(sprintf('~/AISP/categorizationWill/pars/pars_ibs_Freq3_%d_%d.mat',iSubj,iRep),'pars','likelihood')
    case 'var'
        [pars,likelihood] = bads(fun_handle,[2,1,0,-1,-1,-1,-2,0,0.01],[-5,-5,-5,-5,-5,-5,-10,-2,eps],[5,5,5,5,5,5,10,5,0.75],[-3,-3,-3,-3,-3,-3,-4,-2,0.01],[5,4,3,3,3,2,10,5,0.1],options);
        save(sprintf('~/AISP/categorizationWill/pars/pars_ibs_Var_%d_%d.mat',iSubj,iRep),'pars','likelihood')
end