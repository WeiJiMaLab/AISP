function fit_cluster(iRep,iSubj,type)

data = getData();

addpath(genpath('../bads/'))

options = bads;
datSubj = data(data(:,1)==iSubj,:);
fun_handle = @(pars) likelihood_optim(datSubj,pars,type);
switch type
    case 'bayes'
        [pars,likelihood] = bads(fun_handle,[2,1,0,-1,-1,-1,-4,0,0.01],[-3,-3,-3,-3,-3,-3,-10,-2,eps],[6,6,6,6,6,6,10,5,0.25],[-1,-1,-1,-1,-1,-1,-3,-2,0.01],[5,4,3,3,3,2,10,5,0.1]);
        save(sprintf('~/AISP/categorizationWill/pars/parsBayes_%d_%d.mat',iSubj,iRep),'pars','likelihood')
    case 'freq'
        [pars,likelihood] = bads(fun_handle,[2,1,0,-1,-1,-1,-2,0,0.01],[-5,-5,-5,-5,-5,-5,-10,-2,eps],[5,5,5,5,5,5,10,5,0.75],[-3,-3,-3,-3,-3,-3,-4,-2,0.01],[5,4,3,3,3,2,10,5,0.1]);
        save(sprintf('~/AISP/categorizationWill/pars/parsFreq_%d_%d.mat',iSubj,iRep),'pars','likelihood')
    case 'freq2'
        [pars,likelihood] = bads(fun_handle,[2,1,0,-1,-1,-1,-2,0,0.01],[-5,-5,-5,-5,-5,-5,-10,-2,eps],[5,5,5,5,5,5,10,5,0.75],[-3,-3,-3,-3,-3,-3,-4,-2,0.01],[5,4,3,3,3,2,10,5,0.1]);
        save(sprintf('~/AISP/categorizationWill/pars/parsFreq2_%d_%d.mat',iSubj,iRep),'pars','likelihood')
    case 'freq3'
        [pars,likelihood] = bads(fun_handle,[2,1,0,-1,-1,-1,-2,0,0.01],[-5,-5,-5,-5,-5,-5,-10,-2,eps],[5,5,5,5,5,5,10,5,0.75],[-3,-3,-3,-3,-3,-3,-4,-2,0.01],[5,4,3,3,3,2,10,5,0.1]);
        save(sprintf('~/AISP/categorizationWill/pars/parsFreq3_%d_%d.mat',iSubj,iRep),'pars','likelihood')
    case 'var'
        [pars,likelihood] = bads(fun_handle,[2,1,0,-1,-1,-1,-2,0,0.01],[-5,-5,-5,-5,-5,-5,-10,-2,eps],[5,5,5,5,5,5,10,5,0.75],[-3,-3,-3,-3,-3,-3,-4,-2,0.01],[5,4,3,3,3,2,10,5,0.1]);
        save(sprintf('~/AISP/categorizationWill/pars/parsFreq3_%d_%d.mat',iSubj,iRep),'pars','likelihood')end