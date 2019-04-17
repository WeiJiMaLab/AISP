data = getData();

addpath(genpath('../bads/'))

Nreps = 10;

options = bads;
pars = zeros(numel(unique(data(:,1))),9,Nreps);
likelihoods = zeros(numel(unique(data(:,1))),Nreps);
k = 0;
for iSubj = unique(data(:,1))'
    k = k + 1;
    datSubj = data(data(:,1)==iSubj,:);
    for iRep = 1:Nreps
        fun_handle = @(pars) likelihood_optim(datSubj,pars,'bayes');
        [parSub,l] = bads(fun_handle,[2,1,0,-1,-1,-1,-4,0,0.01],[-3,-3,-3,-3,-3,-3,-10,-2,eps],[6,6,6,6,6,6,10,5,0.25],[-1,-1,-1,-1,-1,-1,-3,-2,0.01],[5,4,3,3,3,2,10,5,0.1]);
        pars(iSubj,:,iRep) = parSub;
        likelihoods(iSubj,iRep) = l;
    end
end

save('pars/parsBayes.mat','pars','likelihoods')


options = bads;
pars = zeros(numel(unique(data(:,1))),9,Nreps);
likelihoods = zeros(numel(unique(data(:,1))),Nreps);
k = 0;
for iSubj = unique(data(:,1))'
    k = k + 1;
    datSubj = data(data(:,1)==iSubj,:);
    for iRep = 1:Nreps
        fun_handle = @(pars) likelihood_optim(datSubj,pars,'freq');
        [parSub,l] = bads(fun_handle,[2,1,0,-1,-1,-1,-2,0,0.01],[-5,-5,-5,-5,-5,-5,-10,-2,eps],[5,5,5,5,5,5,10,5,0.75],[-3,-3,-3,-3,-3,-3,-4,-2,0.01],[5,4,3,3,3,2,10,5,0.1]);
        pars(iSubj,:,iRep) = parSub;
        likelihoods(iSubj,iRep) = l;
    end
end

save('pars/parsFreq.mat','pars','likelihoods')



options = bads;
pars = zeros(numel(unique(data(:,1))),9,Nreps);
likelihoods = zeros(numel(unique(data(:,1))),Nreps);
k = 0;
for iSubj = unique(data(:,1))'
    k = k + 1;
    datSubj = data(data(:,1)==iSubj,:);
    for iRep = 1:Nreps
        fun_handle = @(pars) likelihood_optim(datSubj,pars,'freq2');
        [parSub,l] = bads(fun_handle,[2,1,0,-1,-1,-1,-2,0,0.01],[-5,-5,-5,-5,-5,-5,-10,-2,eps],[5,5,5,5,5,5,10,5,0.75],[-3,-3,-3,-3,-3,-3,-4,-2,0.01],[5,4,3,3,3,2,10,5,0.1]);
        pars(iSubj,:,iRep) = parSub;
        likelihoods(iSubj,iRep) = l;
    end
end

save('pars/parsFreq2.mat','pars','likelihoods')


options = bads;
pars = zeros(numel(unique(data(:,1))),9,Nreps);
likelihoods = zeros(numel(unique(data(:,1))),Nreps);
k = 0;
for iSubj = unique(data(:,1))'
    k = k + 1;
    datSubj = data(data(:,1)==iSubj,:);
    for iRep = 1:Nreps
        fun_handle = @(pars) likelihood_optim(datSubj,pars,'freq3');
        [parSub,l] = bads(fun_handle,[2,1,0,-1,-1,-1,-2,0,0.01],[-5,-5,-5,-5,-5,-5,-10,-2,eps],[5,5,5,5,5,5,10,5,0.75],[-3,-3,-3,-3,-3,-3,-4,-2,0.01],[5,4,3,3,3,2,10,5,0.1]);
        pars(iSubj,:,iRep) = parSub;
        likelihoods(iSubj,iRep) = l;
        fprintf('finished run %d of subject %d\n',iRep,iSubj)
    end
end

save('pars/parsFreq3.mat','pars','likelihoods')