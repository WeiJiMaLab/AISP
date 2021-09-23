function cluster_fcn_ibs(dataDir, job_id, index)

fprintf('Job %d:\n',job_id)
fprintf('Started Job #%d\n',index)

[DSet, Nptpnts] = getData(dataDir);
Config = produceConfig();
Nreps = Config.Nreps;

if index <= Nptpnts*Nreps
    iSubj = mod(index-1,Nptpnts)+1;
    iRep = floor((index-1)/Nptpnts)+1;
    fprintf('Fitting Bayes: Subject #%d, repetition #%d\n',iSubj,iRep)
    fit_cluster_ibs(iRep,iSubj,'Bayes',DSet)
elseif index <= 2*Nptpnts*Nreps
    index = index-Nptpnts*Nreps;
    iSubj = mod(index-1,Nptpnts)+1;
    iRep = floor((index-1)/Nptpnts)+1;
    fprintf('Fitting PE: Subject #%d, repetition #%d\n',iSubj,iRep)
    fit_cluster_ibs(iRep,iSubj,'PE',DSet)
elseif index <= 3*Nptpnts*Nreps
    index = index - 2*Nptpnts*Nreps;
    iSubj = mod(index-1,Nptpnts)+1;
    iRep = floor((index-1)/Nptpnts)+1;
    fprintf('Fitting PE-optimal: Subject #%d, repetition #%d\n',iSubj,iRep)
    fit_cluster_ibs(iRep,iSubj,'PE2',DSet)
end