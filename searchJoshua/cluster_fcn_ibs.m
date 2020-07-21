function cluster_fcn_ibs(dataDir, job_id, index)

fprintf('Job %d:\n',job_id)
fprintf('Started Job #%d\n',index)

DSet = getData(dataDir);

Nsubjects = length(DSet.P);
Nreps = 20;

if index <= Nsubjects*Nreps
    iSubj = mod(index-1,Nsubjects)+1;
    iRep = floor((index-1)/Nsubjects)+1;
    fprintf('Fitting Bayes: Subject #%d, repetition #%d\n',iSubj,iRep)
    fit_cluster_ibs(iRep,iSubj,'bayes',DSet)
elseif index <= 2*Nsubjects*Nreps
    index = index-Nsubjects*Nreps;
    iSubj = mod(index-1,Nsubjects)+1;
    iRep = floor((index-1)/Nsubjects)+1;
    fprintf('Fitting PE: Subject #%d, repetition #%d\n',iSubj,iRep)
    fit_cluster_ibs(iRep,iSubj,'PE',DSet)
elseif index <= 3*Nsubjects*Nreps
    index = index - 2*Nsubjects*Nreps;
    iSubj = mod(index-1,Nsubjects)+1;
    iRep = floor((index-1)/Nsubjects)+1;
    fprintf('Fitting PE-optimal: Subject #%d, repetition #%d\n',iSubj,iRep)
    fit_cluster_ibs(iRep,iSubj,'PE2',DSet)
end