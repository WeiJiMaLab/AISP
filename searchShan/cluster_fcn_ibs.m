function cluster_fcn_ibs(job_id,index)
fprintf('Job %d:\n',job_id)
fprintf('Started Job #%d\n',index)

Nsubjects = 9;
Nreps = 20;
if index <= Nsubjects*Nreps
    iSubj = mod(index-1,Nsubjects)+1;
    iRep = floor((index-1)/Nsubjects)+1;
    fprintf('Fitting Bayes: Subject #%d, repetition #%d\n',iSubj,iRep)
    fit_cluster_ibs(iRep,iSubj,'bayes')
elseif index <= 2*Nsubjects*Nreps
    index = index-Nsubjects*Nreps;
    iSubj = mod(index-1,Nsubjects)+1;
    iRep = floor((index-1)/Nsubjects)+1;
    fprintf('Fitting PE: Subject #%d, repetition #%d\n',iSubj,iRep)
    fit_cluster_ibs(iRep,iSubj,'freq')
end