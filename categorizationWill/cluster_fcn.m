function cluster_fcn(job_id,index)
fprintf('Job %d:\n',job_id)
fprintf('Started Job #%d\n',index)

Nsubjects = 11;
Nreps = 10;
if index <= Nsubjects*Nreps
    iSubj = mod(index-1,Nsubjects)+1;
    iRep = floor((index-1)/Nsubjects)+1;
    fprintf('Fitting Bayes: Subject #%d, repetition #%d\n',iSubj,iRep)
    fit_cluster(iRep,iSubj,'bayes')
elseif index <= 2*Nsubjects*Nreps
    index = index-Nsubjects*Nreps;
    iSubj = mod(index-1,Nsubjects)+1;
    iRep = floor((index-1)/Nsubjects)+1;
    fprintf('Fitting PE: Subject #%d, repetition #%d\n',iSubj,iRep)
    fit_cluster(iRep,iSubj,'freq')
elseif index <= 3*Nsubjects*Nreps
    index = index-2*Nsubjects*Nreps;
    iSubj = mod(index-1,Nsubjects)+1;
    iRep = floor((index-1)/Nsubjects)+1;
    fprintf('Fitting PE2: Subject #%d, repetition #%d\n',iSubj,iRep)
    fit_cluster(iRep,iSubj,'freq2')
elseif index <= 4*Nsubjects*Nreps
    index = index-3*Nsubjects*Nreps;
    iSubj = mod(index-1,Nsubjects)+1;
    iRep = floor((index-1)/Nsubjects)+1;
    fprintf('Fitting PE3: Subject #%d, repetition #%d\n',iSubj,iRep)
    fit_cluster(iRep,iSubj,'freq3')
elseif index <= 5*Nsubjects*Nreps
    index = index-4*Nsubjects*Nreps;
    iSubj = mod(index-1,Nsubjects)+1;
    iRep = floor((index-1)/Nsubjects)+1;
    fprintf('Fitting Variational: Subject #%d, repetition #%d\n',iSubj,iRep)
    fit_cluster(iRep,iSubj,'var')
end