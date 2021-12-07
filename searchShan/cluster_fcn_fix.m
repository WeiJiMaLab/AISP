function cluster_fcn_fix(job_id, check)
fprintf('Job %d:\n', job_id)

if ~exist('check','var') || isempty(check)
    check = false;
end

runs = csvread('runs.csv');

count = 0;
for idx = randperm(size(runs,1))
    iSubj = runs(idx, 1);
    iRep = runs(idx, 2);
    iType = runs(idx, 3);
    fprintf('Run %d: Subject #%d, repetition %d, Type %d\n', idx, iSubj, iRep, iType)
    switch iType
        case 1
            file = sprintf('pars/pars_ibs_Bayes_%d_%d.mat', iSubj, iRep);
        case 2
            file = sprintf('pars/pars_ibs_Freq_%d_%d.mat', iSubj, iRep);
        case 3
            file = sprintf('pars/pars_ibs_Sample_%d_%d.mat', iSubj, iRep);
        case 4
            file = sprintf('pars/pars_ibs_cssample_%d_%d.mat', iSubj, iRep);
    end
    if ~exist(file, 'file')
        fprintf('---started fitting---\n')
        if ~check
            if iType == 1
                fit_cluster_ibs(iRep, iSubj, 'bayes')
            elseif iType == 2
                fit_cluster_ibs(iRep, iSubj, 'freq')
            elseif iType == 3
                fit_cluster_ibs(iRep, iSubj, 'sample')
            elseif iType == 4
                fit_cluster_ibs(iRep, iSubj, 'cssample')
            end
        end
        count = count + 1;
    end
end

fprintf('Done after fitting %d times', count)
