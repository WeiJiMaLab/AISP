function cluster_fcn_fix(job_id, check)
fprintf('Job %d:\n', job_id)

if ~exist('check','var') || isempty(check)
    check = false;
end

runs = csvread('runs.csv');

for idx = randperm(size(runs,1))
    iSubj = runs(idx, 1);
    iRep = runs(idx, 2);
    iType = runs(idx, 3);
    fprintf('Run %d: Subject #%d, repetition %d, Type %d\n', idx, iSubj, iRep, iType)
    switch iType
        case 1
            file = sprintf('pars/pars_Bayes_%d_%d.mat', iSubj, iRep);
        case 2
            file = sprintf('pars/pars_Freq_%d_%d.mat', iSubj, iRep);
        case 3
            file = sprintf('pars/pars_Freq2_%d_%d.mat', iSubj, iRep);
        case 4
            file = sprintf('pars/pars_Sample_%d_%d.mat', iSubj, iRep);
        case 5
            file = sprintf('pars/pars_cssample_%d_%d.mat', iSubj, iRep);
    end
    if ~exist(file, 'file')
        fprintf('---started fitting---\n')
        if ~check
            if iType == 1
                fit_cluster_ibs(iRep, iSubj, 'bayes')
            elseif iType == 2
                fit_cluster_ibs(iRep, iSubj, 'PE')
            elseif iType == 3
                fit_cluster_ibs(iRep, iSubj, 'PE2')
            elseif iType == 4
                fit_cluster_ibs(iRep, iSubj, 'sample')
            elseif iType == 5
                fit_cluster_ibs(iRep, iSubj, 'cssample')
            end
        end
    end
end

