function cluster_fcn_fancy(job_id, idx)
fprintf('Job %d:\n', job_id)
fprintf('Started Job #%d\n', idx)

Ntype = 3;
Nsubject = 9;
Nmin = 20;
slack = 2;
Ntarget = 5;

if ~exist('runs.csv','file')
    subj = 1:Nsubject;
    rep = 1:Nmin;
    type = 1:Ntype;
    [subj, rep, type] = meshgrid(subj, rep, type);
    runs = [subj(:), rep(:), type(:)];
    csvwrite('runs.csv', runs);
end

while exist('block_runs', 'file')
    pause(5)
end

runs = csvread('runs.csv');

while idx > size(runs, 1)
    % adding more lines to runs
    f = fopen('block_runs', 'w');
    for iSubj = 1:Nsubject
        for iType = 1:Ntype
            switch iType
                case 1
                    files = dir(sprintf('pars/pars_ibs_Bayes_%d_*', iSubj));
                case 2
                    files = dir(sprintf('pars/pars_ibs_Freq_%d_*', iSubj));
                case 3
                    files = dir(sprintf('pars/pars_ibs_Var_%d_*', iSubj));
            end
            evals = zeros(length(files), 1);
            for iFile = 1:length(files)
                data = load(fullfile(files(iFile).folder, files(iFile).name));
                evals(iFile) = data.likelihood;
            end
            best = min(evals);
            n_good = sum(evals < best + slack);
            if n_good < Ntarget
                iRep = max(runs(runs(:,1)==iSubj & runs(:,3)==iType,2)) + 1;
                runs = [runs; [iSubj, iRep, iType]];
            end
        end
    end
    csvwrite('runs.csv', runs);
    fclose(f);
    delete('block_runs');
end

iSubj = runs(idx, 1);
iRep = runs(idx, 2);
iType = runs(idx, 3);
fprintf('Subject #%d, repetition %d, Type %d\n', iSubj, iRep, iType)
if iType == 1
    fit_cluster_ibs(iRep, iSubj, 'bayes')
elseif iType == 2
    fit_cluster_ibs(iRep, iSubj, 'PE')
elseif iType == 3
    fit_cluster_ibs(iRep, iSubj, 'PE2')
end
