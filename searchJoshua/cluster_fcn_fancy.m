function cluster_fcn_fancy(dataDir, job_id, idx)

% Settings
slack = 2;
N_target = 5;


fprintf('Job %d:\n', job_id)
fprintf('Started Job #%d\n', idx)

[DSet, Nptpnts] = getData(dataDir);
Config = load('Config.mat');
NMinReps = Config.Nreps;
Ntype = length(Config.ModelList);

mkdir('tmp')
if ~exist('./tmp/runs.csv','file')
    ptpnt = 1:Nptpnts;
    rep = 1:NMinReps;
    type = 1:Ntype;
    [ptpnt, rep, type] = meshgrid(ptpnt, rep, type);
    runs = [ptpnt(:), rep(:), type(:)];
    csvwrite('./tmp/runs.csv', runs);
end

while exist('./tmp/block_runs', 'file')
    pause(5)
end

runs = csvread('./tmp/runs.csv');

while idx > size(runs, 1)
    % adding more lines to runs
    f = fopen('./tmp/block_runs', 'w');
    for iPtpnt = 1:Nptpnts
        for iType = 1:Ntype
            files = dir(sprintf('pars/pars_%s_%d_*', ...
                Config.ModelList{iType}, iPtpnt));
            
            evals = nan(length(files), 1);
            for iFile = 1:length(files)
                data = load(fullfile(files(iFile).folder, files(iFile).name));
                evals(iFile) = data.nLogL;
            end
   %         disp(['***** Evaluations found for participant ' num2str(iPtpnt) ' model ' num2str(iType) ': '])
   %         disp(evals)
            best = min(evals);
            n_good = sum(evals < (best + slack));
            disp(['***** Num good runs for participant ' num2str(iPtpnt) ' model ' num2str(iType) ': ' num2str(n_good)])
            if n_good < N_target
                iRep = max(runs(runs(:,1)==iPtpnt & runs(:,3)==iType,2)) + 1;
                runs = [runs; [iPtpnt, iRep, iType]];
            end
        end
    end
    csvwrite('./tmp/runs.csv', runs);
    fclose(f);
    delete('./tmp/block_runs');
end

iPtpnt = runs(idx, 1);
iRep = runs(idx, 2);
iType = runs(idx, 3);
fprintf('Participant #%d, Repetition %d, Type %d\n', iPtpnt, iRep, iType)
fit_cluster_ibs(iRep, iPtpnt, Config.ModelList{iType}, DSet, idx)

