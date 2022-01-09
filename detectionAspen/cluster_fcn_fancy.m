function cluster_fcn_fancy(idx)
fprintf('Started job idx %d\n', idx)

% fitting settings
load('fittingsettings.mat')

% get all minimum runs, if it doesn't already exist
if ~exist('runs.csv','file')
    subj = 1:nSubjs;
    rep = 1:nMinReps;
    type = 1:nModels;
    [subj, rep, type] = meshgrid(subj, rep, type);
    runs = [subj(:), rep(:), type(:)];
    csvwrite('runs.csv', runs);
end

% infinite loop if this file exists
while exist('block_runs', 'file')
    fprintf('block_runs exists')
    pause(5)
end

% load runs file
runs = csvread('runs.csv');

% create more jobs if your idx is larger than current size of runs
% will only create more jobs if prev didn't hit nTargets within slackLL of minLL
while idx > size(runs, 1)
    f = fopen('block_runs', 'w');
    for isubj = 1:nSubjs
        subjid = subjidVec{isubj}
        
        for imodel = 1:nModels
            model = modelVec{imodel};
            
            % get all files for this model and subjid
            files = dir(sprintf('fits/model%s_subj%s_rep*',model,subjid));

            LLVec = zeros(length(files), 1);
            for ifile = 1:length(files)
                data = load(sprintf('fits/%s',files(ifile).name));
                LLVec(ifile) = data.LL;
            end
            bestLL = min(LLVec);
            nGood = sum(LLVec < bestLL+slackLL);
            if nGood < nTargets
                irep = max(runs(runs(:,1)==isubj & runs(:,3)==imodel,2)) + 1;
                runs = [runs; [isubj, irep, imodel]];
            else
                sprintf('satisfied!')
            end
        end
    end
    csvwrite('runs.csv', runs);
    fclose(f);
    delete('block_runs');
end

% fit current idx
imodel = runs(idx, 3);
isubj = runs(idx, 1);
irep = runs(idx, 2);

fprintf('Subject: %d (%s), Model: %s, Repetition: %d\n', isubj, subjidVec{isubj}, modelVec{imodel}, irep)
fit_cluster_ibs(imodel,isubj,irep)