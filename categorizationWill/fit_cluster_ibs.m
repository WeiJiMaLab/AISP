function fit_cluster_ibs(iRep,iSubj,type)

data = getData();

addpath(genpath('../bads/'))
addpath(genpath('../ibs/'))

var_limit = 4;

options = bads;
options.NoiseFinalSamples = 100;
options.NoiseSize = sqrt(var_limit);
datSubj = data(data(:,1)==iSubj,:);
opt_ibs = ibslike;
opt_ibs.Nreps = 1;
opt_ibs.MaxIter = 20000;
FUN = @(pars,data) ibs_fun(data,pars,type);
fun_handle = @(pars) ibslike_var_par(FUN,pars,datSubj(:,4),datSubj,opt_ibs, var_limit);

[X0, LB, UB, PLB, PUB] = get_bads_bounds();

switch type
    case 'bayes'
        fname = sprintf('~/AISP/categorizationWill/pars/pars_ibs_Bayes_%d_%d.mat',iSubj,iRep);
        if ~exist(fname, 'file')
            [pars,likelihood] = bads(fun_handle,X0,LB,UB,PLB,PUB,options);
            save(fname,'pars','likelihood')
        end
    case 'freq'
        fname = sprintf('~/AISP/categorizationWill/pars/pars_ibs_Freq_%d_%d.mat',iSubj,iRep);
        if ~exist(fname, 'file')
            [pars,likelihood] = bads(fun_handle,X0,LB,UB,PLB,PUB,options);
            save(fname,'pars','likelihood')
        end
    case 'freq2'
        fname = sprintf('~/AISP/categorizationWill/pars/pars_ibs_Freq2_%d_%d.mat',iSubj,iRep);
        if ~exist(fname, 'file')
            [pars,likelihood] = bads(fun_handle,X0,LB,UB,PLB,PUB,options);
            save(fname,'pars','likelihood')
        end
    case 'freq3'
        fname = sprintf('~/AISP/categorizationWill/pars/pars_ibs_Freq3_%d_%d.mat',iSubj,iRep);
        if ~exist(fname, 'file')
            [pars,likelihood] = bads(fun_handle,X0,LB,UB,PLB,PUB,options);
            save(fname ,'pars','likelihood')
        end
    case 'var'
        fname = sprintf('~/AISP/categorizationWill/pars/pars_ibs_Var_%d_%d.mat',iSubj,iRep);
        if ~exist(fname, 'file')
            [pars,likelihood] = bads(fun_handle,X0,LB,UB,PLB,PUB,options);
            save(fname,'pars','likelihood')
        end
    case 'sample'
        fname = sprintf('~/AISP/categorizationWill/pars/pars_ibs_sample_%d_%d.mat',iSubj,iRep);
        if ~exist(fname, 'file')
            [pars,likelihood] = bads(fun_handle,[X0,10],[LB,1],[UB,1000],[PLB,1],[PUB, 100],options);
            save(fname, 'pars','likelihood')
        end
    case 'cssample'
        fname = sprintf('~/AISP/categorizationWill/pars/pars_ibs_cssample_%d_%d.mat',iSubj,iRep);
        if ~exist(fname, 'file')
            [pars,likelihood] = bads(fun_handle,[X0,10],[LB,1],[UB,1000],[PLB,1],[PUB, 100],options);
            save(fname,'pars','likelihood')
        end
end