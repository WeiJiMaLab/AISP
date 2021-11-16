function [xbest, LL] = fit_cluster_ibs(imodel,isubj,irep)

% load subject, model, fitting options and bounds
load('fittingsettings.mat')

% ========= DATA/MODEL INFO ========

% fitting settings (determined by index)
subjid = subjidVec{isubj};
model = modelVec{imodel};

% load data
load(sprintf('../../changedetection/multi_item/data/fitting_data/%s_Ellipse_simple.mat',subjid))
%load(sprintf('/Volumes/GoogleDrive/My Drive/Research/VSTM/Aspen Luigi - Reliability in VWM/Exp 5 - Keshvari replication and extension/data/fitting_data/%s_Ellipse_simple.mat', subjid))
% load(sprintf('../data/fitting_data/%s_Ellipse_simple.mat',subjid))

% data in ibs format
dMat = data.Delta;
rels = unique(data.rel);
blah = data.rel;
for irel = 1:length(rels)
    blah(blah == rels(irel)) = irel;
end
dMat = [dMat blah];

% ========= FITTING INFORMATION ========

% getting model fitting settings
% logflag = logflag.(model);
% LB = LB.(model);
% UB = UB.(model);
% PLB = PLB.(model);
% PUB = PUB.(model);

% generate x0s for all reps
nvars = numel(PLB);
% x0_list = lhs(nReps,nvars,PLB,PUB,[],1e3);

% ============ FIT THE DATA =========
rng(irep);
x0 = (PUB-PLB).*rand(1,nvars)+PLB
% x0 = x0_list(irep,:);

% var_limit = 40;
fun = @(x,dMat) simulate_responses(x,model,dMat,logflag);
fun_handle = @(x) ibslike_var(fun,x,data.resp,dMat,options_ibs,var_limit);
[xbest,LL] = bads(fun_handle,x0,LB,UB,PLB,PUB,[],options)


xbest(logflag) = exp(xbest(logflag)); % getting parameters back into natural units

if ~nargout
save(sprintf('fits/model%s_subj%s_rep%d.mat',model,subjid,irep),...
    'options','xbest','LL','var_limit')
end