function resp = cssample_simulate3(x,dMat,logflag, tempp)
%function RESP = cssample_simulate(X,MODEL,DMAT,LOGFLAG) simulates responses of
%joint posterior sampling observer
%
% ============ INPUT VARIABLES ============
% X: parameter values. vector of length seven
% DMAT: data. nTrials x 8 matrix
%       first four columns correspond to amount of orientation change for
%       each item. second four columns correspond to the reliability of
%       each item.
% LOGFLAG: log flag. binary vector of length six
%       indicates which parameters are in log scaling
% tempp: struct of a look-up table for sampling precisions
%       needs to be passed to avoid persistent variables
%
% ============ OUTPUT VARIABLES ============
% RESP: length nTrials vector of simulated responses

%persistent k_range
%persistent J_lin
%persistent highest_J

if nargin < 3; logflag = []; end

x(logflag) = exp(x(logflag));

% start off with lapse
% lapserate = x(end);
% nTrials = size(dMat,1);
% islapse = rand(1,nTrials) < lapserate;
% lapserespVec = rand(1,sum(islapse)) > 0.5;      % for lapse trials, flip a coin
% resp = nan(length(islapse),1);
% resp(islapse) = lapserespVec;
%
% if sum(~islapse) % if there are any trials that did not lapse
%
%     % reduce dMat to only include nonlapse trials
%     dMat = dMat(~islapse,:);

% define data stuff
nTrials = size(dMat,1);
nItems = 4;
Delta = dMat(:,1:nItems);           % amount change for each of four items
Rels = dMat(:,(nItems+1):end);      % reliabilities for each item (1: low, 2: high)
nRelsVec = sum(Rels==2,2);

% ===== GET PARAMETER VALUES ======
Jbar_high = x(1);       % mean precision of high rel ellipse
Jbar_low = x(2);        % mean precision of low rel ellipse
tau = x(3);             % scale parameter of ellipse
beta = x(4);            % softmax temperature parameter (deciison noise)
beta0 = x(5);           % bias 
lambda = x(6);          % lapse rate
nSamples = x(7);       % n_samples

% random averaging between n_samples values
if rand > mod(nSamples, 1)
    nSamples = floor(nSamples);
else
    nSamples = ceil(nSamples);
end
    
% ====== CALCULATE P(\HAT{C}==1|\Theta) FOR nSamples SAMPLES =====

% make CDF for interpolating J to Kappa
%if isempty(k_range)
    % tempp = load('cdf_table.mat');
    % K_interp = tempp.K_interp;
    % cdf = tempp.cdf;
k_range = tempp.k_range;
J_lin = tempp.J_lin;
highest_J = tempp.highest_J;
%end

% calculate actual kappa and noisy representations
Jbar_mat = Rels;
Jbar_mat(Rels==1) = Jbar_low;
Jbar_mat(Rels==2) = Jbar_high;

J_x_mat = gamrnd(Jbar_mat./tau,tau);
J_y_mat = gamrnd(Jbar_mat./tau,tau);

% set kappas too high to highest J (alternatively can resample, as
% keshvari did)
J_x_mat(J_x_mat > highest_J) = highest_J;
J_y_mat(J_y_mat > highest_J) = highest_J;

% convert J to kappa
tempvar = 1/diff(J_lin(1:2))*J_x_mat+1;
kappa_x = k_range(round(tempvar));
tempvar = 1/diff(J_lin(1:2))*J_y_mat+1;
kappa_y = k_range(round(tempvar));

if size(kappa_x,2) ~= nItems
    kappa_x = kappa_x';
    kappa_y = kappa_y';
end

% stimuli
% xi = 0;
xi = rand(nTrials,nItems)*2*pi-pi;

% measurement noise
x = circ_vmrnd(0,kappa_x);
y = circ_vmrnd(0,kappa_y);
x = x + xi;
y = y + xi + Delta;

% first sample of s and C
prop_xi = circ_vmrnd(x, kappa_x);
DeltaMat = zeros(nItems,nTrials);
DeltaVec = (rand(nTrials,1).*2*pi)-pi; % all Deltas
C_samp = rand(nTrials,1)<0.5;
DeltaVec(C_samp)=0;
idx = randi(4,[1 nTrials])+(0:4:(nItems*nTrials-nItems));
DeltaMat(idx) = DeltaVec;
DeltaMat = DeltaMat';
C = nan(nTrials,nSamples);
log_p = sum(kappa_y .* cos(y - prop_xi - DeltaMat), 2);

p_accept = zeros(nSamples, 1);
for isamp = 1:nSamples
    
    % proposal sample
    prop_xi = circ_vmrnd(x, kappa_x);
    
    prop_DeltaMat = zeros(nItems,nTrials);
    DeltaVec = (rand(nTrials,1).*2.*pi)-pi; % all Deltas
    C_samp_new = rand(nTrials,1)<0.5;
    DeltaVec(C_samp_new)=0; % set with 0.5 some trials to 0 delta
    idx = randi(4,[1 nTrials])+(0:4:(nItems*nTrials-nItems));
    prop_DeltaMat(idx) = DeltaVec;
    prop_DeltaMat = prop_DeltaMat';
    
    % calculate acceptance p
    log_p_new = sum(kappa_y .* cos(y - prop_xi + prop_DeltaMat), 2);
    accept = rand(nTrials,1) < exp(log_p_new - log_p);
    
    % update relevant trial samples
    C_samp(accept) = C_samp_new(accept);
    log_p(accept) = log_p_new(accept);
    % DeltaMat(accept,:) = prop_DeltaMat(accept,:);
    % xi(accept,:) = prop_xi(accept,:);

    C(:,isamp) = C_samp;
    p_accept(isamp) = mean(accept);
end

d = log(sum(C, 2) + 1) - log(nSamples - sum(C, 2) + 1);

p = lambda/2 + (1-lambda)./(1+exp(beta0+beta.*d));
resp = (rand(nTrials, 1) < p);

