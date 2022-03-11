function resp = sample_simulate(theta,dMat,logflag)
% function resp = sample_simulate(x,dMat,logflag,tempp)
%function RESP = sample_simulate(X,MODEL,DMAT,LOGFLAG) simulates responses of
%sampling observer
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
Jbar_high = theta(1);       % mean precision of high rel ellipse
Jbar_low = theta(2);        % mean precision of low rel ellipse
tau = theta(3);             % scale parameter of ellipse
beta = theta(4);            % softmax temperature parameter (deciison noise)
beta0 = theta(5);           % bias 
lambda = theta(6);          % lapse rate
nSamples = theta(7);       % n_samples

% random averaging between n_samples values
if rand > mod(nSamples, 1)
    nSamples = floor(nSamples);
else
    nSamples = ceil(nSamples);
end
    
% ====== CALCULATE P(\HAT{C}==1|\Theta) FOR nSamples SAMPLES =====

% make CDF for interpolating J to Kappa
%if isempty(k_range)
%    tempp = load('cdf_table.mat');
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

% replicate noise matrices based on number of samples
Jbar_mat = repmat(Jbar_mat,[1 1 nSamples]); %ntrials, nitems, nsamples

% sample first display and second display stimuli for each C
s1_DeltaMat = zeros(nItems,nTrials,nSamples);
DeltaVec = (rand(nTrials,1,nSamples).*2*pi)-pi; % all Deltas
idx = randi(4,[1 nTrials*nSamples])+(0:4:(nItems*nTrials*nSamples-nItems));
s1_DeltaMat(idx) = DeltaVec;
s1_DeltaMat = permute(s1_DeltaMat,[2 1 3]);

J_x_mat = gamrnd(Jbar_mat./tau,tau);
J_y_mat = gamrnd(Jbar_mat./tau,tau);

% set kappas too high to highest J (alternatively can resample, as
% keshvari did)
J_x_mat(J_x_mat > highest_J) = highest_J;
J_y_mat(J_y_mat > highest_J) = highest_J;

% convert J to kappa
xi = 1/diff(J_lin(1:2))*J_x_mat+1;
kappa_x = k_range(round(xi));
xi = 1/diff(J_lin(1:2))*J_y_mat+1;
kappa_y = k_range(round(xi));

if size(kappa_x,2) ~= nItems
    kappa_x = kappa_x';
    kappa_y = kappa_y';
end

% measurement noise
x = circ_vmrnd(0,kappa_x);
y_0 = circ_vmrnd(0,kappa_y);
y_1 = bsxfun(@plus,y_0,Delta);

% terms for multiplication of von mises, for both sampled category conditions
mu_0 = x + atan2(sin(y_0-x),kappa_x./kappa_y + cos(y_0-x));
mu_1 = x + atan2(sin(y_1-x),kappa_x./kappa_y + cos(x-(y_1-s1_DeltaMat)));

kappa_0 = sqrt(kappa_x.^2 + kappa_y.^2 + (2*kappa_x.*kappa_y.*cos(x-y_0)));
kappa_1 = sqrt(kappa_x.^2 + kappa_y.^2 + (2*kappa_x.*kappa_y.*cos(x-(y_1-s1_DeltaMat))));

% inside_0 = besseli(0,kappa_0,1)./(2*pi*besseli(0,kappa_x,1).*besseli(0,kappa_x,1)).*circ_vmpdf(mu_0,0,kappa_0);
% inside_1 = besseli(0,kappa_1,1)./(2*pi*besseli(0,kappa_x_1,1).*besseli(0,kappa_x_1,1)).*circ_vmpdf(mu_1,0,kappa_1);

insideexp_0 = sum(-log(besseli(0,kappa_x,1).*besseli(0,kappa_y,1)),2)+sum(kappa_0.*cos(mu_0),2);
insideexp_1 = sum(-log(besseli(0,kappa_x,1).*besseli(0,kappa_y,1)),2)+sum(kappa_1.*cos(mu_1),2);

numerator = logsumexp(insideexp_1,3);
denom = logsumexp(insideexp_0,3);

d = numerator - denom;

p = lambda/2 + (1-lambda)./(1+exp(beta0+beta.*d));
resp = (rand(nTrials, 1) < p);

