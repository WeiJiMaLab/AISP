function resp = sample_simulate2(theta,dMat,logflag,tempp)
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
%
% ============ OUTPUT VARIABLES ============
% RESP: length nTrials vector of simulated responses

% persistent k_range
% persistent J_lin
% persistent highest_J

if nargin < 3; logflag = []; end

theta(logflag) = exp(theta(logflag));

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
% if isempty(k_range)
%     tempp = load('cdf_table.mat');
%     % K_interp = tempp.K_interp;
%     % cdf = tempp.cdf;
%     k_range = tempp.k_range;
%     J_lin = tempp.J_lin;
%     highest_J = tempp.highest_J;
% end
k_range = tempp.k_range;
J_lin = tempp.J_lin;
highest_J = tempp.highest_J;

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
xi = 1/diff(J_lin(1:2))*J_x_mat+1;
kappa_x = k_range(round(xi));
xi = 1/diff(J_lin(1:2))*J_y_mat+1;
kappa_y = k_range(round(xi));

if size(kappa_x,2) ~= nItems
    kappa_x = kappa_x';
    kappa_y = kappa_y';
end


% observed variables x, y
% instead of sampling the true orientations randomly
% we always set them to 1
xi = rand(nTrials, nItems) * 2 * pi - pi;
x = circ_vmrnd(xi, kappa_x);
y = circ_vmrnd(xi + Delta, kappa_y);

% replicate noise matrices based on number of samples
kappa_x = repmat(kappa_x, 1, 1, nSamples);
kappa_y = repmat(kappa_y, 1, 1, nSamples);

% sample first display and second display stimuli for each C
x_sample = circ_vmrnd(repmat(x, 1, 1, nSamples), kappa_x);
s1_DeltaMat = zeros(nItems, nTrials, nSamples);
DeltaVec = (rand(nTrials, 1, nSamples) .* 2 .* pi) - pi; % all Deltas
idx = randi(4, [nTrials * nSamples, 1]) ...
    + (0:4:(nItems * nTrials * nSamples - nItems))';
s1_DeltaMat(idx) = s1_DeltaMat(idx) + DeltaVec(:);
s1_DeltaMat = permute(s1_DeltaMat,[2 1 3]);
s1_DeltaMat = s1_DeltaMat + circ_vmrnd(repmat(x, 1, 1, nSamples), kappa_x);

d = logsumexp(sum(kappa_y .* cos(y - s1_DeltaMat), 2), 3)...
    -logsumexp(sum(kappa_y .* cos(y - x_sample), 2), 3);

p = lambda/2 + (1 - lambda)./(1 + exp(beta0 + beta .* d));
resp = (rand(nTrials, 1) < p);

