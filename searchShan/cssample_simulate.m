function responses = cssample_simulate(stim, pars)
% function responses = bayes_simulate(stim, pars)
% simulates responses for the Bayes optimal observer
% parameter passed as pars are:
%       sigmaNoise = pars(1)
%       beta = pars(2)
%       bias = pars(3)
%       lapse = pars(4)
%       n_samp = pars(5)
%
% Note: In contrast to the original paper we use the log evidence ratio,
% not the difference in probability.

sigmaNoise = pars(1);
beta = pars(2);
bias = pars(3);
lapse = pars(4);
if length(pars) > 4
    n_samples = floor(pars(5));
    if n_samples > 10000
        n_samples = 10000;
        warning('n_samples reduced to 10000 to avoid extreme computation')
    elseif rand > mod(n_samples, 1)
        n_samples = floor(n_samples);
    else
        n_samples = ceil(n_samples);
    end
else
    n_samples = 1;
end

n_trials = size(stim, 1);

sigma = 9.06;

range = (1:n_trials)';

s = [stim(:,1), stim(:,2), stim(:,2), stim(:,2)];
x = s + sigmaNoise*randn(size(s));

% initialize sampling to prior samples
s_samp = sigma * randn(n_trials, 2);
C_samp = randi(2, n_trials, 1);
s_samp(:, 1) = (2*C_samp-3) .* abs(s_samp(:,1));
Loc_samp = randi(4, n_trials, 1);
% s_1 = nan(n_trials, n_samples, 2);
% s_2 = nan(n_trials, n_samples, 2);
C = nan(n_trials, n_samples);

idx = (Loc_samp-1) .* n_trials + range;

% evaluate assumed target position
l_samp = -(x(idx) - s_samp(:,1)) .^2 ./ sigmaNoise.^2 ./2;

% evaluate all other positions
for iL = 1:3
    idx = mod(Loc_samp-1 + iL, 4) .* n_trials + range;
    l_samp = l_samp -(x(idx) - s_samp(:,2)) .^2 ./ sigmaNoise.^2 ./2;
end


for i = 1:n_samples
    s_samp_new = sigma * randn(n_trials, 2);
    C_samp_new = randi(2, n_trials, 1);
    s_samp_new(:, 1) = (2*C_samp_new-3) .* abs(s_samp(:,1));
    Loc_samp_new = randi(4, n_trials, 1);
    idx = (Loc_samp_new-1) .* n_trials + range;
    % evaluate assumed target position
    l_samp_new = -(x(idx) - s_samp_new(:,1)) .^2 ./ sigmaNoise.^2 ./2;

    % evaluate all other positions
    for iL = 1:3
        idx = mod(Loc_samp_new-1 + iL, 4) .* n_trials + range;
        l_samp_new = l_samp_new -(x(idx) - s_samp_new(:,2)) .^2 ./ sigmaNoise.^2 ./2;
    end
    accept = rand(n_trials,1)< exp(l_samp_new-l_samp);
    % s1_samp(accept, :) = s1_samp_new(accept, :);
    l_samp(accept) = l_samp_new(accept);
    C_samp(accept) = C_samp_new(accept);
    %s_samp(accept) = s_samp_new(accept);
    
    C(:,i) = C_samp;
    % s(:,i,:) = s_samp;
end

d = log(sum(C==1, 2) + 1) - log(n_samples - sum(C==1, 2) + 1);
p = lapse/2 + (1-lapse)./(1+exp(bias+beta*d));
responses = (rand(n_trials, 1) < p);