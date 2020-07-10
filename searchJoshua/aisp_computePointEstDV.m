function d = aisp_computePointEstDV(percept, nItems, kappa_x, kappa_s, mu_s)
% Compute the decision variable for the point estimate observer

% INPUT
% percept   [numTrials x setSize x numSamples] array of stimulus percepts
% nItems
% kappa_x   Observer's belief about the concetration parameter of the 
%           measurement noise
% kappa_s   Overserver's belief about the concentration parameter of the 
%           distractor distribution
% mu_s      Center of the distractor distribution

%% Setup
% Check input
if size(percept, 2) > 8; error('Bug'); end

% We use implicit expansion below so it is very important that all input
% vectors are the expected shape.
inputVectors = {nItems, kappa_s, kappa_x};

for iInputVec = 1 : length(inputVectors)
    vecSize = size(inputVectors{iInputVec});
    if (length(vecSize) ~= 2) || (vecSize(2) ~= 1)
        error('Bug')
    end
end

% If kappa_s is passed as a single value, expand it into a vector so that
% the calculations for 'termC' below, work.
if size(kappa_s, 1) == 1
    kappa_s = repmat(kappa_s, size(percept, 1), 1);
end

% Same for kappa_x
if size(kappa_x, 1) == 1
    kappa_x = repmat(kappa_x, size(percept, 1), 1); 
end


%% Main calculations
% TODO check maths behind this
error('Check maths')

% TODO not sure if the arguments to atan2 are in the correct order
cosPercept = cos(percept);
mu_d = percept + atan2(-sim(percept), (kappa_x./kappa_s)+cosPercept);
kappa_d = sqrt((kappa_x.^2) + (kappa_s.^2) + (2*kappa_x.*kappa_s.*cosPercept));

varphi = exp(kappa_d .* (cos(mu_d) -1));


% Find the number of elements, tildeL, in Omega, that maximise J (see
% derivations). From the derivations we know that if tildeL items are included
% they will be the items with the greatest values of varphi. 
varphi = sort(varphi, 2, 'descend');
varphiProd = cumprod(varphi, 2);
tildeL = nan(1, size(varphi, 2), 1);
tildeL(:) = 1 : size(varphi, 2);

% Compute tildeL * product of varphi, for all possibel numbers of tildeL, and
% pick the maximum value
maxProduct = max(varphiProd.*tildeL, [], 2);
if any(isnan(maxProduct(:))); error('Bug'); end

% TODO especially check maths behind these next few lines
if length(mu_s) == 1 && mu_s == 0
    logBesseli = aisp_computeLogBesseliForDuplicatedValues(kappa_s)
    logVmTerm = log(2*pi) - kappa_s + logBesseli;
    % TODO check that this matches a von mises evaluated at x=0, and with mean
    % mu=0. (Not sure if got the right bessel function so important to check.
else
    error('Not coded up yet')
end
d = log( (1./(2*pi*nItems)) .* maxProduct ) + logVmTerm;

% Check assumptions made in the derivations 
if any(mu_d(:)==0); error('Derivations assume all mu_d are non-zero'); end




