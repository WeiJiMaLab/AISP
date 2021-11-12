function samples = sampVm(mu, kappa, nSamples, varargin)
% Using an eclectic mix of methods, including sampling methoods, try to 
% produce (approximate) von Mises samples as fast as possible.
%
% When kappa = 0 use rand
% When mu = 0 and kappa = 1.5 use either (a) samples from a large pool of 
%   predrawn von Mises samples, or (b) sampling importance resampling
% Otherwise use sampling importance resampling

% INPUT
% mu: scalar. Centre of the distribution
% kappa: scalar. Concentration parameter
% nSamples: scalar. How many samples to draw?
% varargin{1}: bool. If true, for the special case of mu = 0 and 
%   kappa = 1.5 use approach (a) described above, otherwise use approach 
%   (b). Default is true

% JCT, 2021

if len(varargin) >= 1
    fromPool = varargin{1};
    assert(isboolean(fromPool))
else
    fromPool = true;
end

assert(length(mu) == 1)
assert(length(kappa) == 1)
isSpecialCase = kappa == 1.5 && mu == 0;

persistent evenSpaced
persistent kappaSpecialProbs
persistent kappaSpecialPool

if isempty(evenSpaced)
    evenSpaced = -pi : 0.001 : pi;
end

if isSpecialCase
    if fromPool && isempty(kappaSpecialPool)
        kappaSpecialPool = qrandvm(mu, kappa, 50000);
    elseif (~fromPool) && isempty(kappaSpecialProbs)
        kappaSpecialProbs = circ_vmpdf(evenSpaced, mu, kappa);
    end
end

if kappa == 0
    samples = (rand(nSamples, 1) * 2 * pi) - pi;
else
    % WORKING HERE --- use the kappaSpecialPool
    if isSpecialCase
        assocProbs = kappaSpecialProbs;
    else
        assocProbs = circ_vmpdf(evenSpaced, mu, kappa);
    end
    samples = randsample(evenSpaced, nSamples, true, assocProbs);
end


