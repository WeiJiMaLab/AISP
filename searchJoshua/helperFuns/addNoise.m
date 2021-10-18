function noisy = addNoise(noiseless, kappas, varargin)
% Add VM noise to the values in the matrix noiseless. Values in
% noiseless that are nan remain as nan. Values are mapped back into the 
% range -2pi to pi

% INPUT
% noiseless: matrix (may be 3D)
% kappas: vector as long as the first dimension of noiseless, and 
%   determines the kappa of the noise used for each row. 
% varargin: str. Only provide when want to debug/time the function. Does
%   lots of random things.

% JCT, 2021

if ismepty(varargin)
    mode = 'standard';
else
    mode = varargin{1};
end

assert(ndims(kappas)==2)
assert(size(kappas, 2)==1)
assert(size(kappas, 1)==size(noiseless, 1))

nlShape = size(noiseless);
repFactors = [1, nlShape(2:end)];
kappaShaped = repmat(kappas, repFactors);
assert(isequal(size(kappaShaped), nlShape))

inactive = isnan(noiseless);
kappaShaped(inactive) = nan;

if strcmp(mode, 'standard')
    addedNoise = qrandvm(0, kappaShaped, nlShape);
elseif strcmp(mode, 'useNormDist')
    addedNoise = kappaShaped .* randn(nlShape);
end 

noisy = addedNoise + noiseless;
noisy = vS_mapBackInRange(noisy, -pi, pi);

assert(isequal(isnan(noisy), isnan(noiseless)))