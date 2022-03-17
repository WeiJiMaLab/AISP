function noisy = addNoise(noiseless, kappas, varargin)
% Add VM noise to the values in the matrix noiseless. Values in
% noiseless that are nan remain as nan. Values are mapped back into the 
% range -pi to pi

% INPUT
% noiseless: matrix (may be 3D)
% kappas: vector as long as the first dimension of noiseless, and 
%   determines the kappa of the noise used for each row, or for 3D 
%   noiseless, determines the noise uses for [i, :, :] for each i. 
% varargin{1}: str. 'standard' to use qrandvm to add noise, 'efficientSamp'
%   to use an approach based on sampling importance resampling. Several
%   other options that do very random/nonsensical things.

% JCT, 2021

if isempty(varargin)
    mode = 'standard';
else
    mode = varargin{1};
end

assert(ndims(kappas)==2)
assert(size(kappas, 2)==1)
assert(size(kappas, 1)==size(noiseless, 1))

nlShape = size(noiseless);

if any(strcmp(mode, {'standard', 'useNormDist', 'randInt', 'impResamp'}))
    repFactors = [1, nlShape(2:end)];
    kappaShaped = repmat(kappas, repFactors);
    assert(isequal(size(kappaShaped), nlShape))
    
    inactive = isnan(noiseless);
    kappaShaped(inactive) = nan;
    
    if strcmp(mode, 'standard')
        addedNoise = qrandvm(0, kappaShaped, nlShape);
    elseif strcmp(mode, 'useNormDist')
        addedNoise = kappaShaped .* randn(nlShape);
    elseif strcmp(mode, 'randInt')
        addedNoise = randi(6000, nlShape);
    elseif strcmp(mode, 'impResamp')
        error('Option removed')
    else
        error('Uknknown setting')
    end

elseif strcmp(mode, 'efficientSamp')
    active = ~isnan(noiseless);
    uniKappas = unique(kappas);
    addedNoise = nan(nlShape);
    
    for iKappa = 1 : length(uniKappas)
        thisKappa = uniKappas(iKappa);
        match = thisKappa == kappas;
        
        thisActive = false(nlShape);
        thisActive(match, :, :) = active(match, :, :);
        numReq = sum(thisActive, 'all');
        samples = sampVm(0, thisKappa, numReq);
        addedNoise(thisActive) = samples(:);
    end
else
    disp('Setting is...')
    disp(mode)
    error('This setting is unknown')
end

noisy = addedNoise + noiseless;
noisy = vS_mapBackInRange(noisy, -pi, pi);

assert(isequal(isnan(addedNoise), isnan(noiseless)))
assert(isequal(isnan(noisy), isnan(noiseless)))

