function percepts = aisp_addNoiseToStim(kappa_x, stim, varargin)
% Add noise to observed stimulus. The level of noise added depends on the 
% set size (the number of stimuli in the display). The code treats all 
% orieations with the same first index, ie,
% stim(i, :, :) as coming from the same trial.

% INPUT
% kappa_x: [num trials X 1] vector describing the value of kappa_x in each
%   trial, or single value
% stim: [num trials X num locations] or 
%   [num trials X num locations X num repitions] array describing the 
%   orientation of the Gabor patches
% varargin{1}: Passed to addNoise as varargin{1} for that function

% JCT

if ~isempty(varargin)
    passOn = varargin{1};
else
    passOn = [];
end 

% If kappa_x is a single value expland so it becomes a vector with one entry per
% trial
if all(size(kappa_x) == [1, 1])
    kappa_x = repmat(kappa_x, size(stim, 1), 1);
end

if isempty(passOn)
    percepts = addNoise(stim, kappa_x);
else
    percepts = addNoise(stim, kappa_x, passOn);
end
