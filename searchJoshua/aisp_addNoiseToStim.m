function percepts = aisp_addNoiseToStim(kappa_x, stim)
% Add noise to observed stimulus. The level of noise added depends on the 
% set size (the number of stimuli in the display). The code treats all 
% orieations with the same first index, ie,
% stim(i, :, :) as coming from the same trial.

% INPUT
% kappa_x: [num trials X 1] vector describing the value of kappa_x in each
% trial, or single value
% stim   [num trials X num locations] array describing the orientation of the
%        Gabor patches

% Joshua Calder-Travis 
% j.calder.travis@gmail.com
% GitHub: jCalderTravis

% If kappa_x is a single value expland so it becomes a vector with one entry per
% trial
if all(size(kappa_x) == [1, 1])
    kappa_x = repmat(kappa_x, size(stim, 1), 1);
end

percepts = addNoise(stim, kappa_x);
