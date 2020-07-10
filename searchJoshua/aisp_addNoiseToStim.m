function percepts = aisp_addNoiseToStim(kappa_x, stim)
% Add noise to observed stimulus. The level of noise added depends on the 
% set size (the number of stimuli in the display). The code treats all 
% orieations with the same first index, ie,
% stim(i, :, :) as coming from the same trial.

% INPUT
% kappa_x: [num trials X 1] vector describing the value of kappa_x in each
% trial, or single value
% setSizeCond: [num trials X 1] vector describing the set size condition. This
% is a number between 1 and the total number of possible set sizes. It is not
% the set size itself.
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

% Produce array of simulated von Mises noise. Set up an array where rows are
% trials, and the third dimention indexes trial simulations (draws).
noise = NaN(size(stim));


% Only need to simulate noise for locations in which a stimulus was presented.
presentedLocs = ~isnan(stim);


% Find the relevant Kappa_x value. This depends on the set size.
relKappaXShaped = repmat(kappa_x, 1, size(stim, 2));

if length(size(stim)) == 3
    relKappaXShaped = repmat(relKappaXShaped, 1, 1, size(stim, 3));
end

relKappaXShaped(~presentedLocs) = NaN;

noise ...
    = qrandvm(0, relKappaXShaped, size(noise));


% Add to the stimuli and then map them back in range if they are outside
% [-pi pi]
percepts = stim + noise;

percepts = vS_mapBackInRange(percepts, -pi, pi);


% Check we have NaNs in the same place as we started
if any(any(any(isnan(stim) ~= isnan(percepts)))); error('bug'); end

