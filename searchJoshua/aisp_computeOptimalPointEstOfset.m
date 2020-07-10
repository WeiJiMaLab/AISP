function ofset = aisp_computeOptimalPointEstOfset(nItems, kappa_x, kappa_s, mu_s)
% Compute the offset required to optimally compensate for changes in kappa_x and
% kappa_s

% INPUT
% nItems: Vector of unique values of number of items in the display
% kappa_x: Vector same length as nItems. Values of kappa_x for each value of
% nItems
% kappa_s: Vector of unique values of kappa_s to use

% NOTE
% nItems, kappa_x, kappa_s do not have to be unique but this would be a waste of
% computation time
if all(unique(nItems) ~= nItems); error('See note'); end
if all(unique(kappa_x) ~= kappa_x); error('See note'); end
if all(unique(kappa_s) ~= kappa_s); error('See note'); end

% OUTPUT
% ofset: [length(nItems), length(kappa_s)] array, giving the
% optimal offset for each combination of kappa_x and kappa_s


persistent optimalCrit
persistent savedNItems
persistent savedKappa_x
persistent savedKappa_s

nSim = 500000;

% Only recompute the optimal criteria if these have changed
if isempty(savedNItems) || isempty(savedKappa_x) || isempty(savedKappa_s) || ...
        (savedNItems ~= nItems) || ...
        (savedKappa_x ~= kappa_x) || (savedKappa_s ~= kappa_s)
    optimalCrit = nan(length(nItems), length(kappa_s));
    
    for iN = 1 : length(nItems)
       for iK = 1 : length(kappa_s)
           
           % Simulate stimuli with and without target
           s0 = qrandvm(mu_s, kappa_s(iK), [nSim, nItems(iN)]);
           s1 = qrandvm(mu_s, kappa_s(iK), [nSim, (nItems(iN)-1)]);
           s1 = [repmat(mu_s, nSim, 1), s1]; % Adding the target
           
           % Simulate a value of the associated decision variable, for a point
           % estimate observer
           x0 = aisp_addNoiseToStim(kappa_x(iN), s0);
           d0 = aisp_computePointEstDV(x0, nItems(iN), kappa_x(iN), ...
               kappa_s(iK), mu_s);
           
           x1 = aisp_addNoiseToStim(kappa_x(iN), s1);
           d1 = aisp_computePointEstDV(x1, nItems(iN), kappa_x(iN), ...
               kappa_s(iK), mu_s);
           
           
           % TODO
           % bisection method
           
           optimalCrit(iN, iK) = nan; %TODO result here
       end
    end
    
    % Update persistant variables to reflect that have recomputed the ofset
    savedNItems = nItems;
    savedKappa_x = kappa_x;
    savedKappa_s = kappa_s;
end


ofset = optimalCrit;