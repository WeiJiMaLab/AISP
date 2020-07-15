function [ofset, assocNItems, assocKappa_x, assocKappa_s] ...
    = aisp_computeOptimalPointEstOfset(nItems, kappa_x, kappa_s, mu_s)
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
if any(size(nItems) ~= size(kappa_x)); error('These are meant to correspond.'); end
if all(unique(nItems) ~= nItems); error('See note'); end
if all(unique(kappa_x) ~= kappa_x); error('See note'); end
if all(unique(kappa_s) ~= kappa_s); error('See note'); end

% OUTPUT
% ofset: array, giving the optimal offset for each combination of nItems and kappa_s
% assocNItems: Array of same size as ofset indicating the corresponding nItems
% for each entry in ofset
% assocKappa_s: Same as assocNItems but for kappa_s

persistent optimalCrit
persistent savedNItems
persistent savedKappa_x
persistent savedKappa_s

nSim = 500000;

% Only recompute the optimal criteria if these have changed
if isempty(savedNItems) || isempty(savedKappa_x) || isempty(savedKappa_s) || ...
        (savedNItems ~= nItems) || ...
        (savedKappa_x ~= kappa_x) || (savedKappa_s ~= kappa_s)
    
    [assocNItems, assocKappa_s] = meshgrid(nItems, kappa_s);
    assocKappa_x = nan(size(assocNItems)); 
    optimalCrit = nan(size(assocNItems));
    
    for iCase = 1 : length(optimalCrit(:))
        thisNItems = assocNItems(iCase);
        thisKappa_x = kappa_x(nItems == thisNItems);
        assocKappa_x(iCase) = thisKappa_x;
        thisKappa_s = assocKappa_s(iCase);
        
        % Simulate stimuli with and without target
        s0 = qrandvm(mu_s, thisKappa_s, [nSim, thisNItems]);
        s1 = nan(nSim, thisNItems);
        s1(:, 1 : (thisNItems-1)) ...
            = qrandvm(mu_s, thisKappa_s, [nSim, (thisNItems-1)]);
        s1(:, end) = repmat(mu_s, nSim, 1); % Adding the target
        
        % Simulate a value of the associated decision variable, for a point
        % estimate observer
        x0 = aisp_addNoiseToStim(thisKappa_x, s0);
        d0 = aisp_computePointEstDV(x0, thisNItems, thisKappa_x, ...
            thisKappa_s, mu_s);
        
        x1 = aisp_addNoiseToStim(thisKappa_x, s1);
        d1 = aisp_computePointEstDV(x1, thisNItems, thisKappa_x, ...
            thisKappa_s, mu_s);
        
        if (length(size(d0))~=2) || (size(d0, 2)~=1); error('Bug'); end
        if (length(size(d1))~=2) || (size(d1, 2)~=1); error('Bug'); end
        
        % bisection method
        crit0 = min(d1);
        crit2 = max(d0);
        crit1 = crit0 + (crit2-crit0) / ( 3 + sqrt(5) ) * 2;
        f1 = sum(d0 < crit1) + sum(d1 > crit1);
        while (crit2-crit0) > 0.001
            if abs(crit2-crit1) > abs(crit1-crit0)
                crit_new = crit1 + (crit2-crit1) / ( 3 + sqrt(5) ) * 2;
                f_new = sum(d0 < crit_new) + sum(d1 > crit_new);
                if f_new < f1
                    crit2 = crit_new;
                else
                    crit0 = crit1;
                    f1 = f_new;
                    crit1 = crit_new;
                end
            else
                crit_new = crit1 + (crit0-crit1) / ( 3 + sqrt(5) ) * 2;
                f_new = sum(d0 < crit_new) + sum(d1 > crit_new);
                if f_new < f1
                    crit0 = crit_new;
                else
                    crit2 = crit1;
                    f1 = f_new;
                    crit1 = crit_new;
                end
            end
        end
        
        % If the distributions were already completely seperated, then we
        % would never have entered the while loop, and something potentially
        % strange is going on.
        if crit0 > crit2
            warning('Complete seperation of the distributions')
            % Use midpoint between ends of the two distributions
            crit1 = (crit0 + crit2)/2;
        end
        
        optimalCrit(iN, iK) = crit1;
    end
    
    % Update persistant variables to reflect that have recomputed the ofset
    savedNItems = nItems;
    savedKappa_x = kappa_x;
    savedKappa_s = kappa_s;
end


ofset = optimalCrit;

