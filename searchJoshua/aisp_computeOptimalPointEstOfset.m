function [ofset, finalAssocNItems, finalAssocKappa_x, finalAssocKappa_s] ...
    = aisp_computeOptimalPointEstOfset(nItems, kappa_x, kappa_s, mu_s, varargin)
% Compute the offset required to optimally compensate for changes in kappa_x and
% kappa_s

% INPUT
% nItems: Vector of unique values of number of items in the display
% kappa_x: Vector same length as nItems. Values of kappa_x for each value of
% nItems
% kappa_s: Vector of unique values of kappa_s to use
% varargin: If set to true, recompute the ofset regardless of whether the imput
% options have changed or not. (Normally does not recompute if inputs are the
% same.)

% NOTE
% nItems, kappa_x, kappa_s do not have to be unique but this would be a waste of
% computation time

% OUTPUT
% ofset: array, giving the optimal offset for each combination of nItems and kappa_s
% assocNItems: Array of same size as ofset indicating the corresponding nItems
% for each entry in ofset
% assocKappa_s: Same as assocNItems but for kappa_s

if ~isequal(unique(nItems), sort(nItems)); error('See note'); end
if any(size(nItems) ~= size(kappa_x)); error('These are meant to correspond.'); end
if ~isequal(unique(kappa_s), sort(kappa_s)); error('See note'); end

if ~isempty(varargin)
    forceRecompute = varargin{1};
else
    forceRecompute = false;
end

persistent optimalCrit
persistent savedNItems
persistent savedKappa_x
persistent savedKappa_s
persistent assocNItems
persistent assocKappa_x
persistent assocKappa_s

nSim = 50000; %500000;


% Work out if the requested calulations are the same as ones done previously,
% a subset, or if they are new
if forceRecompute
    calc = 'new';
elseif isempty(savedNItems) || isempty(savedKappa_x) || isempty(savedKappa_s)
    calc = 'new';
elseif isequal(savedNItems, nItems) ...
        && isequal(savedKappa_x, kappa_x) ...
        && isequal(savedKappa_s, kappa_s)
    calc = 'prev';
elseif all(ismember(nItems(:), savedNItems(:))) ...
        && all(ismember(kappa_x(:), savedKappa_x(:))) ...
        && all(ismember(kappa_s(:), savedKappa_s(:)))
    calc = 'subset';
else
    calc = 'new';
end

% Only recompute the optimal criteria if these have changed
if strcmp(calc, 'new')
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
            thisKappa_s, mu_s, 'stimAndTarg');
        
        x1 = aisp_addNoiseToStim(thisKappa_x, s1);
        d1 = aisp_computePointEstDV(x1, thisNItems, thisKappa_x, ...
            thisKappa_s, mu_s, 'stimAndTarg');
        
        if (length(size(d0))~=2) || (size(d0, 2)~=1); error('Bug'); end
        if (length(size(d1))~=2) || (size(d1, 2)~=1); error('Bug'); end
        
        % bisection method
        crit0 = min(d1);
        crit2 = max(d0);
        crit1 = crit0 + (crit2-crit0) / ( 3 + sqrt(5) ) * 2;
        %%
%         % For debugging
%         figure; hold on
%         histogram(d0, 'Normalization', 'countdensity')
%         histogram(d1, 'Normalization', 'countdensity')
%         plot([crit0, crit0], get(gca,'YLim'))
%         plot([crit2, crit2], get(gca,'YLim'))
%         %%%
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
%         %%%
%         % For debugging
%         plot([crit1, crit1], get(gca,'YLim'))
%         %%%
        
        % If the distributions were already completely seperated, then we
        % would never have entered the while loop, and something potentially
        % strange is going on.
        if crit0 > crit2
            warning('Complete seperation of the distributions')
            % Use midpoint between ends of the two distributions
            crit1 = (crit0 + crit2)/2;
        end
        
        optimalCrit(iCase) = crit1;
    end
    
    % Update persistant variables to reflect that have recomputed the ofset
    savedNItems = nItems;
    savedKappa_x = kappa_x;
    savedKappa_s = kappa_s;
end

finalAssocNItems = assocNItems;
finalAssocKappa_x = assocKappa_x;
finalAssocKappa_s = assocKappa_s;
ofset = optimalCrit;

