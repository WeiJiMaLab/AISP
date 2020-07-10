function d = aisp_computeOptimalPointEstDV(percept, nItems, kappa_x, kappa_s, mu_s)
% Compute the decision variable for the optimal point estimate observer

% First compute the standard point estiamte
d = aisp_computePointEstDV(percept, nItems, kappa_x, kappa_s, mu_s);

% Then calculate the optimal offset
aisp_computeOptimalPointEstOfset();