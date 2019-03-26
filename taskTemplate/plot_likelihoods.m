function plot_likelihoods(penalty)
%function plot_likelihoods(penalty)
% uses the result of evaluate to generate a plot of the performance of the
% different models allowing for a switch "penalty" which allows to switch
% between raw likelihood(0), AIC(1) and BIC(1) 
% for consistency accross plots this should probably only be writting once!

if ~exist('penalty','var') || isempty(penalty)
    penalty=1;
end