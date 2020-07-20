function x = vS_mapBackInRange(x, lowerLim, upperLim)
% Maps the variables in the array x, back into the range lowerLim, and
% upperLim, assuming that they are circular variables and that lowerLim and
% upperLim are at the same location on the circle.

if lowerLim >= upperLim; error('incorrect use of inputs'); end


x = mod(x - lowerLim, upperLim - lowerLim) + lowerLim;