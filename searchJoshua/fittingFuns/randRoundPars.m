function parVec = randRoundPars(parVec, pos)
% Rounds a parameter up on down, with a probability determined by 
% (one minus) the distance to the nearest integer

% INPUT
% parVec: Vector or parameters
% pos: The index of the parameter to round up or down

relParam = parVec(pos);

distDown = relParam - floor(relParam);
distUp = ceil(relParam) - relParam;
assert(round(distDown + distUp, 10) == 1)

u = rand(1);
if u > distDown
    newVal = floor(relParam);
else
    newVal = ceil(relParam);
end

parVec(pos) = newVal;
   
