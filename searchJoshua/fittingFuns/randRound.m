function newVal = randRound(value)
% Rounds a value up on down, with a probability determined by 
% (one minus) the distance to the nearest integer

% INPUT
% value: scalar.

if floor(value) == ceil(value)
    newVal = value;
else
    distDown = value - floor(value);
    distUp = ceil(value) - value;

    if ~(round(distDown + distUp, 10) == 1)
        error('Bug')
    end

    u = rand(1);
    if u > distDown
        newVal = floor(value);
    else
        newVal = ceil(value);
    end
end
   
