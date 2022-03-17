function checkDerivations()
% Checks certain lines from the mathematcial derivations in the paper

% Specify number of repititions for each check
nReps = [1000, 1000, 1000];

% Key functions
Compute.Kappa_d = @(kappa_x, kappa_s, X) ...
    sqrt((kappa_x.^2) + (kappa_s.^2) + (2*kappa_x.*kappa_s.*cos(X)));
Compute.Rho = @(kappa_d, kappa_x, kappa_s) ...
    besseli(0, kappa_d) ./ (2*pi*besseli(0, kappa_x).*besseli(0, kappa_s));
Compute.NonTargTerm = @(X, S, kappa_x, kappa_s) ...
    circ_vmpdf(X', S', kappa_x)'.*circ_vmpdf(S', 0, kappa_s)';


%% CHECK 1
% Checking the steps from the first expression for B (eq. 126) to the
% expression for the maximum of B over s (eq. 136)

maxLogB = nan(nReps(1), 1);
expectedMaxLogB = nan(nReps(1), 1);
for iRep = 1 : nReps(1)
    
    % Randomly draw some values
    [nStim, X, kappa_x, kappa_s] = randomlyDrawVals();
    
    % Evalute the maximum of the expression for B in eq. 126 by hand
    maxLogB(iRep) = findByHandMaxB(X, kappa_x, kappa_s, Compute);
    
    % Evaluate the expression derived in the paper for this maximum
    expectedMaxLogB(iRep) = findByDerivationMaxB(X, kappa_x, kappa_s, Compute);
    
    if mod(iRep, 100) == 0
        disp('100 done')
    end
    
    if round(maxLogB(iRep), 8) > round(expectedMaxLogB(iRep), 8)
        disp('Unexpected result')
    end
end

figure
scatter(expectedMaxLogB, maxLogB);
refline(1, 0);


%% CHECK 2
% Going from the expression for A when L = i (eq. 138) to the expression for the
% maximum of A over s (eq. 142). Not that in this case, s_i (s at the target
% location) must equal to zero so don't have to maximise over this

maxLogA = nan(nReps(2), 1);
expectedMaxLogA = nan(nReps(2), 1);
for iRep = 1 : nReps(2)
    
    % Randomly draw some vals
    [nStim, X, kappa_x, kappa_s] = randomlyDrawVals();
    L = randi(nStim, [1, 1]);
    X_T = X(L);
    X_nonT = X;
    X_nonT(L) = [];
    
    % Evalute the maximum of the expression for A in eq. 138 by hand
    maxLogA(iRep) = findByHandMaxA(X_T, X_nonT, kappa_x, kappa_s, Compute);
    
    % Evaluate the expression derived in the paper for this maximum
    expectedMaxLogA(iRep) = findByDerivationMaxA(X_T, X_nonT, kappa_x, kappa_s, Compute);
    
    if mod(iRep, 100) == 0
        disp('100 done')
    end
    
    if round(maxLogA(iRep), 8) > round(expectedMaxLogA(iRep), 8)
        disp('Unexpected result')
    end
end

figure
scatter(expectedMaxLogA, maxLogA);
refline(1, 0);


%% CHECK 3
% Compute through maximisations by hand, the point estimate DV. Compare
% this to the point estimate DV computed using the derivations

dv_byHand = nan(nReps(3), 1);
dv_byDeriv = nan(nReps(3), 1);
for iRep = 1 : nReps(3)
    
    % Randomly draw some vals
    [nStim, X, kappa_x, kappa_s] = randomlyDrawVals();
    
    % Evalute the maximum of the expression for B in eq. 126 by hand
    maxLogB = findByHandMaxB(X, kappa_x, kappa_s, Compute);
    maxB = exp(maxLogB);
    
    % For each possible target location, evaluate the maximum of the
    % expression for A in eq. 138 by hand
    maxAtThisL = nan(nStim, 1);
    for L = 1 : nStim

        X_T = X(L);
        X_nonT = X;
        X_nonT(L) = [];
        
        % Evalute
        maxAtThisL(L) = findByHandMaxA(X_T, X_nonT, kappa_x, kappa_s, Compute);
    end
    
    % Find the maximum value of A over locations
    maxA = max(exp(maxAtThisL));
    
    dv_byHand(iRep) = log( (maxA / nStim) / maxB);
    
    
    % Also evaluate the DV using the derived expressions
    dv_byDeriv(iRep) = aisp_computePointEstDV(X, nStim, kappa_x, ...
        kappa_s, 0, 'stimAndTarg', true);
    
    if mod(iRep, 100) == 0
        disp('100 done')
    end
end


figure
scatter(dv_byHand, dv_byDeriv);
refline(1, 0);

end


function [nStim, X, kappa_x, kappa_s] = randomlyDrawVals

nStim = randi(5, [1, 1]) +1;
X = rand(nStim, 1)'*2*pi;
logKappa_x = (rand(1, 1) * 3)-0.5;
kappa_x = exp(logKappa_x);
possKappa_s = [0, 1.5];
kappa_s = possKappa_s(randi(2, 1, 1));

end


function vmTerm = computeVmPdfZeroZero(kappa_d)
% Compute the PDF of the von Misses distribution at value = 0, mean = 0,
% for different values of kappa_d

vmTerm = nan(1, length(kappa_d));
for iS = 1 : length(kappa_d)
    vmTerm(iS) = circ_vmpdf(0, 0, kappa_d(iS));
end

end


function maxLogA = findByHandMaxA(X_T, X_nonT, kappa_x, kappa_s, Compute)
% Evalute the maximum of the expression for A in eq. 138 by hand

logA = @(X_T, X_nonT, S_nonT) ...
    sum(log(Compute.NonTargTerm(X_nonT, S_nonT, kappa_x, kappa_s))) ...
    + log(circ_vmpdf(X_T, 0, kappa_x));

minS = repmat(-pi, 1, length(X_nonT));
maxS = repmat(pi, 1, length(X_nonT));
startS = repmat(0, 1, length(X_nonT));
[paramVals, minNegLogA] = fmincon(@(S_nonT) -logA(X_T, X_nonT, S_nonT), ...
    startS, [], [], [], [], minS, maxS);
maxLogA = -minNegLogA;

end


function maxLogA = findByDerivationMaxA(X_T, X_nonT, kappa_x, kappa_s, Compute)
% Evaluate the expression derived in the paper for this maximum

kappa_d = Compute.Kappa_d(kappa_x, kappa_s, X_nonT);
rhoTerm = Compute.Rho(kappa_d, kappa_x, kappa_s);
vmTerm = computeVmPdfZeroZero(kappa_d);

maxLogA = log(circ_vmpdf(X_T, 0, kappa_x)) + ...
    sum(log(rhoTerm.*vmTerm));

end


function maxLogB = findByHandMaxB(X, kappa_x, kappa_s, Compute)
% Evalute the maximum of the expression for B in eq. 126 by hand

logB = @(X, S) sum(log(Compute.NonTargTerm(X, S, kappa_x, kappa_s)));

minS = repmat(-pi, 1, length(X));
maxS = repmat(pi, 1, length(X));
startS = repmat(0, 1, length(X));
[paramVals, minNegLogB] = fmincon(@(S) -logB(X, S), startS, [], [], [], [], minS, maxS);
maxLogB = -minNegLogB;

end

function maxLogB = findByDerivationMaxB(X, kappa_x, kappa_s, Compute)
% Evaluate the expression derived in the paper for this maximum

mu_d = X + atan2(-sin(X), (kappa_x / kappa_s) + cos(X));
kappa_d = Compute.Kappa_d(kappa_x, kappa_s, X);
rhoTerm = Compute.Rho(kappa_d, kappa_x, kappa_s);
vmTerm = computeVmPdfZeroZero(kappa_d);

maxLogB = sum(log(rhoTerm.*vmTerm));

end



