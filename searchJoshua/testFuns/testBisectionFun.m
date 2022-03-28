function testBisectionFun()
% Runs the bisection function and produces some plots that can be used to
% see if it is working

for i = 1 : 3
    
    % Draw some random vals
    nItems = [2 : 6]';
    logKappa_x = (rand(5, 1) * 3)-0.5;
    kappa_x = exp(logKappa_x);
    kappa_s = [0, 1.5];
    
    % Run the bisection
    aisp_computeOptimalPointEstOfset(nItems, kappa_x, kappa_s, 0, ...
        true, false, true);
    
end