function ibsExample

responses = ones(50, 1);

for nReps = [1, 2, 3, 10, 50]
    opt_ibs = ibslike;
    opt_ibs.Nreps = nReps;
    opt_ibs.Vectorized = true;
    
    value = ibslike(@(pars, data)randomResponses(pars, data), ...
        0, responses, [], opt_ibs);
    disp(['Estimated nLogL with ' num2str(nReps) ' Nreps: ' num2str(value)])
end

end


function resp = randomResponses(~, data)
    resp = randi(2, size(data, 1), 1);
end

