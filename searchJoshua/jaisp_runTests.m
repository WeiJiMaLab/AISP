function jaisp_runTests(dataDir)

% INPUT
% dataDir: Directory containing the original unfitted dataset

% JCT, 2020 - 2022

%% Fitting with checks on

[DSet, ~] = getData(dataDir);
Config = produceConfig();

repNum = 1;
for iModel = 1 : length(Config.ModelList)
    thisModel = Config.ModelList{iModel};
    
    for thisPtpnt = [3, 9]
        fit_cluster_ibs(repNum, thisPtpnt, thisModel, DSet, repNum, ...
            false, true)
        repNum = repNum +1;
    end
end


%% imResampVm
mu = pi/8;
kappa = 1;
nSamples = 2000000; 

imResampSamples = imResampVm(mu, kappa, nSamples);
standardSamples = qrandvm(mu, kappa, nSamples);

figure; hold on
histogram(imResampSamples)
histogram(standardSamples)


%% vmpdf
x = -pi : 0.01 : pi;
mu = randn(size(x));
kappa = abs(randn(size(x))*10);

expected = nan(size(x));
for i = 1 : size(x, 1)
    for j = 1 : size(x, 2)
        expected(i, j) = circ_vmpdf(x(i, j), mu(i, j), kappa(i, j));
    end
end

obtained = vmpdf(x, mu, kappa);

figure
scatter(expected, obtained)
refline(1, 0)

error = abs(obtained - expected);
assert(all(error(:)<0.000001))
