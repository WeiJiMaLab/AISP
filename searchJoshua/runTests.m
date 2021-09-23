function runTests()

addpath('./circstat-matlab')

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

