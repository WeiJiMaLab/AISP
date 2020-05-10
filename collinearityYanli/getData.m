function data = getData()

dat = load('data/newsubjdata.mat');

data = zeros(0, 6);
for i = 1:length(dat.newsubjdataC)
    datSub = dat.newsubjdataC{i};
    data = cat(1, data, [i * ones(size(datSub(:,1))), datSub(:, [1,3,4,6,7])]);
end