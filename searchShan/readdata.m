function [stimulus, response, performance] = readdata(subjid, half)
% input argument half specifies which half of the data should
% be returned, 0 means returning the whole data set, 1 means
% returning the half with odd indices, 2 means returning the half
% of data with even indices

if ~exist('half','var')
    half = 0;
end
if ~ismember(half, [0,1,2])
    error('Input of half is not valid, please input 0,1,2.')
end
if isnumeric(subjid)
    subs = {'MBC','MG','RC','WYZ','XLM','YC','YL','YMH','YZ'};
    subs2 = {'AK', 'GG','SJ','TQ','XZ'};
    subs = [subs,subs2];
    subjid = subs{subjid};
end

% get list of filenames
if ismember(subjid, {'MBC','MG','RC','WYZ','XLM','YC','YL','YMH','YZ'})
    datapath = 'data/Experiment 1/';
else
    datapath = 'data/Experiment 2/';
end
files = dir([datapath subjid '*.mat']);

% read every file
stimulus = [];
response = [];
performance = [];
for ii=1:length(files)
    load([datapath files(ii).name]);
    stimulus = [stimulus; data(:,1:2)];
    response = [response; data(:,3)];
    performance = [performance; data(:,4)];
end
idx = length(stimulus);
if half == 1
    stimulus = stimulus(1:2:idx-1,:);
    response = response(1:2:idx-1,:);
    performance = performance(1:2:idx-1,:);
elseif half == 2
    stimulus = stimulus(2:2:idx,:);
    response = response(2:2:idx,:);
    performance = performance(2:2:idx,:);
end

end