function data = getData()

dataFile = 'D:\Research data _ BACKED UP\Visual search\Main Study\In use\StandardFormat_participantExcluded.mat';

Loaded = load(dataFile);
data = Loaded.DSet;