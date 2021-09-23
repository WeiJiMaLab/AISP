function debugRun()

Dirs = produceDirs(5);
[DSet, Nptpnts] = getData(Dirs.dataDir);

Config = produceConfig();

iRep = 1;
iPtpnt = 1;
iType = 5;
idx = 1;
fit_cluster_ibs(iRep, iPtpnt, Config.ModelList{iType}, DSet, idx)