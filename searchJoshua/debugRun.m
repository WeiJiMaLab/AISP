function debugRun(modelNum)

addReqPaths()

Dirs = produceDirs(6);
[DSet, Nptpnts] = getData(Dirs.dataDir);

Config = produceConfig();

iRep = 1;
iPtpnt = 1;
idx = 1;
fit_cluster_ibs(iRep, iPtpnt, Config.ModelList{modelNum}, DSet, idx, true)