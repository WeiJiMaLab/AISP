function outData = aisp_loadBestFits(dataDir, parsDir, ...
    model, varargin)
% Load information on the parameter fits for one model, find the best
% fits, and put into a cell array of structs.

% INPUT
% varargin{1}: 'array' | 'cellOfStruct'. If 'array' return a 
%   [num participants x num params] array. If 'cellOfStruct' return a 
%   num participant long cell array, with each element being a structure
%   representing the parameters for that participant. Default is 
%   'cellOfStruct'.

if (length(varargin)>0) && (~isempty(varargin{1}))
    outMode = varargin{1};
else
    outMode = 'cellOfStruct';
end

[~, Nptpnts] = getData(dataDir);

fname = [parsDir '/pars_' model '.mat'];
f = load(fname);

bestPars = aisp_collectBestFittingParams(f.nLogLs, f.pars);

if strcmp(outMode, 'cellOfStruct')
    allParamStructs = cell(Nptpnts, 1);
    for iP = 1 : Nptpnts
        allParamStructs{iP} = paramVec2Struct(bestPars(iP, :), model, ...
            'to struct');
    end
    outData = allParamStructs;
    
elseif strcmp(outMode, 'array')
    outData = bestPars;
    
else
    error('Unknown option selected')
end
    