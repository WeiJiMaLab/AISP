function timeAddNoise(varargin)

nTrials = 2000;
noiseless = rand(nTrials, 6, 200);
noiseless(noiseless > 0.5) = nan;
kappas = randi(6, [nTrials, 1]);

nReps = 20;
tic
for iRep = 1 : nReps
   addNoise(noiseless, kappas, varargin{1});
end
totalTime = toc;

disp(['Average duration of call to addNoise: ' num2str(totalTime/nReps)])

