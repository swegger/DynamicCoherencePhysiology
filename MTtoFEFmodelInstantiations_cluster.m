function MTtoFEFmodelInstantiations_cluster(simi,varargin)
%% gainNoiseNeuralModelParameterSweeps
%
%   [ws,sigGs,Gs] = gainNoiseNeuralModelParameterSweeps()
%
%   Sweeps through a default set of parameters and measures the resulting w
%   and sigG of circuit model output for comparison to behavior.
%
%%

%% Defaults
saveOpts_default.On = true;
saveOpts_default.location = ...
    ['model_' datestr(now,'yyyymmdd') '_0'];
saveOpts_default.Figs = false;


%% Parse inputs
Parser = inputParser;

addRequired(Parser,'simi')
addParameter(Parser,'fefN',1000)
addParameter(Parser,'leakRate',NaN)
addParameter(Parser,'fefBaseline',10)
addParameter(Parser,'randWeight',10)
addParameter(Parser,'structuredWeight',0)
addParameter(Parser,'fractTuned',0)
addParameter(Parser,'saveOpts',saveOpts_default)

parse(Parser,simi,varargin{:})

simi = Parser.Results.simi;
fefN = Parser.Results.fefN;
leakRate = Parser.Results.leakRate;
fefBaseline = Parser.Results.fefBaseline;
randWeight = Parser.Results.randWeight;
structuredWeight = Parser.Results.structuredWeight;
fractTuned = Parser.Results.fractTuned;
saveOpts = Parser.Results.saveOpts;


%% Run this instantiation

[modelRval,mtWeights,gainPrediction,spref] = simpleMTtoFEFmodel(...
    'fefN',fefN,'leakRate',leakRate,'fefBaseline',fefBaseline,...
    'randWeight',randWeight,'structuredWeight',structuredWeight,...
    'fractTuned',fractTuned,'pltFlg',false);

%% Save results
if saveOpts.On
    save([saveOpts.location(1:end-1) num2str(simi)])
end