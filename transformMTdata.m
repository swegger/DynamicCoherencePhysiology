function data = transformMTdata(file,sname,varargin)
%% transformMTdata
%
%
%
%%

%% Defaults
saveOpts_default.On = false;

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'file')
addRequired(Parser,'sname')
addParameter(Parser,'fileDirectory','~/Projects/DynamicCoherencePhysiology/ar/mtdata')
addParameter(Parser,'saveOpts',saveOpts_default)

parse(Parser,file,sname,varargin{:})


file = Parser.Results.file;
sname = Parser.Results.sname;
fileDirectory = Parser.Results.fileDirectory;
saveOpts = Parser.Results.saveOpts;

%% Load data into matlab
mtdata = load(file);

%% For each neuron create an mtobj
neuronIDs = unique(mtdata.sethdata.neuron);
clear mtdata

for neuroni = 1:length(neuronIDs)
    mt{neuroni} = mtObj(sname,...
        [fileDirectory '/' file],...
        char(neuronIDs(neuroni)));
end

%% Save to destination
if saveOpts.On
    save([saveOpts.name '_' datestr(now,'yyyymmdd')],'mt')
end