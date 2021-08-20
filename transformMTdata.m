function data = transformMTdata(file,varargin)
%% transformMTdata
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'file')

parse(Parser,file,varargin{:})

file = Parser.Results.file;

%% for each neuron create an mtobj