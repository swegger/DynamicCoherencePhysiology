function dcp = dcpPrelim(subject,FileList,extractSpikes)
%% dcpPrelim
%
%
%%

%% Perform preliminary analysis
for listi = 1:length(FileList)
    data = FileList{listi};
    dataShort = data(3:end);
    datapath = ['/home/seth/Projects/DynamicCoherencePhysiology/' subject '/' data];
    if str2num(dataShort(1:end-1)) > 191114
        plxfile = [subject dataShort '.pl2'];
    else
        plxfile = [subject dataShort '.plx'];
    end
    plxpath = [datapath(1:end-1) 'plx'];

    dcp{listi} = dcpObj(subject,datapath);
    if extractSpikes && ~strcmp(FileList{listi},'20191021a') && ~strcmp(FileList{listi},'20191022a') && ~strcmp(FileList{listi},'20191022b')
        dcp{listi} = extractSpikingData(dcp{listi},plxfile,plxpath,dcp{listi}.datapath);
    end
    dcp{listi} = tableImport(dcp{listi});
    dcp{listi} = unitsIndex(dcp{listi});                  % Finds the indices to the units
end