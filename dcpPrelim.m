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
    if strcmp(subject,'ar')
        datapath2 = ['/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle/' data];
    elseif strcmp(subject,'re')
        datapath2 = data;
    end
    if str2num(dataShort(1:end-1)) > 191114
        plxfile = [subject dataShort '.pl2'];
        kkfile = [subject dataShort '.kwik'];
        kkpath = [datapath2(1:end-1) 'kk'];
        kkflg = exist([kkpath '/' kkfile],'file');
    else
        kkflg = false;
        plxfile = [subject dataShort '.plx'];
    end
    if strcmp(subject,'ar')
        plxpath = ['/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle/' data(1:end-1) 'plx/'];
    end
%    plxpath = [datapath(1:end-1) 'plx'];

    dcp{listi} = dcpObj(subject,datapath);
    if extractSpikes && ~strcmp(FileList{listi},'20191021a') && ~strcmp(FileList{listi},'20191022a') && ~strcmp(FileList{listi},'20191022b')
        if kkflg
            dcp{listi} = extractKKData(dcp{listi},kkfile,kkpath,plxfile,plxpath,dcp{listi}.datapath);
        else
            dcp{listi} = extractSpikingData(dcp{listi},plxfile,plxpath,dcp{listi}.datapath);
        end
        dcp{listi} = unitsIndex(dcp{listi});                  % Finds the indices to the units
        dcp{listi} = tableImport(dcp{listi});
    end
end