function out = dcpChamberMap(subject,varargin)
%% ChamberMapGUI
%
%   out = ChamberMapGUI(d)
%
%   Graphical user interphase for examining neural data for a project based
%   on location in the chamber.
%
%   Date    Initials    Comment
%   150520  swe         Initial commit
%%

% Defaults
ChamberDimensions_default.radius = 12;
ChamberDimensions_default.x = 0;
ChamberDimensions_default.y = 0;
ChamberDimensions_default.Properties.Color = [0 0 0];
ChamberDimensions_default.Properties.LineWidth = 2;
excludeList_default = {'20191111a','20191115a','20191118t','20191119a',...
    '20200210a','20200213a','20200218a','20200221a','20200309b'};

% Parse input
Parser = inputParser;

addRequired(Parser,'subject')
addParameter(Parser,'dcpObjectsFile',[])
addParameter(Parser,'PlotType','Raster')
addParameter(Parser,'CustomParameters',[]);
addParameter(Parser,'ChamberDimensions',ChamberDimensions_default);
addParameter(Parser,'UserSuppliedObjects',[]);
addParameter(Parser,'SaveLocation','default');
addParameter(Parser,'extractSpikes',true)
addParameter(Parser,'trList',1:5000)
addParameter(Parser,'boxCarWidth',30)
addParameter(Parser,'excludeList',excludeList_default)

parse(Parser,subject,varargin{:})

subject = Parser.Results.subject;
dcpObjectsFile = Parser.Results.dcpObjectsFile;
PlotType = Parser.Results.PlotType;
CustomParameters = Parser.Results.CustomParameters;
ChamberDimensions = Parser.Results.ChamberDimensions;
UserSuppliedObjects = Parser.Results.UserSuppliedObjects;
SaveLocation = Parser.Results.SaveLocation;
extractSpikes = Parser.Results.extractSpikes;
trList = Parser.Results.trList;
boxCarWidth = Parser.Results.boxCarWidth;
excludeList = Parser.Results.excludeList;

% Validate inputs
if strcmp('Custom',PlotType)
    if isempty(CustomParameters)
        error('Custom plot types require custom parameter inputs!')
    end
end

      
%% Initalize variables
UnitsString = {''};
FigureTypesForFile = {'ts v tp'};
FiguresForFile = {''};
axesForFigure = {''};
temphandle = [];
colors = [     0    0.4470    0.7410;...
          0.8500    0.3250    0.0980;...
          0.9290    0.6940    0.1250;...
          0.4940    0.1840    0.5560;...
          0.4660    0.6740    0.1880;...
          0.3010    0.7450    0.9330;...
          0.6350    0.0780    0.1840];
      
% Build FileList
potentialFiles = dir(['~/Projects/DynamicCoherencePhysiology/' subject]);
regOut = regexpi({potentialFiles.name},'[0-9]{8}[a-z]{1,3}','match');
ind = 0;
for listi = 1:length(regOut)
    if ~isempty(regOut{listi}) && ~strcmp(regOut{listi}{1}(end-2:end),'plx') ...
            && ~strcmp(regOut{listi}{1}(end-1:end),'kk') ...
            && ~any(strcmp(regOut{listi}{1},excludeList))
        ind = ind+1;
        FileList{ind} = regOut{listi}{1};
    end
end

%% Run preliminary analyses
if isempty(dcpObjectsFile)
    dcp = dcpPrelim(subject,FileList,extractSpikes);
else
    load(dcpObjectsFile,'dcp')
end

%% Generate GUI

% Set up main figure
h = figure;
h.Units = 'normalized';
h.Position = [0.05 0.25 0.90 0.50];
% vbox = uix.VBox( 'Parent', h );
% hbox1 = uix.HBox( 'Parent', h, 'Padding', 1 );

% Dynamic figure
% FigureAx = axes( 'Parent', hbox1, 'Title', 'DummyTitle' );

% Chamber Figure
% ChamberAx = axes( 'Parent', hbox1 ,'Title', ['Chamber map for ' subject]);
% hbox1.Widths = [-1 -1];

% Selection GUI
hbox2 = uix.HBox( 'Parent', h, 'Padding', 2 );
ChamberAx = axes( 'Parent', hbox2 ,'Title', ['Chamber map for ' subject]);
fileSelectionBox = uix.BoxPanel( 'Parent', hbox2, 'Title', 'Maestro file' );
unitSelectionBox = uix.BoxPanel( 'Parent', hbox2, 'Title', 'Units' );
analysisSelectionBox = uix.BoxPanel( 'Parent', hbox2, 'Title', 'Analyses' );
selectFile = uicontrol('Parent',fileSelectionBox,'Style','listbox','String',FileList,'Max',1,'Callback',@fileSelectionCallback);
selectUnit = uicontrol('Parent',unitSelectionBox,'Style','listbox','String',UnitsString,'Max',1,'Callback',@unitSelectionCallback);
selectAnalysis = uicontrol('Parent',analysisSelectionBox,'Style', 'listbox',...
    'String',{'dirPref','speedPref','initiateCoh','dynamicCoh'},...
    'Max',1,'Callback',@analysisCallback);
% vbox.Heights = [-2 -1]; 

%% Plot recording sites
axes(ChamberAx)
% plot(0,0,'k.','MarkerSize',10)
hold on
xhandle = plot3(NaN,NaN,NaN,'x');
xhandle.Visible = 'off';
for i = 1:length(dcp)
    if isempty(dcp{i}.location)
        siteHandles(i) = plot3(NaN,NaN,NaN,'k.');
        ringHandles(i) = plot3(NaN,NaN,NaN,'ko');
        set(ringHandles(i),'Visible','off');
        %sites(i,:) = [NaN,NaN,NaN];
    else
        siteHandles(i) = plot3(dcp{i}.location.x,dcp{i}.location.y,-dcp{i}.location.depth,'k.');
        ringHandles(i) = plot3(dcp{i}.location.x,dcp{i}.location.y,-dcp{i}.location.depth,'ko');
        set(ringHandles(i),'Visible','off');
        %sites(i,:) = [dcp{i}.location.x,dcp{i}.location.y,dcp{i}.location.depth];
    end
    figsOpen(i) = 0;
    
    % Plot the Chamber, if provided
    if ~isempty(ChamberDimensions)
        [tempx,tempy] = ellipse(ChamberDimensions.radius,ChamberDimensions.radius,...
            ChamberDimensions.x,ChamberDimensions.y,pi/360);
        temp = plot3(tempx,tempy,zeros(size(tempx)),'k');
        if isfield(ChamberDimensions,'Properties')
            set(temp,ChamberDimensions.Properties)
        end
    end
    
    % Plot other objects, if provided
    if ~isempty(UserSuppliedObjects)
        for j = 1:length(UserSuppliedObjects)
            temp = plot(UserSuppliedObjects{j}.x,UserSuppliedObjects{j}.y,'k');
            if isfield(UserSuppliedObjects{j},'Properties')
                set(temp,UserSuppliedObjects{j}.Properties)
            end
        end
    end
end
% axis equal
xlabel('Medial-lateral (mm)')
ylabel('Anterior-posterior (mm)')
zlabel('Depth from cortical surface (um)')
view([45,10])
grid on
% set(gca,'ButtonDownFcn', @mouseclick_callback)
disp('')

%% Callbacks

% File selection callback
    function fileSelected = fileSelectionCallback(gcbo,eventdata)
        % Get the index of the file selected
        index = selectFile.Value;
        
        % Turn on outer ring and set color of dot to be the same color
        for i = 1:length(dcp) %d.runs
            if any(index == i)
                colorindex = size(colors,1)+1;
                iter = 0;
                while any(colorindex > size(colors,1))
                    colorindex = index - iter*size(colors,1);
                    iter = iter+1;
                end
                siteHandles(i).Color = colors(colorindex,:);
                ringHandles(i).Visible = 'on';
                ringHandles(i).Color = colors(colorindex,:);
                
            else
                siteHandles(i).Color = [0 0 0];
                ringHandles(i).Visible = 'off';
                ringHandles(i).Color = [0 0 0];
            end
        end
        xhandle.Visible = 'off';
        fileSelected = FileList{selectFile.Value};
        
        % Generate unit selection callback
        UnitsString = {};
        if isempty(dcp{selectFile.Value}.unitIndex)
            UnitsString{1} = 'None';
        else
            totalunits = length(dcp{selectFile.Value}.unitIndex);
            for i = 1:totalunits
                UnitsString{i} = num2str(dcp{selectFile.Value}.unitIndex(i));
            end
        end
        selectUnit.String = UnitsString;
        selectUnit.Value = 1;
        
    end

% Unit selection callback
    function unitSelected = unitSelectionCallback(gcbo,eventdata)
        
        unitSelected = UnitsString{selectUnit.Value};
        index = selectFile.Value;
        
        if length(dcp{index}.location.x) > 1
            siteIndex = floor(dcp{index}.chansIndex(selectUnit.Value)/4)+1;
            
            % Place x on correct location based on siteIndex
            tempIndex = find(~isnan(siteHandles(index).XData));
            xhandle.XData = siteHandles(index).XData(tempIndex(siteIndex));
            xhandle.YData = siteHandles(index).YData(tempIndex(siteIndex));            
            xhandle.ZData = siteHandles(index).ZData(tempIndex(siteIndex));
            xhandle.Color = siteHandles(index).Color;     
            xhandle.Visible = 'on';
        end
        
        % Generate list of analyses
        selectAnalysis.String = {'dirPref','speedPref','initiateCoh','dynamicCoh'};
    end

% Analysis selection callback
    function analysisType = analysisCallback(gcbo,eventdata)
        % Generate temporary vboxes
        
        analysisType = selectAnalysis.String{selectAnalysis.Value};
        
        switch analysisType
            case 'dirPref'
                dirPref = dirPrefObj(subject,...
                    dcp{selectFile.Value}.datapath);         % Builds direction preference object
                dirPref = assertSpikesExtracted(dirPref,...
                    dcp{selectFile.Value}.spikesExtracted);  % Asserts that spiking data has (not) been extracted already
                dirPref = unitsIndex(dirPref);                  % Finds the indices to the units
                dirPref = dirPrefTrials(dirPref,trList);        % Finds dirPref trial data
                
                % Sort by direction and plot rasters
                unitSelected= str2double(UnitsString{selectUnit.Value});
                dirPrefRaster(dirPref,0:45:315,unitSelected);
                
                % Plot grand average PSTH across all trials
                figure;
                r = dirPref.calcRates(boxCarWidth);                      % Calculates rates based with 50 ms boxcar
                plot(-100:1600,mean(r(:,:,selectUnit.Value)*1000,2));
                hold on
                plotVertical(0);
                xlabel('Time since motion onset (ms)')
                ylabel('Sp/s')
                mymakeaxis(gca);
                
            case 'speedPref'
                speedPref = speedPrefObj(subject,...
                    dcp{selectFile.Value}.datapath);
                speedPref = assertSpikesExtracted(speedPref,...
                    dcp{selectFile.Value}.spikesExtracted);
                speedPref = unitsIndex(speedPref);
                speedPref = speedPrefTrials(speedPref,trList);
                
                % Sort by speed and plot rasters
                unitSelected= str2double(UnitsString{selectUnit.Value});
                speedPrefRaster(speedPref,unique(speedPref.directions),...
                    unique(speedPref.speeds),unitSelected)
                
                % Plot grand average PSTH across all trials
                figure;
                r = speedPref.calcRates(boxCarWidth);                      % Calculates rates based with 50 ms boxcar
                plot(-100:1600,mean(r(:,:,selectUnit.Value)*1000,2));
                hold on
                plotVertical(0);
                xlabel('Time since motion onset (ms)')
                ylabel('Sp/s')
                mymakeaxis(gca);
                
            case 'initiateCoh'
                initCoh = initiateCohObj(subject,...
                    dcp{selectFile.Value}.datapath);
                initCoh = assertSpikesExtracted(initCoh,...
                    dcp{selectFile.Value}.spikesExtracted);
                initCoh = unitsIndex(initCoh);
                initCoh = initiateCohTrials(initCoh,trList);
                dirs = [0];
                
                % Plot behavioral data
                initiateCohMeanEye(initCoh,dirs);
                
                % Sort by coherence and plot rates
                if ~isempty(initCoh.coh)
                    [R,Rste] = cohConditionedRates(initCoh,'width',boxCarWidth,'dirs',dirs);
                    
                    figure('Name','Speed/coherence conditioned rates','Position',[342 449 1972 420])
                    t = -100:1600;
                    for cohi = 1:length(unique(initCoh.coh))
                        subplot(1,length(unique(initCoh.coh)),cohi)
                        %     for speedi = 1:length(unique(initCoh.speeds))
                        plot(t(t<300),R(t<300,:,cohi,selectUnit.Value)*1000)
                        %     end
                        hold on
                        axis tight
                        ax(cohi,:) = axis;
                    end
                    
                    for cohi = 1:length(unique(initCoh.coh))
                        subplot(1,length(unique(initCoh.coh)),cohi)
                        axis([min(ax(:,1)) max(ax(:,2)) min(ax(:,3)) max(ax(:,4))])
                        plotVertical(0);
                        xlabel('Time (ms)')
                        ylabel('sp/s')
                        mymakeaxis(gca)
                    end
                    
                    figure('Name','Speed/coherence conditioned rates 2','Position',[342 449 1972 420])
                    t = -100:1600;
                    for speedi = 1:length(unique(initCoh.speeds))
                        subplot(1,length(unique(initCoh.speeds)),speedi)
                        %     for speedi = 1:length(unique(initCoh.speeds))
                        plot(t(t<300),squeeze(R(t<300,speedi,:,selectUnit.Value)*1000))
                        %     end
                        hold on
                        axis tight
                        ax(speedi,:) = axis;
                    end
                    
                    for speedi = 1:length(unique(initCoh.speeds))
                        subplot(1,length(unique(initCoh.speeds)),speedi)
                        axis([min(ax(:,1)) max(ax(:,2)) min(ax(:,3)) max(ax(:,4))])
                        plotVertical(0);
                        xlabel('Time (ms)')
                        ylabel('sp/s')
                        mymakeaxis(gca)
                    end
                end
                
            case 'dynamicCoh'
                dynamicCoh = dynamicCohObj(subject,...
                    dcp{selectFile.Value}.datapath);
                dynamicCoh = assertSpikesExtracted(dynamicCoh,...
                    dcp{selectFile.Value}.spikesExtracted);
                dynamicCoh = unitsIndex(dynamicCoh);
                dynamicCoh = dynamicCohTrials(dynamicCoh,trList);
                dirs = [0];
                
                % Plot behavioral data
                dynamicCohMeanEyeSeq(dynamicCoh,dirs);
                ax = axis;
                
                % Sort by speed and plot rates
                t = -100:1600;
                if ~isempty(dynamicCoh.sequences)
                    [R,Rste] = dynamicCohSeqConditionedRates(dynamicCoh,'width',boxCarWidth,...
                        'dirs',dirs,'t',t);
                    
                    controlInd = 5;
                    figure;
                    subplot(2,1,1)
                    plot(t,R(:,:,selectUnit.Value))
                    hold on
                    plotVertical([150 150+0:300:1500]);
                    xlim(ax(1:2))
                    subplot(2,1,2)
                    plot(t,R(:,:,selectUnit.Value) - repmat(R(:,controlInd,selectUnit.Value),[1,size(R,2)]))
                    hold on
                    plotVertical([150 150+0:300:1500]);
                    xlim(ax(1:2))
                end
        end
    end




%% Set output on close of window
disp('Close GUI to continue...')
waitfor(h)
close(temphandle)
out = [];

end % End of ChamberMapGUI