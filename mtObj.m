%% mtObj
%
%   Defines properties and methods of MT neuron data used to compare to FEF
%   data.
%
%   Defined as a subclass of dcpObj
%
%%

classdef mtObj < dcpObj
    
    properties (SetAccess = private)
        objType = 'mtObj';
        neuronName;
        coh;
        r;
        R;
        Rste;
        neuron_t;
    end
    
    
    methods
        %% Core methods
        function obj = mtObj(sname,datapath,neuronName)
        % Constructor
            obj = obj@dcpObj(sname,datapath);
            obj.neuronName = neuronName;
            obj.unitIndex = 1;
        end
        
        %% mtObjTrials
        function obj = mtObjTrials(obj,trialNs)
            
            ind = 0;
            
            % Find matches in table
            mtdata = load(obj.datapath);
            fnames = fieldnames(mtdata);
            
            if any(strcmp(fnames,'sethdata'))
                matches = find(mtdata.sethdata.neuron == obj.neuronName);
            
            
                obj.neuron_t = (-400:100:800) + 50;
                
                neuronIDs = unique(mtdata.sethdata.neuron);
                
                for ti = trialNs
                    
                    if length(matches) > ti
                        ind = ind+1;
                        obj.coh(ind,1) = mtdata.sethdata.coh(matches(ti));
                        obj.speeds(ind,1) = mtdata.sethdata.speed(matches(ti));
                        if mtdata.sethdata.prefdir(matches(ti)) == 'true'
                            obj.directions(ind,1) = 0;
                        else
                            obj.directions(ind,1) = 180;
                        end
                        obj.r(:,ind) = str2num(mtdata.sethdata.bins{matches(ti)}(5:end-1))';
                    end
                    
                end
                
            elseif any(strcmp(fnames,'MTdataTable'))
                matches = find(strcmp(mtdata.MTdataTable.neuron,obj.neuronName));
                
                datal = nan(size(mtdata.MTdataTable,1),1);
                for i = 1:size(mtdata.MTdataTable,1)
                    datal(i) = length(mtdata.MTdataTable.Var7{1});
                end
                
                obj.neuron_t = mtdata.xaxis(1:min(datal));
            
                neuronIDs = unique(mtdata.MTdataTable.neuron);
                
                for ti = trialNs
                    
                    if length(matches) > ti && mtdata.MTdataTable.speedpulse{matches(ti)} == 0 && mtdata.MTdataTable.cohpulse{matches(ti)} == 0
                        ind = ind+1;
                        obj.coh(ind,1) = mtdata.MTdataTable.coherence(matches(ti));
                        obj.speeds(ind,1) = mtdata.MTdataTable.speed{matches(ti)};
                        if mtdata.MTdataTable.prefdir(matches(ti))
                            obj.directions(ind,1) = 0;
                        else
                            obj.directions(ind,1) = 180;
                        end
                        obj.r(:,ind) = mtdata.MTdataTable.Var7{matches(ti)}(1:min(datal))';
                    end
                    
                end
                
            else
                error('Fields not recognized!')
            end
            
        end
        
        
        %% Neural analysis methods
        
        function C = findCovariance(obj,binT,speeds,cohs,dirs)
            res = nan(length(binT),1000,length(obj.unitIndex));
            ind = 1;
            for di = 1:length(dirs)
                for si = 1:length(speeds)
                    for ci = 1:length(cohs)
                        [~,condLogical] = trialSort(obj,dirs(di),speeds(si),NaN,cohs(ci));
                        for bi = 1:length(binT)
                            counts = obj.r(binT(bi) == obj.neuron_t,condLogical,:);
                            res(bi,ind:ind+sum(condLogical)-1,:) = counts - mean(counts);
%                             if length(obj.unitIndex) == 1
%                                 resTemp(bi,:) = counts*obj.filterWidth*2 - mean(counts*obj.filterWidth*2);
%                             else
%                                 resTemp(bi,:,:) = counts*obj.filterWidth*2 - mean(counts*obj.filterWidth*2);
%                             end
                        end
%                         res(:,ind:ind+size(resTemp,2)-1,:) = resTemp;
                        ind = ind+sum(condLogical);
                    end
                end
            end
            res = res(:,1:ind-1,:);
            for uniti = 1:length(obj.unitIndex)
                C(:,:,uniti) = cov(res(:,:,uniti)');
            end
            
        end
        
        function varCE = findVarCE(obj,speeds,cohs,dirs)
            rtemp = obj.r*obj.filterWidth*2;
            res = nan(size(rtemp));
            res = permute(res,[1,3,2]);
            ind = 1;
            for di = 1:length(dirs)
                for si = 1:length(speeds)
                    for ci = 1:length(cohs)
                        [~,condLogical] = trialSort(obj,dirs(di),speeds(si),NaN,cohs(ci));
                        m(:,:,di,si,ci) = mean(rtemp(:,condLogical,:),2);
                        res(:,:,ind:ind+sum(condLogical)-1) = permute(rtemp(:,condLogical,:),[1,3,2]) - ...
                            repmat(m(:,:,di,si,ci),[1,1,sum(condLogical)]);
                        n(1,1,di,si,ci) = sum(condLogical);
                        ind = ind+sum(condLogical);
                    end
                end
            end
            res = res(:,:,1:ind-1);
            M = sum(repmat(n,[size(m,1),size(m,2),1,1,1]).*m/sum(n(:)),[3,4,5]);
            V = var(res,[],3);
            FF = V./M;
            phi = min(FF,[],1);
            varCE = V - repmat(phi,[size(M,1),1]).*M;
        end
    end
end