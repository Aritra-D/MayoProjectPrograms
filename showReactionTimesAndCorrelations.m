% Plot Reaction Times
zscoreFlag=1;

% folderSourceString = '/Volumes/Mojave/Users/supratimray/Supratim/Projects/MayoProject/';
folderSourceString = 'C:\Users\Supratim Ray\OneDrive - Indian Institute of Science\Supratim\Projects\Surya_MayoProject';
folderNameOut = fullfile(folderSourceString,'Data','segmentedData','behavior');
fileNameStringList=getAttentionExperimentDetails;
fileNameStringListAll=cat(2,fileNameStringList{1},fileNameStringList{2});

%%%%%%%%%%%%% Get the reaction times for all 12 conditions %%%%%%%%%%%%%%%%
oriChangeList = [2 3];
numConditions=12;
numSessions = length(fileNameStringListAll);

for i=1:numSessions
    
    fileNameSave=fullfile(folderNameOut,[fileNameStringListAll{i} '_reactionTimes.mat']);
    load(fileNameSave);
    
    for k=1:numConditions
        x = orientationChangeDeg(goodIndexList{k});
        goodPos = [];
        for j=1:length(oriChangeList)
            goodPos = cat(2,goodPos,find(x==uniqueOrientationChangeDeg(oriChangeList(j))));
        end
        goodPos = sort(goodPos);
        
        allReactionTimes{i,k} = reactionTimeMS(goodIndexList{k}(goodPos));
        xx = x(goodPos);
        
        orientationPos=zeros(1,length(xx));
        for j=1:length(oriChangeList)
            orientationPos(xx==uniqueOrientationChangeDeg(oriChangeList(j)))=j;
        end
        allOrientationIDs{i,k}=orientationPos;
    end
end

%%%%%%%%%%%%%%%%%%%%%% Classification of conditions %%%%%%%%%%%%%%%%%%%%%%%
indexList{1} = [1 2]; % Valid
indexList{2} = [9 10]; % Neutral
indexList{3} = [3 4]; % Invalid
    
numConditionsType = length(indexList);

for i=1:numConditionsType
    
    allRTs=[];
    for j=1:length(indexList{i})
        for s=1:numSessions
            allRTs = cat(2,allRTs,allReactionTimes{s,indexList{i}(j)});
        end
    end
    allRTVIN{i} = allRTs;
    allmeanRTVNI(i) = mean(allRTs);
    allsemRTVNI(i) = std(allRTs)/sqrt(length(allRTs));
end
h=subplot(221); 
errorbar(allmeanRTVNI,allsemRTVNI,'ko');
hold on;
plot(1:3,allmeanRTVNI,'color','k');
set(h,'Xtick',1:3,'Xticklabel',[{'V'} {'N'} {'I'}]);
xlim([0.5 3.5]);
ylabel('Reaction time (ms)'); ylim([300 350]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Projections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorNames = 'brg'; %spike, gamma and alpha
numMeasures=3;
measureNames = [{'Spikes'} {'Gamma'} {'Alpha'}];

subplot(222);
numFolds=5;
for i=1:numMeasures
    %[~,~,allProjections{i}]=displayHvsMAnalysisPopulation2(i,3,0,15,1,1,0,[],colorNames(i));
    [~,~,allProjections{i}]=displayHvsMAnalysisPopulation2(i,3,0,15,1,numFolds,0,[],colorNames(i));
    hold on;
end
ylabel('d-prime');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indexListCorrelation = [1:4 9 10];
indexListCorrelationNames = [{'H0V'} {'H1V'} {'H0I'} {'H1I'} {'H0N'} {'H1N'}];
numConditionCorrelations = length(indexListCorrelation);

pCutoff = 0.05/6;
for i=1:numMeasures % For each measure
    
    for j=1:numConditionCorrelations
        allRTs=[]; allProjs=[];
        
        allCorr2 = zeros(1,numSessions);
        for s=1:numSessions
            if zscoreFlag
                xTMP = zscore(allProjections{i}{s,indexListCorrelation(j)}');
                yTMP = zscore(allReactionTimes{s,indexListCorrelation(j)});
            else
                xTMP = (allProjections{i}{s,indexListCorrelation(j)}');
                yTMP = (allReactionTimes{s,indexListCorrelation(j)});
            end
            allProjs = cat(2,allProjs,xTMP);
            allRTs = cat(2,allRTs,yTMP);
            r=corrcoef(xTMP,yTMP);
            allCorr2(s) = r(1,2);
        end
        
        mCorr2(j) = mean(allCorr2);
        sCorr2(j) = std(allCorr2)/sqrt(numSessions);
        [~,allPVals2(j)] = ttest(allCorr2,0);
        
        [r,p] = corrcoef(allProjs,allRTs);
        allCorrs(j) = r(1,2);
        allPVals(j) = p(1,2);
        
        allRTsFull{i,j} = allRTs;
        allProjsFull{i,j} = allProjs;
    end
    
    h3=subplot(223);
    plot(allCorrs,'color',colorNames(i),'marker','o'); hold on;
    plot(allCorrs,'color',colorNames(i)); hold on;
    sigpts = find(allPVals<pCutoff);
    plot(sigpts,allCorrs(sigpts),'color',colorNames(i),'marker','o','markerFaceColor',colorNames(i),'markerSize',10,'linestyle','none'); hold on;
    ylabel('Correlation Coefficient');
    plot(1:numConditionCorrelations,zeros(1,numConditionCorrelations),'k');
    
    h4=subplot(224);
    errorbar(mCorr2,sCorr2,[colorNames(i) 'o']); hold on;
    plot(mCorr2,colorNames(i));
    sigpts2 = find(allPVals2<pCutoff);
    plot(sigpts2,mCorr2(sigpts2),colorNames(i),'marker','o','markersize',10,'markerFaceColor',colorNames(i),'linestyle','none');
    plot(1:numConditionCorrelations,zeros(1,numConditionCorrelations),'k');
    
%     h4=subplot(224);
%     plot(log10(allPVals),'color',colorNames(i)); hold on;
%     plot(1:numConditionCorrelations,log10(pCutoff)+zeros(1,numConditionCorrelations),'k');
%     ylabel('log10(pValue)');

end

set(h3,'XTick',1:numConditionCorrelations,'XTickLabel',indexListCorrelationNames);
set(h4,'XTick',1:numConditionCorrelations,'XTickLabel',indexListCorrelationNames);

figure;
for i=1:numMeasures
    if zscoreFlag==1
        allRTsVNI{1} = [zscore(allRTsFull{i,1}) zscore(allRTsFull{i,2})];
        allProjsVNI{1} = [-zscore(allProjsFull{i,1}) zscore(allProjsFull{i,2})];
        
        allRTsVNI{2} = [zscore(allRTsFull{i,5}) zscore(allRTsFull{i,6})];
        allProjsVNI{2} = [-zscore(allProjsFull{i,5}) zscore(allProjsFull{i,6})];
        
        allRTsVNI{3} = [(allRTsFull{i,3}) (allRTsFull{i,4})];
        allProjsVNI{3} = [-(allProjsFull{i,3}) (allProjsFull{i,4})];
    else
        allRTsVNI{1} = [(allRTsFull{i,1}) (allRTsFull{i,2})];
        allProjsVNI{1} = [-(allProjsFull{i,1}) (allProjsFull{i,2})];
        
        allRTsVNI{2} = [(allRTsFull{i,5}) (allRTsFull{i,6})];
        allProjsVNI{2} = [-(allProjsFull{i,5}) (allProjsFull{i,6})];
        
        allRTsVNI{3} = [(allRTsFull{i,3}) (allRTsFull{i,4})];
        allProjsVNI{3} = [-(allProjsFull{i,3}) (allProjsFull{i,4})];
    end
    
    for j=1:3
        subplot(3,3,3*(i-1)+j);
        plot(allRTsVNI{j},allProjsVNI{j},[colorNames(i) '.']); axis([-3 3 -3 3]);
        [r,p] = corrcoef(allRTsVNI{j},allProjsVNI{j});
        title(['r=' num2str(r(1,2),2) ', p=' num2str(p(1,2),2)]);
    end
        
end