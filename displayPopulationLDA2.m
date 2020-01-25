% This function saves the LDA or Regularized LDA projections of 12 cue conditions for different population size / features and Gamma
% Parameter

% dataType: 1:Firing Rate; 2:Gamma Power; 3: Alpha Power; 4:Firing Rate + Gamma Power; 5:Firing Rate + Gamma Power + Alpha Power   
% regFlag: 0: LDA ;      1: Regularized LDA ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similar to displayPopulationLDA but written semi-independently for cross
% verification. Only projections for H0V and H1V are used here.

% We also implement several ways of doing regularization. regFlag = 0 (no
% regularization), 1 (regularize only gamma), 2 (regularize both gamma and
% delta), 3 (use hypoeroptimzation)

function displayPopulationLDA2(dataType,regFlag,numFolds,useEqualStimRepsFlag,displayFlag,colorName)

if ~exist('regFlag','var');             regFlag=0;                      end
if ~exist('numFolds','var');            numFolds=5;                     end
if ~exist('useEqualStimRepsFlag','var'); useEqualStimRepsFlag=0;        end
if ~exist('displayFlag','var');         displayFlag=1;                  end
if ~exist('colorName','var');           colorName='b';                  end

oriChangeList = [2 3]; tpStr = '_TargetOnset'; timePeriod = [-0.5 0]; populationType = 'Stimulated'; 
tapers = [2 3]; alphaRangeHz = [8 12]; gammaRangeHz = [42 78];

folderSourceString = 'C:\Users\Supratim Ray\OneDrive - Indian Institute of Science\Supratim\Projects\Surya_MayoProject';
%folderSourceString = 'E:\Mayo';
folderNameData = fullfile(folderSourceString,'Data','savedDataSummary');
folderNameSave = fullfile(folderSourceString,'Data','populationSizeLDA');
fileNameStr = 'dataOri_';

for i=1:length(oriChangeList)
    fileNameStr = cat(2,fileNameStr,num2str(oriChangeList(i)));
end

fileNameStr = cat(2,fileNameStr,tpStr,num2str(timePeriod(1)),'_',num2str(timePeriod(2)),...
    '_tapers', num2str(tapers(1)),num2str(tapers(2)),'_alpha',num2str(alphaRangeHz(1)),'_',num2str(alphaRangeHz(2)),...
    '_gamma',num2str(gammaRangeHz(1)),'_',num2str(gammaRangeHz(2)));

fileNameData = fullfile(folderNameData,[fileNameStr '.mat']);
load(fileNameData); %#ok<*LOAD>

fileNameSaveElectrodes = fullfile(folderNameData,['electrodeArrayList' populationType '.mat']);
load(fileNameSaveElectrodes);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eList = cellfun(@(x) cat(1,x{1},x{2}),electrodeArrayList,'un',0);

if dataType==1
    data{1} = cellfun(@(x) x./diff(timePeriod),spikeData,'un',0); % converting spike count data into firing rate
    dataStr = 'FR';
elseif dataType==2
    data{1} = gammaData;
    dataStr = 'Gamma';
elseif dataType==3
    data{1} = alphaData;
    dataStr = 'Alpha';
elseif dataType==4
    data{1} = cellfun(@(x) x./diff(timePeriod),spikeData,'un',0); % converting spike count data into firing rate
    data{2} = gammaData;
    dataStr = 'FRGamma';
elseif dataType==5
    data{1} = cellfun(@(x) x./diff(timePeriod),spikeData,'un',0); % converting spike count data into firing rate
    data{2} = gammaData;
    data{3} = alphaData;
    dataStr = 'FRGammaAlpha';
end

fileNameStr2 = cat(2,dataStr,'_reg',num2str(regFlag),'_','Folds',num2str(numFolds),'equalReps',num2str(useEqualStimRepsFlag),tpStr,num2str(timePeriod(1)), '_',num2str(timePeriod(2)),...
    '_tapers', num2str(tapers(1)),num2str(tapers(2)),'_alpha',num2str(alphaRangeHz(1)),'_',num2str(alphaRangeHz(2)),...
    '_gamma',num2str(gammaRangeHz(1)),'_',num2str(gammaRangeHz(2)));
fileNameSave=fullfile(folderNameSave,[fileNameStr2,'.mat']);

if exist(fileNameSave,'file')
    load(fileNameSave);
else
    dPrimeList = getProjection(data,eList,regFlag,numFolds,useEqualStimRepsFlag);
    makeDirectory(folderNameSave);
    save(fileNameSave,'dPrimeList');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if displayFlag
    numSessions = length(dPrimeList);
    
    numElectrodeList = zeros(1,numSessions);
    for s=1:numSessions
        numElectrodeList(s) = length(dPrimeList{s});
    end
    maxElectrodes = max(numElectrodeList);
    
    numSessionPerN = zeros(1,maxElectrodes);
    mData = zeros(1,maxElectrodes);
    sData = zeros(1,maxElectrodes);
    
    for i=1:maxElectrodes
        allVals=[];
        for s=1:numSessions
            if i<=numElectrodeList(s)
                allVals = cat(2,allVals,dPrimeList{s}(i));
            end
        end
        numSessionPerN(i) = length(allVals);
        mData(i) = mean(allVals);
        sData(i) = std(allVals)/sqrt(numSessionPerN(i));
    end
    
    numSessionCutoff = 15;
    plotIndices = (numSessionPerN>=numSessionCutoff);
    errorbar(mData(plotIndices),sData(plotIndices),'color',colorName,'marker','o');
    hold on;
    plot(mData(plotIndices),'color',colorName);
    
end

end


function dPrimeList = getProjection(data,eList,regFlag,numFolds,useEqualStimRepsFlag)

[nSessions,~,~] = size(data{1});
dPrimeList = cell(1,nSessions);

for s=1:nSessions
    disp(['Working on session ' num2str(s) ' of ' num2str(nSessions)]);
    nSet = length(eList{s});

    eSet = getElectrodeSet(eList{s});
    N1 = length(data{1}{s,1,1}); % Number of stimulus repeats for HOV
    N2 = length(data{1}{s,2,1}); % Number of stimulus repeats for H1V
    [testingIndices1,testingIndices2] = getIndices(N1,N2,numFolds,useEqualStimRepsFlag);

    dPrimeTMP1 = zeros(1,nSet);
    for poolSize=1:nSet
        disp(['PoolSize ' num2str(poolSize) ' of ' num2str(nSet)]);

        numSetSize = size(eSet{poolSize},1);
        dPrimeTMP2 = zeros(1,numSetSize);
        for nIter= 1:numSetSize
            data1 = [];   data2 = [];
            for p=1:length(data)
                data1 = cat(1,data1,cell2mat(squeeze(data{p}(s,1,eSet{poolSize}(nIter,:))))); %H0V
                data2 = cat(1,data2,cell2mat(squeeze(data{p}(s,2,eSet{poolSize}(nIter,:))))); %H1V
            end
            
            [~,p1,p2,allIndices1,allIndices2] = getProjectionsAndWeightVector(data1,data2,testingIndices1,testingIndices2,3,regFlag);
            dPrimeTMP2(nIter) = getDPrime(p1(allIndices1),p2(allIndices2));
        end
        dPrimeTMP1(poolSize) = mean(dPrimeTMP2);
    end
    dPrimeList{s} = dPrimeTMP1;
end
end

function dPrime = getDPrime(p1,p2)
dPrime = (mean(p1)-mean(p2))/std(p2(:));
end
function eSet = getElectrodeSet(eList)
nSet = length(eList);

for i=1:nSet
    if i==1
        eSet{i} = eList; %#ok<*AGROW>
    elseif i==nSet
        eSet{i} = eList';
    else
        n=1;
        while n<=nSet
            eSet{i}(n,:) = sort(randsample(eList,i));
            eSet{i} = unique(eSet{i},'rows','stable');
            n=size(eSet{i},1)+1;
        end
    end
end
end
function [testingIndices1,testingIndices2] = getIndices(N1,N2,numFolds,useEqualStimRepsFlag)

allIndices1 = randperm(N1);
allIndices2 = randperm(N2);

if useEqualStimRepsFlag
    N1=min(N1,N2);
    N2=N1;
end

allIndices1 = sort(allIndices1(1:N1));
allIndices2 = sort(allIndices2(1:N2));
allIndices = [allIndices1 allIndices2];

testingIndices1 = cell(1,numFolds);
testingIndices2 = cell(1,numFolds);

if numFolds==1 % No cross Validation
    testingIndices1{1} = allIndices1;
    testingIndices2{1} = allIndices2;
else
    Y = [zeros(N1,1) ; ones(N2,1)];
    cvp = cvpartition(Y,'KFold',numFolds); % Partition data
    
    for i=1:numFolds
        testingIDs = find(cvp.test(i)==1);
        testingIndices1{i} = allIndices(testingIDs(testingIDs<=N1));
        testingIndices2{i} = allIndices(testingIDs(testingIDs>N1));
    end
end
end
function [weightVector,projections1,projections2,fullSetIndices1,fullSetIndices2] = getProjectionsAndWeightVector(data1,data2,testingIndices1,testingIndices2,transformType,regFlag)

numFolds = length(testingIndices1);
fullSetIndices1=[]; fullSetIndices2=[];
for i=1:numFolds
    fullSetIndices1 = cat(2,fullSetIndices1,testingIndices1{i});
    fullSetIndices2 = cat(2,fullSetIndices2,testingIndices2{i});
end
fullSetIndices1 = sort(fullSetIndices1);
fullSetIndices2 = sort(fullSetIndices2);

projections1 = zeros(size(data1,2),1);
projections2 = zeros(size(data2,2),1);
D = size(data1,1);

weightVectorTMP = zeros(D,numFolds);

for i=1:numFolds
    t1 = testingIndices1{i};
    t2 = testingIndices2{i};
    
    if numFolds==1 % No cross validation. Train and test on the same data
        train1 = t1;
        train2 = t2;
    else
        train1 = setdiff(fullSetIndices1,t1);
        train2 = setdiff(fullSetIndices2,t2);
    end
    d1 = data1(:,train1); d2 = data2(:,train2);
    m1 = size(d1,2); m2 = size(d2,2);
    
    if transformType==2 % LDA for uncorrelated case
        meanDiff = mean(d1,2) - mean(d2,2);
        var1 = var(d1,[],2);
        var2 = var(d2,[],2);
        pooledVar = ((m1-1)*var1 + (m2-1)*var2)/(m1+m2-2); %Pooled Variance
        
        weightVectorTMP(:,i) = meanDiff./pooledVar;
        
    elseif transformType==3 % LDA
        
        label1 = repmat({'H0V'},m1,1);
        label2 = repmat({'H1V'},m2,1);
        labelAll = cat(1,label1,label2);
        dataAll = cat(2,d1,d2);
        
        if regFlag==0 % No regularization
            Mdl = fitcdiscr(dataAll',labelAll);
            
        elseif regFlag==1 % Only optimize gamma
            Mdl = fitcdiscr(dataAll',labelAll);
            [err,gamma,~,~] = cvshrink(Mdl,'NumGamma',20);
            [~,minIndex] = min(err);
            if minIndex>1
                Mdl.Gamma = gamma(minIndex); % Changing gamma changes the model weights also
            end
            
         elseif regFlag==2 % Optimize gamma and delta
             Mdl = fitcdiscr(dataAll',labelAll);
             [err,gamma,delta,~] = cvshrink(Mdl,'NumGamma',20,'numDelta',20);
             minerr = min(err(:));
             [x,y] = find(err==minerr);
             if x(1)>1 || y(1)>1
                Mdl.Gamma = gamma(x(1)); % Take the smallest gamma
                Mdl.Delta = delta(x(1),y(1)); % Take the smallest delta
             end
             
        elseif regFlag==3 % Hyper optimize gamma and delta
            rng(1)
            myOpts.AcquisitionFunctionName = 'expected-improvement-plus';
            myOpts.ShowPlots = 0;
            myOpts.Verbose = 0;
            Mdl = fitcdiscr(X,Y,'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',myOpts);
            
        end
        weightVectorTMP(:,i) = Mdl.Coeffs(1,2).Linear;
    end
    
    projections1(t1) = data1(:,t1)' *  weightVectorTMP(:,i);
    projections2(t2) = data2(:,t2)' *  weightVectorTMP(:,i);
end

weightVector = mean(weightVectorTMP,2);
end