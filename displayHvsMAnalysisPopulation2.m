% This Program is used in displayFigure5, displayFigure5v2 to save the
% output of this function.
% The difference between this function and
% displayHvsMPopulation is that here pooled variance is used for uncorrelated LDA case

% dataTypeNum - 1: spikes, 2: gamma, 3: alpha
% transformType - 1: simple averaging across sides, 2: LDA-uncorrelated,
% 3: LDA-covariance

% Adding an option to do cross-validation. Set numFolds to 1 if you do not want to do cross-validation
% Also adding an option to equalize the stimulus repeats for H0V and H1V

function [mData,typeList1] = displayHvsMAnalysisPopulation2(dataTypeNum,transformType,trialCutoff,normalizeFlag,numFolds,useEqualStimRepsFlag)

if ~exist('transformType','var');       transformType=2;                end
if ~exist('trialCutoff','var');         trialCutoff=15;                 end
if ~exist('normalizeFlag','var');       normalizeFlag=1;                end
if ~exist('numFolds','var');            numFolds=5;                     end
if ~exist('useEqualStimRepsFlag','var'); useEqualStimRepsFlag=0;        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oriChangeList = [2 3]; tpStr = '_TargetOnset'; timePeriod = [-0.5 0]; populationType = 'Stimulated';
tapers = [2 3]; alphaRangeHz = [8 12]; gammaRangeHz = [42 78];

folderSourceString = 'C:\Users\Supratim Ray\OneDrive - Indian Institute of Science\Supratim\Projects\Surya_MayoProject';
% folderSourceString = 'E:\Mayo';
folderNameSave = fullfile(folderSourceString,'Data','savedDataSummary');
fileNameStr = 'dataOri_';
for i=1:length(oriChangeList)
    fileNameStr = cat(2,fileNameStr,num2str(oriChangeList(i)));
end

fileNameStr = cat(2,fileNameStr,tpStr,num2str(timePeriod(1)),'_',num2str(timePeriod(2)),...
    '_tapers', num2str(tapers(1)),num2str(tapers(2)),'_alpha',num2str(alphaRangeHz(1)),'_',num2str(alphaRangeHz(2)),...
    '_gamma',num2str(gammaRangeHz(1)),'_',num2str(gammaRangeHz(2)));

fileNameSave = fullfile(folderNameSave,[fileNameStr '.mat']);
load(fileNameSave); %#ok<*LOAD>

fileNameSaveElectrodes = fullfile(folderNameSave,['electrodeArrayList' populationType '.mat']);
load(fileNameSaveElectrodes);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

typeList1 = [{'H0V'} {'H1V'} {'H0I'} {'H1I'} {'M0V'} {'M1V'} {'M0I'} {'M1I'} {'H0N'} {'H1N'} {'M0N'} {'M1N'}];

colorNameList = 'rbg';
colorName = colorNameList(dataTypeNum);
if dataTypeNum==1
    dataTMP = transformData(spikeData,transformType,dataTypeNum,timePeriod,electrodeArrayList,numFolds,useEqualStimRepsFlag);
elseif dataTypeNum==2
    dataTMP = transformData(gammaData,transformType,dataTypeNum,timePeriod,electrodeArrayList,numFolds,useEqualStimRepsFlag);
elseif dataTypeNum==3
    dataTMP = transformData(alphaData,transformType,dataTypeNum,timePeriod,electrodeArrayList,numFolds,useEqualStimRepsFlag);
end

if normalizeFlag
    dataTMP = normalizeData(dataTMP);
end

[numSessions,numConditions] = size(dataTMP);
mData = cell(1,numConditions);
nTrials = zeros(numSessions,numConditions);

for c=1:numConditions
    pos=1;
    for s=1:numSessions
        tmp = dataTMP{s,c};
        nTrials(s,c) = length(tmp);
        if nTrials(s,c)>trialCutoff
            mData{c}(pos) = mean(tmp);
            pos = pos+1;
        else
            continue
        end
    end
end

% Plot Data
mDataTMP = zeros(1,numConditions);
conditionOrder = [4 1 11 10 3 2 12 9 8 5 7 6]; 

for c=1:numConditions
    
    x = mData{c};
    mtmp = mean(x);
    ntmp = length(x);
    stmp = std(x)/sqrt(ntmp);
    errorbar(conditionOrder(c),mtmp,stmp,[colorName 'o']); hold on;
    text(conditionOrder(c)-0.2,mtmp+stmp+0.005,num2str(ntmp));
    mDataTMP(conditionOrder(c)) = mtmp;
end
plot(mDataTMP,colorName);
h=gca;
set(h,'XTick',1:12,'XTickLabel',typeList1);

end

function newData = transformData(data,transformType,dataTypeNum,timePeriod,electrodeArrayList,numFolds,useEqualStimRepsFlag)

[numSessions,numConditions,~] = size(data);
newData=cell(numSessions,numConditions);

for s=1:numSessions
    eList = [electrodeArrayList{s}{1};electrodeArrayList{s}{2}]';
    if dataTypeNum==1
        data1 = cell2mat(squeeze(data(s,1,eList)))./diff(timePeriod); % H0V
        data2 = cell2mat(squeeze(data(s,2,eList)))./diff(timePeriod); % H1V
    else
        data1 = cell2mat(squeeze(data(s,1,eList)));
        data2 = cell2mat(squeeze(data(s,2,eList)));
    end
    
    [testingIndices1,testingIndices2] = getIndices(size(data1,2),size(data2,2),numFolds,useEqualStimRepsFlag);

    % Get the weightVector and projections of data1 and data2
    n1=length(electrodeArrayList{s}{1}); n2=length(electrodeArrayList{s}{2}); % Needed only for TransformType==1
    [weightVector,p1,p2,allIndices1,allIndices2] = getProjectionsAndWeightVector(data1,data2,testingIndices1,testingIndices2,transformType,n1,n2);
    newData{s,1} = p1(allIndices1);
    newData{s,2} = p2(allIndices2);
    
    for c=3:numConditions
        if dataTypeNum==1
            dataTMP = cell2mat(squeeze(data(s,c,eList)))./diff(timePeriod);
        else
            dataTMP = cell2mat(squeeze(data(s,c,eList)));
        end
        newData{s,c} = dataTMP' * weightVector;
    end
end
end
function newData = normalizeData(data)

[numSessions,numConditions] = size(data);
newData=cell(numSessions,numConditions);

for s=1:numSessions
    mX = mean(data{s,2});   % Normalized w.r.t H1V condition
    %   stdX = sqrt((var(data{s,1}) + var(data{s,2}))/2);
    stdX= std(data{s,2});
    
    
    for c=1:numConditions
        newData{s,c} = (data{s,c}-mX)/stdX;
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
function [weightVector,projections1,projections2,fullSetIndices1,fullSetIndices2] = getProjectionsAndWeightVector(data1,data2,testingIndices1,testingIndices2,transformType,n1,n2)

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
    
    if transformType==1 % Simple Averaging of sides
        weightVectorTMP(:,i) = [ones(n1,1)/n1; -ones(n2,1)/n2];
        
    else
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
            
            % LDA using fitcdiscr function of MATLAB
            Mdl = fitcdiscr(dataAll',labelAll);
            weightVectorTMP(:,i) = Mdl.Coeffs(1,2).Linear;
        end
    end
    
    projections1(t1) = data1(:,t1)' *  weightVectorTMP(:,i);
    projections2(t2) = data2(:,t2)' *  weightVectorTMP(:,i);
end

weightVector = mean(weightVectorTMP,2);
end