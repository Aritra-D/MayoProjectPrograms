% June 9, 2020: Adding functionality to also run the attention axis using
% codes provided by Patrick J Mayo.

% This Program is used in displayFigure5, displayFigure5v2 to save the
% output of this function.
% The difference between this function and
% displayHvsMPopulation is that here pooled variance is used for uncorrelated LDA case

% dataTypeNum - 1: spikes, 2: gamma, 3: alpha
% transformType - 1: simple averaging across sides, 2: Mean-difference, 3:
% LDA-uncorrelated, 4: LDA-covariance, 5: Attention Axis 

% Adding an option to do cross-validation. Set numFolds to 1 if you do not want to do cross-validation
% Also adding an option to equalize the stimulus repeats for H0V and H1V

% Adding an option to work with a smaller set of electrodes

function [mData,typeList1,dataTMP] = displayHvsMAnalysisPopulation2(dataTypeNum,transformType,regFlag,trialCutoff,normalizeFlag,numFolds,useEqualStimRepsFlag,numElectrodesToUse,colorToUse)

if ~exist('transformType','var');       transformType=4;                end
if ~exist('regFlag','var');             regFlag=0;                      end
if ~exist('trialCutoff','var');         trialCutoff=15;                 end
if ~exist('normalizeFlag','var');       normalizeFlag=1;                end
if ~exist('numFolds','var');            numFolds=5;                     end
if ~exist('useEqualStimRepsFlag','var'); useEqualStimRepsFlag=0;        end
if ~exist('numElectrodesToUse','var');  numElectrodesToUse=[];          end
if ~exist('colorToUse','var');          colorToUse=[];                  end

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

if isempty(colorToUse)
    colorNameList = 'rbg';
    colorName = colorNameList(dataTypeNum);
else
    colorName = colorToUse;
end

if dataTypeNum==1
    dataTMP = transformData(spikeData,transformType,dataTypeNum,timePeriod,electrodeArrayList,numFolds,useEqualStimRepsFlag,numElectrodesToUse,regFlag);
elseif dataTypeNum==2
    dataTMP = transformData(gammaData,transformType,dataTypeNum,timePeriod,electrodeArrayList,numFolds,useEqualStimRepsFlag,numElectrodesToUse,regFlag);
elseif dataTypeNum==3
    dataTMP = transformData(alphaData,transformType,dataTypeNum,timePeriod,electrodeArrayList,numFolds,useEqualStimRepsFlag,numElectrodesToUse,regFlag);
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
typeList1ordered = cell(1,numConditions);
conditionOrder = [4 1 11 10 3 2 12 9 8 5 7 6]; 

for c=1:numConditions
    x = mData{c};
    mtmp = mean(x);
    ntmp = length(x);
    stmp = std(x)/sqrt(ntmp);
    errorbar(conditionOrder(c),mtmp,stmp,[colorName 'o']); hold on;
    text(conditionOrder(c)-0.2,mtmp+stmp+0.005,num2str(ntmp));
    mDataTMP(conditionOrder(c)) = mtmp;
    typeList1ordered(conditionOrder(c)) = typeList1(c);
end
plot(mDataTMP,colorName);
h=gca;
set(h,'XTick',1:12,'XTickLabel',typeList1ordered);

end

function newData = transformData(data,transformType,dataTypeNum,timePeriod,electrodeArrayList,numFolds,useEqualStimRepsFlag,numElectrodesToUse,regFlag)

[numSessions,numConditions,~] = size(data);
newData=cell(numSessions,numConditions);

for s=1:numSessions
    eList = [electrodeArrayList{s}{1};electrodeArrayList{s}{2}]';
    n1=length(electrodeArrayList{s}{1}); n2=length(electrodeArrayList{s}{2}); % Needed only for TransformType==1
    eSideList = [zeros(1,n1) ones(1,n2)];
    if ~isempty(numElectrodesToUse)
        if n1+n2>numElectrodesToUse % Use a subset of electrodes
            eIndices = randperm(n1+n2);
            eIndices = sort(eIndices(1:numElectrodesToUse));
            eList = eList(eIndices);
            n1 = length(find(eSideList(eIndices)==0));
            n2 = length(find(eSideList(eIndices)==1));
        end
    end
    
    if dataTypeNum==1
        data1 = cell2mat(squeeze(data(s,1,eList)))./diff(timePeriod); % H0V
        data2 = cell2mat(squeeze(data(s,2,eList)))./diff(timePeriod); % H1V
    else
        data1 = cell2mat(squeeze(data(s,1,eList)));
        data2 = cell2mat(squeeze(data(s,2,eList)));
    end
    
    [testingIndices1,testingIndices2] = getIndices(size(data1,2),size(data2,2),numFolds,useEqualStimRepsFlag);

    % Get the weightVector and projections of data1 and data2
    [weightVector,p1,p2,allIndices1,allIndices2,multiDimMean1,multiDimMean2] = getProjectionsAndWeightVector(data1,data2,testingIndices1,testingIndices2,transformType,n1,n2,regFlag);
    newData{s,1} = p1(allIndices1);
    newData{s,2} = p2(allIndices2);
    
    for c=3:numConditions
        if dataTypeNum==1
            dataTMP = cell2mat(squeeze(data(s,c,eList)))./diff(timePeriod);
        else
            dataTMP = cell2mat(squeeze(data(s,c,eList)));
        end
        if transformType==5 % Attention Axis
            newData{s,c} = projpointline(dataTMP',multiDimMean1',multiDimMean2');
        else
            newData{s,c} = dataTMP' * weightVector;
        end
    end
end
end
function newData = normalizeData(data)

[numSessions,numConditions] = size(data);
newData=cell(numSessions,numConditions);

for s=1:numSessions
    mX = mean(data{s,2});   % Normalized w.r.t H1V condition
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
function [weightVector,projections1,projections2,fullSetIndices1,fullSetIndices2,multiDimMean1,multiDimMean2] = getProjectionsAndWeightVector(data1,data2,testingIndices1,testingIndices2,transformType,n1,n2,regFlag)

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
multiDimMean1TMP = zeros(D,numFolds); % needed for the Attention Axis
multiDimMean2TMP = zeros(D,numFolds); % needed for the Attention Axis

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
        
        if transformType==2 % Only mean difference
            weightVectorTMP(:,i) = mean(d1,2) - mean(d2,2);
            
        elseif transformType==3 % LDA for uncorrelated case       
            meanDiff = mean(d1,2) - mean(d2,2);
            var1 = var(d1,[],2);
            var2 = var(d2,[],2);
            pooledVar = ((m1-1)*var1 + (m2-1)*var2)/(m1+m2-2); %Pooled Variance
            
            weightVectorTMP(:,i) = meanDiff./pooledVar;
            
        elseif transformType==4 % LDA

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
    end
    
    if transformType==5 % JPM Attention Axis - Done separately since it provides both projections and weights
        multiDimMean1TMP(:,i) = nanmean(d1, 2); % take mean across Hit trials only at Location 0
        multiDimMean2TMP(:,i) = nanmean(d2, 2);
        
        % % First input must be TRIALS (rows) x CHANNELS (columns)
        [projections1(t1),weightVectorTMP(:,i)] = projpointline(data1(:,t1)',multiDimMean1TMP(:,i)',multiDimMean2TMP(:,i)'); % FUNCTION DEFINED BELOW
        projections2(t2) = projpointline(data2(:,t2)',multiDimMean1TMP(:,i)',multiDimMean2TMP(:,i)'); % FUNCTION DEFINED BELOW

    else
        projections1(t1) = data1(:,t1)' *  weightVectorTMP(:,i);
        projections2(t2) = data2(:,t2)' *  weightVectorTMP(:,i);
    end
end

weightVector = mean(weightVectorTMP,2);
multiDimMean1 = mean(multiDimMean1TMP,2);
multiDimMean2 = mean(multiDimMean2TMP,2);
end

% Code provided by Patrick J Mayo.
function [t, wvTEST] = projpointline(p,a,b) % t is projection

% projection of point p onto line between points a and b
% output values range from 0-1
% If input is a matrix, p must be TRIALS (rows) x CHANNELS (columns)

if size(p,1) == 0
    t=[];
elseif size(p,1) > 1 % if data input matrix is more than a single trial

    denom = norm(a-b)^2;  % denominator from above; single value regardless of points
    
    if size(p,1) ~= size(a,1) && size(p,1) ~= size(b,1) % if point matrix is not same size at mean matrices
        ar = repmat(a, size(p,1), 1); % copy point A values
        br = repmat(b, size(p,1), 1); % copy point B values
        numer = dot((ar-br),ar-p, 2); % dot product along the second dimension
        t=numer/denom;
        
        wvTEST = (a-b)/norm(a-b); % normalized vector from b to a. 
        
    else % if input matrices (a, b, p) are already all the same size, no need to use repmat
        numer = dot((a-b),a-p, 2); % dot product along the second dimension
        t=numer/denom;
        
        wvTEST = (a-b)/norm(a-b); % normalized vector from b to a. 
    end
    
else % if only a single trial, use Marlene's original line of code
    t=dot((a-b),(a-p))/norm(a-b)^2;
    wvTEST = (a-b)/norm(a-b); % normalized vector from b to a. 
end

t=-t; % Projections are reversed in sign to be consistent with other measures
end