% This function saves the detection performance(%correct) of the subject in all the sessions.


function saveBehavior(folderSourceString)
if ~exist('folderSourceString','var');   folderSourceString='E:\Mayo';       end

folderNameOut = fullfile(folderSourceString,'Data','segmentedData','behavior');
makeDirectory(folderNameOut);
fileNameStringList=getAttentionExperimentDetails;
fileNameStringListAll=cat(2,fileNameStringList{1},fileNameStringList{2});

for i=1:length(fileNameStringListAll)
    disp(['Working on ' fileNameStringListAll{i}]); 
    [perCorrect,uniqueOrientationChangeDeg,goodIndexList,orientationChangeDeg,reactionTimeMS]=getBehavior(fileNameStringListAll{i},folderSourceString,1);
    
    fileNameSave=fullfile(folderNameOut,[fileNameStringListAll{i} '_perCorrect']);
    save(fileNameSave,'perCorrect');
    
    fileNameSave2=fullfile(folderNameOut,[fileNameStringListAll{i} '_reactionTimes']);
    save(fileNameSave2,'uniqueOrientationChangeDeg','goodIndexList','orientationChangeDeg','reactionTimeMS');
end
end