sbxPath = 'D:\2photon\PNS2\211007\PNS2_211007_030.sbxfix';
%load('D:\2photon\NaV18PNSGcamp8s\NaV18PNSGcamp8s_210827_000.mat'); %
sbxInfo = MakeInfoStruct('D:\2photon\PNS2\211007\PNS2_211007_030.mat');

WriteSbxPlaneTif(sbxPath, sbxInfo, 12, 'dir','D:\2photon\PNS2\211007\', 'name',sbxInfo.exptName, 'verbose',true, 'chan','both', 'binT',1, 'overwrite',true);  

%% Extract all sbx files within a folder
mainDir = 'D:\2photon\NaVAi6-2\211026\';
[sbxName, sbxPath] = FileFinder( mainDir, 'type','sbx' );
for s = 1:numel(sbxPath) % flip(1:numel(sbxPath))
    sbxInfoPath = sbxPath{s};
    sbxInfoPath(end-2:end) = 'mat';
    % Move all files to a subfolder
    tempInfo = MakeInfoStruct( sbxInfoPath );
    subDir = strcat(mainDir, tempInfo.fileName,'\');
    mkdir( subDir );
    [~,tempFiles] = FileFinder(mainDir, 'contains', tempInfo.fileName ); % , '*'
    for f = 1:numel(tempFiles)
        movefile(tempFiles{f}, subDir)
    end
    
    % Unpack the sbx file in the new subfolder
    [~,  newSBXpath] = FileFinder(subDir, 'type','sbx', 'keepExt',true);
    projPath = strcat(subDir, tempInfo.fileName, '_rawProj.tif');
    if ~tempInfo.optotune_used
        WriteSbxPlaneTif(newSBXpath{1}, tempInfo, 1, 'dir',subDir, 'name',sbxName{s}, 'verbose',true, 'chan','both', 'binT',8, 'overwrite',false);  % tempStack =
        WriteSbxProjection(newSBXpath{1}, tempInfo, projPath, 'verbose',true, 'chan','both', 'overwrite',true);
    else
        tempInfo = FixSBX(newSBXpath{1}, tempInfo, false);
        fixPath = strcat(newSBXpath{1},'fix');
        WriteSbxProjection(fixPath, tempInfo, projPath, 'verbose',true, 'chan','both', 'overwrite',true); % 'both'
        %WriteSbxPlaneTif(fixPath, tempInfo, 6, 'dir',subDir, 'name',sbxName{s}, 'verbose',true, 'chan','both', 'binT',1, 'overwrite',false); % fixPath
    end
end