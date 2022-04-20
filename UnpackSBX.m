sbxPath = 'D:\2photon\NaVAi6-2\211210\NaVAi6-2_211210_007\NaVAi6-2_211210_007.sbxfix';
%load('D:\2photon\NaV18PNSGcamp8s\NaV18PNSGcamp8s_210827_000.mat'); %
sbxInfo = MakeInfoStruct('D:\2photon\NaVAi6-2\211210\NaVAi6-2_211210_007\NaVAi6-2_211210_007.mat');

WriteSbxPlaneTif(sbxPath, sbxInfo, 7, 'dir','D:\2photon\NaVAi6-2\211210\NaVAi6-2_211210_007\', 'name',sbxInfo.fileName, 'verbose',true, 'chan','both', 'binT',5, 'overwrite',true);  

%% Extract all sbx files within a folder
mainDir = 'D:\2photon\NaVGCaMP-01\';
[sbxName, sbxPath] = FileFinder( mainDir, 'type','sbx' );
flipZ = true;
for s = flip(1:numel(sbxPath)) %1:numel(sbxPath) % flip(1:numel(sbxPath))
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
        tempInfo = FixSBX(newSBXpath{1}, tempInfo, flipZ, false);
        fixPath = strcat(newSBXpath{1},'fix');       
        WriteSbxProjection(fixPath, tempInfo, projPath, 'verbose',true, 'chan','both', 'overwrite',true, 'RGB',false, 'monochrome',false); % 'both'
        %WriteSbxPlaneTif(fixPath, tempInfo, 6, 'dir',subDir, 'name',sbxName{s}, 'verbose',true, 'chan','both', 'binT',1, 'overwrite',false); % fixPath
    end
end