function sbxInfo = FixSBX(sbxInputPath, sbxInfo, flipZ, overwrite)
if nargin < 3, flipZ = false; end
if nargin < 4, overwrite = false; end
fprintf('\n%s', sbxInputPath)
sbxFixPath = strcat(sbxInputPath, 'fix'); %sbxInputPath; %
[~, chan] = DeterminePMT('both', sbxInfo); % usePMTind
%chan = DetermineChanSBX(sbxInfo);
if ~exist(sbxFixPath, 'file') || overwrite
    % Load the planar data, shift the planes
    tic;
    planeData = cell(1,sbxInfo.Nplane);
    for z = 1:sbxInfo.Nplane
        planeData{z} = WriteSbxPlaneTif(sbxInputPath, sbxInfo, z, 'verbose',true, 'dir','', 'name',sbxInfo.fileName, 'chan',chan, 'RGB',false, 'overwrite',false ); % sbxInfo.dir
    end
    toc
    % reshape data to [chan x row x column x z x scan]
    if sbxInfo.nchan > 1
        planeData = cat( 5, planeData{:} );
        planeData = flip( permute(planeData, [4,1,2,5,3]), 1); %  % switch from RGB to PMT chan order % permute(planeData, [4,1,2,5,3]);
        planeData = circshift( planeData, -1, 4);
        if flipZ, planeData = flip(planeData, 4); end
        planeData = reshape(planeData, [size(planeData,[1,2,3]), prod(size(planeData,[4,5]))] );
        tic
        rw = pipe.io.RegWriter(sbxFixPath, sbxInfo, '.sbxfix', true);
        fprintf('\nWriting sbxfix file');
        rw.write( planeData ); % rw.write(squeeze(uint16(tempScan)));
        toc
        rw.delete;
    else
        planeData = cat( 4, planeData{:} );
        planeData = circshift( planeData, -1, 4);
        planeData = permute(planeData, [1,2,4,3]); % permute(planeData, [5,1,2,4,3]); 
        planeData = reshape(planeData, [size(planeData,[1,2]), prod(size(planeData,[3,4]))] );
        fprintf('\nWriting %s', sbxFixPath);
        rw = pipe.io.RegWriter(sbxFixPath, sbxInfo, '.sbxfix', true);
        tic
        rw.write( planeData );
        toc
        rw.delete;
    end
    toc  
else
    fprintf('%s already exists', sbxFixPath);
end

% Update sbxInfo to point to sbxfix
%sbxInfo.path = sbxFixPath;
%projPath = sprintf('%s%s_fixProj.tif', sbxInfo.dir, sbxInfo.fileName);
%WriteSbxProjection(sbxInfo.path, sbxInfo, projPath, 'dir','', 'name',sbxInfo.fileName, 'type','fix', 'chan','green', 'firstScan',1, 'Nscan',sbxInfo.Nscan, 'overwrite',true);

end