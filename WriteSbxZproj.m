function projData = WriteSbxZproj(sbxPath, sbxInfo, projPath, varargin) % firstScan, Nscan, pmt, 

IP = inputParser;
addRequired( IP, 'sbxPath', @ischar )
addRequired( IP, 'sbxInfo', @isstruct )
addRequired( IP, 'projPath', @ischar )
addParameter( IP, 'type', 'mean', @ischar )
addParameter( IP, 'edge', [100, 100, 120, 120], @isnumeric ) % [60,60,40,40]
addParameter( IP, 'firstScan', 1, @isnumeric )
addParameter( IP, 'lastScan', sbxInfo.totScan, @isnumeric )
addParameter( IP, 'chan', 'green', @ischar )
addParameter( IP, 'z', 1:sbxInfo.Nplane, @isnumeric )
addParameter( IP, 'binT', 1, @isnumeric ) 
addParameter( IP, 'write', true, @islogical )
addParameter( IP, 'overwrite', false, @islogical ) % for scanbox, 1 = green, 2 = red. -1 = both
parse( IP, sbxPath, sbxInfo, projPath, varargin{:} ); % mouse, exptDate,
projType = IP.Results.type;
edge = IP.Results.edge;
zProj = IP.Results.z;
firstScan = IP.Results.firstScan;
lastScan = IP.Results.lastScan;
Nscan = numel(firstScan:lastScan);
projChan = IP.Results.chan;
projPMT = DeterminePMT(projChan, sbxInfo);
binT = IP.Results.binT;
write_tiff = IP.Results.write;
overwrite = IP.Results.overwrite;
if exist(projPath,'file') && ~overwrite, write_tiff = false; end

if nargout > 0 || write_tiff % don't bother if there will be no output
    if projPMT == -1 && sbxInfo.nchan == 2
        Nchan = 2;
        projDim = 4;
    else
        Nchan = 1;
        projDim = 3;
    end
    %load each volume scan (parallelized) and z-project
    tic
    A = cell(Nscan,1); % firstScan+Nscan
    h = parfor_progressbar(Nscan,'z-projecting...');
    parfor scan = firstScan:firstScan+Nscan-1
        vol = readSBX(sbxPath, sbxInfo, scan, 1, projPMT, []);
        if Nchan == 1
            vol = vol(:,:,zProj);
        else
            vol = vol(:,:,:,zProj);
        end
        vol(vol==0) = NaN;
        if strcmpi(projType,'mean')
            slice = squeeze( mean(vol,projDim,'omitnan') );
        elseif strcmpi(projType,'max')
            slice = squeeze( max(vol,[],projDim,'omitnan') );
        elseif strcmpi(projType,'median')
            slice = squeeze( median(vol,projDim,'omitnan') );
        else
            disp('invalid projection type');
        end
        slice(isnan(slice)) = 0;
        A{scan} = slice;
        h.iterate(1);
    end
    delete(h);
    toc

    %concatenate each projected volume
    projData = cat(projDim, A{:});
    clearvars A;
    if ndims(projData) == 4,  projData = permute(projData, [2,3,4,1]); end % permute(projData, [2,3,1,4])
    projData = projData(edge(3)+1:end-edge(4),  edge(1)+1:end-edge(2), :, :);% crop edges
    projLength = size(projData, 3);

    % Temporal binning
    if binT ~= 1
        if rem( size(projData, 3), binT ) ~= 0
            cutScans = projLength-binT*floor(size(projData, 3)/binT);
            fprintf('\nNumber of scans (%i) is not divisible by %i: cutting off first %i scans',projLength, binT, cutScans);
            projData = projData(:,:,cutScans+1:end,:);
        end
        if Nchan == 1
            projData = reshape(projData, size(projData, 1), size(projData, 2), binT, []);
            projData = squeeze(mean(projData, 3));
        else
            projData = reshape(projData, size(projData, 1), size(projData, 2), binT, [], 2); %reshape(projData, size(projData, 1), size(projData, 2), [], binT, 2);
            projData = squeeze(mean(projData, 3));
        end
    end

    %write the projection to a tif file
    if write_tiff 
        fprintf('\nWriting %s\n', projPath);
        if Nchan > 1
            dataMin = min(projData(:)); dataMax = max(projData(:)); 
            projData(end,end,end,3) = 0;  % need a blue channel for RGB data  
            WriteTiff(uint8( rescale(projData, 0, 255, 'InputMin',dataMin, 'InputMax',dataMax)), projPath);
        else
            WriteTiff(uint16(projData), projPath); 
        end 
    end
else
    fprintf('\n%s: No output, nor tiff requested', projPath);
end
end