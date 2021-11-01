function [stackOut, stackChan] = WriteSbxPlaneTif(sbxPath, sbxInfo, z, varargin) 
IP = inputParser;
addRequired( IP, 'sbxPath', @ischar )
addRequired( IP, 'sbxInfo', @isstruct)
addRequired( IP, 'z', @isnumeric)
addParameter( IP, 'chan', 'both', @ischar ) % 'green', 'red', 'both'. for scanbox, 1 = green, 2 = red. -1 = both
addParameter( IP, 'firstScan', 1, @isnumeric )
addParameter( IP, 'Nscan', -1, @isnumeric )
addParameter( IP, 'edges', [0,0,0,0], @isnumeric ) % [left, right, top, bottom]
addParameter( IP, 'scale', 1, @isnumeric ) %
addParameter( IP, 'binT', 1, @isnumeric ) 
addParameter( IP, 'zeros', true, @islogical )
addParameter( IP, 'dir', '', @ischar ) 
addParameter( IP, 'name', '', @ischar )
addParameter( IP, 'type', '', @ischar ) 
addParameter( IP, 'monochrome', true, @islogical ) % should we write individual monochrome tifs from each channel?
addParameter( IP, 'RGB', true, @islogical ) % should we write an RGB tif?
addParameter( IP, 'rescale', false, @islogical )
addParameter( IP, 'verbose', false, @islogical )
addParameter( IP, 'overwrite', false, @islogical )
parse( IP, sbxPath, sbxInfo, z, varargin{:} ); 
tifDir = IP.Results.dir;
writeChan = IP.Results.chan;  %1; %1 = red, 2 = green
firstScan = IP.Results.firstScan;
Nscan = IP.Results.Nscan;
edges = IP.Results.edges;
scaleFactor = IP.Results.scale;
binT = IP.Results.binT;
monochrome = IP.Results.monochrome;
RGBtoggle = IP.Results.RGB;
rescaleIntToggle = IP.Results.rescale;
allowZeros = IP.Results.zeros;
verbose = IP.Results.verbose;
exptName = IP.Results.name;
sbxType = IP.Results.type;
overwrite = IP.Results.overwrite;
outputToggle = nargout > 0;
if verbose, tic; end
if ~isempty(sbxType), sbxType = strcat('_',sbxType); end
if isempty(tifDir) || isempty(exptName)
    if verbose, fprintf('\ntifDir and/or exptName is empty. No tifs will be saved'); end
    monochrome = false;
    RGBtoggle = false;
end

% DETERMINE WHICH CHANNELS TO WRITE
pmtName = {'green','red'};
[usePMT, ~] = DeterminePMT(writeChan, sbxInfo);
Npmt = numel(usePMT);

% DETERMINE # OF SCANS, BINNING, ETC
if Nscan < 0,  Nscan = sbxInfo.totScan;  end
if binT ~= 1  
    [chunkLims, Nchunk, chunkLength] = MakeChunkLims(firstScan, firstScan+Nscan-1, sbxInfo.totScan, 'size',binT );  
end

% GET DATA FROM EACH CHANNEL AND WRITE MONOCHROME TIFFS
stackChan = cell(1,2);
for pmt = usePMT
    tifPath = sprintf('%s%s%s_Z%02d_%s.tif', tifDir, exptName, sbxType, z, pmtName{pmt});
    fileExists = exist(tifPath,'file');
    if outputToggle && ~overwrite && fileExists
        % Load stackChan from tif
        if verbose, fprintf('\nLoading %s', tifPath);  end
        stackChan{pmt} = loadtiff( tifPath );
        if ~allowZeros,  stackChan{pmt}(stackChan{pmt} == 0) = NaN;  end % Set zeros to NaN to avoid averaging over blank frames elsewhere
    elseif ~(~outputToggle && ~overwrite && fileExists)
        if binT == 1
            % Load, crop and resize the sbx data
            if verbose, fprintf('\nLoading %s (plane %i, first scan = %i, Nscan = %i, chan = %i)    ', sbxPath, z, firstScan, Nscan, pmt ); end
            if sbxInfo.Nplane > 1
                stackChan{pmt} = readSBX(sbxPath, sbxInfo, firstScan, Nscan, pmt, z); 
            else
                stackChan{pmt} = readSBX(sbxPath, sbxInfo, firstScan, Nscan, pmt, []); 
            end
        else
            if verbose, fprintf('\nLoading %s (plane %i, first scan = %i, Nscan = %i, chan = %i, averaging every %i frames)    ', sbxPath, z, firstScan, Nscan, pmt, binT ); end 
            w = waitbar(0, sprintf('Binning %s', sbxPath ));
            for c = flip(1:Nchunk)
                if sbxInfo.Nplane > 1
                    tempChunk = readSBX(sbxPath, sbxInfo, chunkLims(c,1), chunkLength(c), pmt, z);  
                else
                    tempChunk = readSBX(sbxPath, sbxInfo, chunkLims(c,1), chunkLength(c), pmt, []); 
                end
                stackChan{pmt}(:,:,c) = mean(tempChunk, 3);
                waitbar((Nchunk-c)/Nchunk, w);
                %
            end
            delete(w);
            stackChan{pmt} = uint16(stackChan{pmt}); 
        end
        if any(edges > 0), stackChan{pmt} = stackChan{pmt}(edges(3)+1:end-edges(4),edges(1)+1:end-edges(2),:,:,:); end  % crop the data
        if scaleFactor ~= 1, stackChan{pmt} = imresize3( stackChan{pmt}, [size(stackChan{pmt},1)/scaleFactor, size(stackChan{pmt},2)/scaleFactor, size(stackChan{pmt},3)] );  end  % resize       
        % Save channel data to monochrome tif (optional)
        if monochrome
            if ~fileExists || overwrite
                if verbose, fprintf('\nWriting %s    ', tifPath ); end
                WriteTiff(stackChan{pmt}, tifPath); %pipe.io.writeTiff(stackChan{pmt}, tifPath);  %saveastiff( uint16(stackChan), tifPath, RGBOpt ); %
            end
        end
    end
end 
% MERGE CHANNELS, SWITCH TO RGB, AND (OPTIONAL) WRITE SINGLE RGB TIFF
if Npmt > 1 %sbxInfo.nchan == 2  %NgoodChan > 1 && 
    % reverse color order to RGB
    pmtName = flip(pmtName);
    stackChan = stackChan([2,1]);
    % Write an RGB movie (optional)
    if RGBtoggle 
        rgbPath = sprintf('%s%s%s_Z%02d_RGB.tif', tifDir, exptName, sbxType, z);
        if (~exist(rgbPath,'file') || overwrite)
            if verbose, fprintf('\nWriting %s', rgbPath); end
            rgbTiff = cell(1,2);
            for pmt = 1:2
                if ~rescaleIntToggle
                    rgbTiff{pmt} = uint8(stackChan{pmt}/256);
                else
                    % Adjust data for uint8 RGB tiff if needed
                    chanLower = prctile(stackChan{pmt}(:), 1);
                    chanUpper = max(stackChan{pmt}(:)); %prctile(stackChan{chan}(:), 1);
                    fprintf('\nRescaling channel %s: [%i, %i] -> [0, 255]', pmtName{pmt}, chanLower, chanUpper);
                    rgbTiff{pmt} = uint8(rescale(stackChan{pmt}, 0, 2^8-1, 'inputMin',chanLower)); % min(stackChan{chan}(:))
                end
            end
            rgbTiff = cat(4, rgbTiff{:});
            rgbTiff(end,end,end,3) = 0; % need to add a blank blue channel for tif format
            if verbose, fprintf('\nWriting %s    ', rgbPath ); end
            WriteTiff(rgbTiff, rgbPath); % pipe.io.writeTiff(rgbTiff, rgbPath);
        end
    end
    stackOut = cat(4, stackChan{:}); 
else
    stackOut = stackChan{usePMT}; 
end
if ~allowZeros,  stackOut(stackOut == 0) = NaN;  end % Set zeros to NaN to avoid averaging over blank frames elsewhere
if verbose, toc, end  
%}
end