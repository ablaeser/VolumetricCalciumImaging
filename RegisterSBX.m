function [reg_transforms, varargout] = RegisterSBX(mov_path, sbxInfo, refImage, params, varargin)
%SBXALIGNTURBOREGCORE aligns a file (given by path) using ImageJ's TurboReg
% adapted from pipe.reg.turboreg (Andy Blaeser, Jan 2021)
%   NOTE: hardcoded path to ImageJ.
IP = inputParser;
addRequired( IP, 'mov_path', @ischar )
addRequired( IP, 'sbxInfo', @isstruct )
addRequired( IP, 'refImage', @isnumeric )
addRequired( IP, 'params', @isstruct )
addOptional(IP, 'z', 1, @isnumeric)  % Which optotune level to align
addParameter(IP, 'method', 'affine', @ischar )
addParameter(IP, 'firstScan', 1)  % The frame to start reading from
addParameter(IP, 'Nscan', 1)  % Number of frames to read 
addParameter(IP, 'overwrite', false, @islogical)
addParameter(IP, 'intTif', false, @islogical)
addParameter(IP, 'finalTif', false, @islogical)
addParameter(IP, 'verbose', true, @islogical)
addParameter(IP, 'dir', '', @ischar)
addParameter(IP, 'name', '', @ischar)
parse( IP, mov_path, sbxInfo, refImage, params, varargin{:} ); 
z = IP.Results.z;
firstScan = IP.Results.firstScan;
Nscan = IP.Results.Nscan;
lastScan = firstScan + Nscan - 1;
saveName = IP.Results.name;
writeIntermediateTif = IP.Results.intTif;
writeFinalTif = IP.Results.finalTif;
verbose = IP.Results.verbose;
overwrite = IP.Results.overwrite;
RegTempDir = sprintf('%sRegTemp\\', sbxInfo.dir); [~,~] = mkdir(RegTempDir);
chunkTempDir = sprintf('%sChunks\\', RegTempDir); [~,~] = mkdir(chunkTempDir);
chanName = {'green','red'};
params.refChanInd = find( strcmpi(params.refChan, chanName));
if isempty(params.name)
    nameStr = '';
else
    nameStr = ['_',params.name];
end
tic
% Crop the reference image
rawRef = uint16(refImage(params.edges(3)+1:end-params.edges(4), params.edges(1)+1:end-params.edges(2)));
preRef = rawRef;

% Break the data into chunk of scans to avoid memory issues and improve TurboReg performance
if params.chunkSize == 0, params.chunkSize = Nscan+1; end %  sbxInfo.totScan
[chunkScan, Nchunk, chunkLength] = MakeChunkLims(firstScan, firstScan+Nscan-1, sbxInfo.totScan, 'size',params.chunkSize );

if nargout > 1 || writeFinalTif
    unregData = zeros(size(rawRef,1), size(rawRef,2), numel(firstScan:lastScan) );
    regData = zeros(size(rawRef,1), size(rawRef,2), numel(firstScan:lastScan) );
end
clearvars chunkLims
%{
chunkLims(:,1) = (firstScan:params.chunkSize:lastScan)'; 
chunkLims(:,2) = chunkLims(:,1) + params.chunkSize - 1;
chunkLims( chunkLims(:,1) > sbxInfo.totScan, : ) = [];
if chunkLims(end,2) > sbxInfo.totScan, chunkLims(end,2) = sbxInfo.totScan; end
if chunkLims(end,2) > lastScan, chunkLims(end,2) = lastScan; end
Nchunk = size(chunkLims, 1);
%}
K = 0; % counter to keep track of regData frames

% Set up registration method
reg_transforms = cell(1, sbxInfo.totScan);
if strcmpi(params.method, 'affine')
    % Run imagej and prepare turboreg
    if verbose, fprintf('\nRunning ImageJ  '); end
    pipe.lab.runimagej(); % Hardcoded path to ImageJ
    trhl = TurboRegHL_();
else %if strcmpi(params.method, 'translation') || strcmpi(params.method, 'rigid')
    [optimizer, metric] = imregconfig('monomodal');
    targets = '';
    targetstr = '';
    szstr = '';
    alignstr = '';
end
% Set up temporal averaging
weightVector = ones(1, params.avgT);
if params.avgT > 0
    if( rem(params.avgT,2) == 0 ),  error('params.avgT should be set to an odd number!');  end 
    avgTpad = round(params.avgT/2) - 1; % how many scans before/after the scan do we average over?
    if params.avgTsigma > 0
        tGauss = -avgTpad:avgTpad; 
        weightVector = exp(-(tGauss).^2/(2*params.avgTsigma^2));  
        weightVector = double(weightVector/sum(weightVector)); %plot(tGauss, gaussWeight)       
    end
else
    avgTpad = 0;
end

if Nscan > 1, w = waitbar(0,'Performing registration'); end
reg_transforms = repmat( {affine2d}, 1, Nscan );
for c = 1:Nchunk
    chunkSavePath = sprintf('%sz%02i%s_chunk%04i.mat', chunkTempDir, z, nameStr, c);
    chunkScans = chunkScan(c,1):chunkScan(c,2);
    NchunkScans = chunkLength(c); %size(chunkScans, 3); 
    padScans = chunkScan(c,1)-avgTpad:chunkScan(c,2)+avgTpad; %
    padScans(padScans < 1) = []; padScans(padScans > sbxInfo.totScan) = [];
    NpadScans = numel(padScans);
    if ~exist(chunkSavePath, 'file') || overwrite
        % LOAD DATA
        if verbose, fprintf('\nchunk %i / %i:  Loading %s (z = %i, scans %i - %i)', c, Nchunk, mov_path, z, chunkScan(c,1), chunkScan(c,2)  ); end
        if sbxInfo.Nplane > 1
            rawData = readSBX(mov_path, sbxInfo, padScans(1), numel(padScans), params.refChanInd, z); % , chunkLims(c,1), diff(chunkLims(c,:))+1
        else
            rawData = readSBX(mov_path, sbxInfo, padScans(1), numel(padScans), params.refChanInd, []);
        end
         
        % PRE-PROCESS DATA AND REFERENCE IMAGE (note, reference image gets unnecessarily processed for each chunk, consider preprocessing once outside chunk loop)
        % Crop edges
        if verbose, fprintf('\nCropping edges: [L, R, T, B] = [%i, %i, %i, %i]', params.edges(1), params.edges(2), params.edges(3), params.edges(4) ); end
        rawData = rawData(params.edges(3)+1:end-params.edges(4), params.edges(1)+1:end-params.edges(2), :);
        preData = rawData;
        %figure; subplot(1,2,1); imshow(rawRef,[]); subplot(1,2,2); imshow(preRef,[]); 
        
        % Low-pass filter 
        if params.lowpass > 0
            if verbose, fprintf('\nLowpass = %2.1f pix   ', params.lowpass ); end
            if c == 1
                preRef = imgaussfilt(preRef, params.lowpass);
            end
            preData = imgaussfilt(preData, params.lowpass);
        end
        
        % Median filter
        if any(params.medFilter)
            if verbose,  fprintf('\nMedian filtering: [%i, %i, %i]', params.medFilter(1), params.medFilter(2), params.medFilter(3)); end
            if c == 1
                preRef = medfilt2( preRef, params.medFilter(1:2) );
            end
            preData = medfilt3( preData, params.medFilter );
        end
        %figure; subplot(1,2,1); imshow(rawRef,[]); subplot(1,2,2); imshow(preRef,[]); impixelinfo;

        % Spatial downsampling
        if params.binXY ~= 1
            if verbose, fprintf('\nSpatial downsampling factor %i', params.binXY ); end
            preData = pipe.proc.binxy(preData, params.binXY);
            if c == 1
                preRef = pipe.proc.binxy(preRef, params.binXY);
            end
        end
        
        % Temporal averaging (not downsampling!)
        if params.avgT > 1 
            if verbose, fprintf('\nTemporal averaging over %i frames (sigma = %2.1f frames) ', params.avgT, params.avgTsigma); end
            [~,padChunkScans] = intersect(padScans, chunkScans); 
            avgData = cell(1,NchunkScans); % zeros(size(preData,1), size(preData,2), NchunkScans); %preData;
            parfor s = 1:NchunkScans % padChunkScans' %find(chunkScans == padScans)
                tempPadScans = padChunkScans(s)-avgTpad:padChunkScans(s)+avgTpad;
                badTempScans = find(tempPadScans< 1 | tempPadScans> NpadScans);
                tempPadScans(badTempScans) = []; 
                tempWeight = weightVector;
                tempWeight(badTempScans) = [];
                %plot(tempPadScans,  tempWeight); xlim([tempPadScans(1)-1, tempPadScans(end)+1]);
                tempWeight = repmat( permute(tempWeight, [1,3,2]), size(preData,1), size(preData,2) );
                avgData{s} = wmean(preData(:,:,tempPadScans), tempWeight, 3);  %mean(preData(:,:,tempPadScans), 3); % 
                %subplot(1,2,1); imshow( preData(:,:,padChunkScans(s)), []); 
                %subplot(1,2,2); imshow(avgData{s}, []); title(sprintf('s = %i', s));
                %pause(0.2);
            end
            avgData = cat(3, avgData{:});
            preData = avgData;
        end
        %figure; subplot(1,2,1); imshow(rawRef,[]); subplot(1,2,2); imshow(preRef,[]); impixelinfo;
        %toc
        
        % Equalize histograms to reference
        if params.histmatch
            if verbose, fprintf('\nHistogram matching enabled'); end
            figure; 
            subplot(2,3,1); imshow(preRef,[]); title('Reference Image');
            subplot(2,3,4); histogram(preRef(:)); set(gca,'Xscale','log', 'yscale','log');
            for i = 1:NchunkScans
                tempIm = preData(:,:,i);
                subplot(2,3,2); imshow(tempIm,[]);  title(sprintf('Non-matched data (i = %i)', i));
                subplot(2,3,5); histogram(tempIm(:)); set(gca,'Xscale','log', 'yscale','log');
                matchIm = imhistmatch(preData(:,:,i), preRef);
                preData(:,:,i) = matchIm;
                subplot(2,3,3); imshow(matchIm,[]);  title('Histogram-matched Image');
                subplot(2,3,6); histogram(matchIm(:));
                impixelinfo;
                pause; cla;
            end
        else
            if verbose, fprintf('\nHistogram matching disabled'); end
        end

        % High-pass filter
        if params.highpass > 0
            if verbose, fprintf('\nHighpass = %2.1f pix   ', params.highpass ); end
            if c == 1
                preRef = preRef - imgaussfilt(preRef, params.highpass);
            end
            preData = preData - imgaussfilt(preData, params.highpass);
        end
        
        % Write tifs of pre-affine-registered data and reference (optional)
        if writeIntermediateTif
            % Determine tif filenames, initialize regData (optional)
            if c == 1 % only need to save reference tif once since reference is the same for all chunks
                rawRefTifPath = sprintf('%s%s_rawRef.tif', RegTempDir, saveName );
                if verbose, fprintf('\nWriting %s', rawRefTifPath); end
                WriteTiff(uint16(rawRef), rawRefTifPath); %pipe.io.writeTiff(uint16(rawRef), rawRefTifPath); %saveastiff( uint16(rawRef), rawRefTifPath, GrayOpt);
                preRefTifPath = sprintf('%s%s_preRef.tif', RegTempDir, saveName );
                if verbose, fprintf('\nWriting %s', preRefTifPath); end
                WriteTiff(uint16(preRef), preRefTifPath); %saveastiff( uint16(preRef), preRefTifPath, GrayOpt);
            end
            if Nchunk > 1
                %rawDataTifPath = sprintf('%s%s_chunk%04i_rawData.tif', RegTempDir, saveName, c  );
                preDataTifPath = sprintf('%s%s_chunk%04i_preData.tif', RegTempDir, saveName, c );
                postDataTifPath = sprintf('%s%s_chunk%04i_postData.tif', RegTempDir, saveName, c );
            else
                %rawDataTifPath = sprintf('%s%s_rawData.tif', RegTempDir, saveName );
                preDataTifPath = sprintf('%s%s_preData.tif', RegTempDir, saveName );
                postDataTifPath = sprintf('%s%s_postData.tif', RegTempDir, saveName );
            end
            % Save raw tifs
            %if verbose, fprintf('\nWriting %s', rawDataTifPath); end
            %WriteTiff(uint16(rawData), rawDataTifPath); %saveastiff( uint16(preData), preDataTifPath, GrayOpt);
            if verbose, fprintf('\nWriting %s', preDataTifPath); end
            WriteTiff(uint16(preData), preDataTifPath);
        end

        % PERFORM REGISTRATION
        % pre-align using DFT (optional)
        if params.prereg > 0
            if verbose, fprintf('\nDFT pre-registration enabled (upsample = %i) ', params.prereg); end
            target_fft = fft2(double(preRef));
            dft_transforms = cell(1,NchunkScans); dft_reg = cell(1,NchunkScans); % data_fft = cell(1,NchunkScans); 
            parfor i = 1:NchunkScans
                [dft_transforms{i}, reg] = pipe.reg.dftcore(target_fft, fft2( double(preData(:,:,i)) ), params.prereg);
                dft_reg{i} = abs(ifft2(reg));
            end
            preData = cat(3, dft_reg{:});
            dft_transforms = cat(1, dft_transforms{:})';
        else
            if verbose, fprintf('\nPre-registration disabled'); end
        end

        % Final registration
        if strcmpi(params.method, 'affine')
            % Use TurboReg to perform affine registration
            if c == 1
                [Nrow, Ncol] = size(preRef);  % Get the dimensions of the preprocessed data chunk
                % Estimate targets the way turboreg does
                targets = [0.5*Ncol 0.15*Nrow 0.5*Ncol 0.15*Nrow 0.15*Ncol 0.85*Nrow 0.15*Ncol 0.85*Nrow  0.85*Ncol 0.85*Nrow 0.85*Ncol 0.85*Nrow];
                targets = round(targets);
                targetstr = sprintf('%i ', targets);
                % Create the text for the ImageJ macro
                szstr = sprintf('0 0 %i %i ', Ncol - 1, Nrow - 1);
                alignstr = sprintf('-align -window data %s -window ref %s -affine %s -hideOutput', szstr, szstr, targetstr); % -rigidBody
            end
            if ~exist('ref', 'var')
                ref = pipe.io.arrtoij(preRef); % 'ref' object cannot be saved
            end
            src = pipe.io.arrtoij(preData);
            if verbose, fprintf('\nRunning TurboReg... '); end
            trhl.runHL(ref, src, alignstr); % Run turboreg
            tform = trhl.getAllSourcePoints();
            targetgeotransform = trhl.getTargetPoints();
            targetgeotransform = targetgeotransform(1:3, 1:2);
            % Generate 2D transform objects for each frame
            if verbose, fprintf('\nGenerating transform objects  '); end
            for i = 1:NchunkScans % consider parallelizing here!
                ftform = reshape(tform(i,:), 2, 3)';
                reg_transforms{chunkScans(i)} = fitgeotrans(ftform, targetgeotransform, 'affine');
            end
        elseif strcmpi(params.method, 'translation') || strcmpi(params.method, 'rigid') 
            tempTrans = cell(1,NchunkScans);
            parfor i = 1:NchunkScans
                tempTrans{i} = imregtform(preData(:,:,i), preRef, params.method, optimizer, metric ); %fitgeotrans(ftform, targetgeotransform, 'affine'); reg_transforms{chunkScans(i)}
            end
            reg_transforms(chunkScans) = tempTrans;
        else
            fprintf('\nNo valid registration method specified');
        end
        
        % Combine DFT-registration and final registration results (keep in mind we're still in downsampled units)
        if params.prereg > 0
            for i = 1:NchunkScans
                reg_transforms{chunkScans(i)}.T(3,1) = reg_transforms{chunkScans(i)}.T(3,1) + dft_transforms(3,i);
                reg_transforms{chunkScans(i)}.T(3,2) = reg_transforms{chunkScans(i)}.T(3,2) + dft_transforms(4,i);
            end
        end

        % Write out the post-registered movie to tif
        if writeIntermediateTif
            postData = zeros(size(preData)); %zeros(size(rawData));
            for i = 1:NchunkScans 
                postData(:,:,i) = imwarp(preData(:,:,i), reg_transforms{chunkScans(i)}, 'OutputView',imref2d(size(preRef)));
            end
            if verbose, fprintf('\nWriting %s ...', postDataTifPath); end
            WriteTiff(uint16(postData), postDataTifPath); %toc %saveastiff( uint16(postData), regDataTifPath, GrayOpt) ;
        end
        
        % Correct reg_transforms for spatial downsampling
        if params.binXY ~= 1
            for i = 1:NchunkScans
                reg_transforms{chunkScans(i)}.T(3,1) = reg_transforms{chunkScans(i)}.T(3,1)*params.binXY;
                reg_transforms{chunkScans(i)}.T(3,2) = reg_transforms{chunkScans(i)}.T(3,2)*params.binXY;
            end
        end

        % Add registered frames to regData
        if nargout > 1 || writeFinalTif
            for i = 1:NchunkScans 
                unregData(:,:,K+i) = rawData(:,:,i);
                regData(:,:,K+i) = imwarp(rawData(:,:,i), reg_transforms{chunkScans(i)}, 'OutputView',imref2d(size(rawRef)));
            end
            K = K + NchunkScans;
        end
        
        % Save temporary chunk results
        chunkTforms = reg_transforms(chunkScans);
        if verbose,  fprintf('\nSaving %s', chunkSavePath); end
        if c == 1
            save(chunkSavePath, 'chunkScans','chunkTforms','params','rawRef','preRef','targets','targetstr','szstr','alignstr', '-v7.3'); % 'ref' object cannot be saved
        else
            save(chunkSavePath, 'chunkScans','chunkTforms');
        end
    else
        if verbose, fprintf('\nLoading %s', chunkSavePath); end
        load(chunkSavePath);
        reg_transforms(chunkScans) = chunkTforms;
    end
    if Nscan > 1, waitbar(c/Nchunk, w); end
end
if Nscan > 1, delete(w); end 
    
% Generate the final registered movie and write to tiff (optional)
if nargout > 1 || writeFinalTif
    varargout{1} = regData;
    if writeFinalTif
        unregDataTifPath = sprintf('%s%s_unregData.tif', RegTempDir, saveName );
        if verbose, fprintf('\nWriting %s ...', unregDataTifPath); end
        WriteTiff(uint16(unregData), unregDataTifPath);
        regDataTifPath = sprintf('%s%s_regData.tif', RegTempDir, saveName );
        if verbose, fprintf('\nWriting %s ...', regDataTifPath); end
        WriteTiff(uint16(regData), regDataTifPath); %toc %saveastiff( uint16(regData), regDataTifPath, GrayOpt) ;
    end
end
if verbose, fprintf('\n'); end
end