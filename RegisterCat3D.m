function RegisterCat3D(sbxInfo, varargin) % mouse, exptDate
% Run scanbox data through stages of rigid 3D registration, z interpolation and affine registration. Optionally, write out tifs showing results from each stage.
IP = inputParser;
addRequired( IP, 'sbxInfo', @isstruct )
addOptional( IP, 'regParams', struct(), @isstruct )
addParameter( IP, 'overwrite', false, @islogical )
addParameter( IP, 'edge', [80,80,20,20], @isnumeric ) % [60,60,40,40]
addParameter( IP, 'chunk', 20, @isnumeric ) %don't go over 20
addParameter( IP, 'writeChan', 'both', @ischar )
addParameter( IP, 'refChan', 'green', @ischar ) % for scanbox, 1 = green, 2 = red. -1 = both
addParameter( IP, 'scale', 2, @isnumeric ) %
addParameter( IP, 'minInt', 1500, @isnumeric )
addParameter( IP, 'fix', false, @islogical )
addParameter( IP, 'flip', true, @islogical )
addParameter( IP, 'preaff', true, @islogical )
addParameter( IP, 'binT', NaN, @isnumeric )
% Tif writing parameters
addParameter( IP, 'writeZ', true, @islogical )
addParameter( IP, 'Zint', 3, @isnumeric )
addParameter( IP, 'Zwrite', [], @isnumeric )
parse( IP, sbxInfo, varargin{:} ); % mouse, exptDate,
regParams = IP.Results.regParams;
regToggle = ~isempty( fieldnames( regParams ) );
fixSbx = IP.Results.fix;
flipZ = IP.Results.flip;
preaffToggle = IP.Results.preaff;
binT = IP.Results.binT;
setEdges = IP.Results.edge;
%fixedEdges = IP.Results.edge; 
scaleFactor = IP.Results.scale;
chunkSize = IP.Results.chunk;
minInt = IP.Results.minInt;
writeChan = IP.Results.writeChan;
refChan = IP.Results.refChan;  %1; %1 = red, 2 = green
refChanInd = find(contains({'red','green'}, refChan)); %find(contains({'green','red'}, refChan));
%[refPMT, refPMTname] = DeterminePMT(refChan, sbxInfo); % PMT1 = green, PMT2 = red
overwrite = IP.Results.overwrite;
writeZ = IP.Results.writeZ;
Zint = IP.Results.Zint;
Zwrite = IP.Results.Zwrite;

% Determine file paths and names
fSep = '\';
sbxInputPath = sbxInfo.path; %sbxPath{runs};
ZtifDir = strcat(sbxInfo.dir, 'Ztifs', fSep); mkdir( ZtifDir );
[fDir, fName] = fileparts(sbxInputPath);
pathTemplate = strcat( fDir, fSep, fName );
sbxFixPath = [pathTemplate, '.sbxfix'];
sbxOptPath = [pathTemplate, '.sbxopt'];
sbxDftPath = [pathTemplate, '.sbxdft'];
sbxZpath = [pathTemplate, '.sbxz'];
interpMatPath = strcat(sbxInfo.dir,sbxInfo.exptName,'_zinterp.mat');
shiftPath = strcat(pathTemplate, '_dftshifts.mat');
rawGreenProjPath = strcat(pathTemplate,'_rawProj_green.tif');
rawRedProjPath = strcat(pathTemplate,'_rawProj_red.tif');
rawProjPath = strcat(pathTemplate,'_rawProj.tif');
%fixProjPath = strcat(pathTemplate,'_rawProj_red.tif');
optProjPath = strcat(pathTemplate,'_optProj.tif');
dftProjPath = strcat(pathTemplate,'_dftProj.tif');
interpProjPath = strcat(pathTemplate,'_interpProj.tif');
if ~isempty(fieldnames(regParams)) %isempty(regParams)
    if isempty(regParams.name)
        regNameStr = '';
    else
        regNameStr = ['_',regParams.name];
    end
    sbxRegPath = sprintf('%s%s.sbxreg', pathTemplate, regNameStr);
    regProjPath = sprintf('%s%s_regProj.tif', pathTemplate, regNameStr); %strcat(pathTemplate,'_regProj_RGB.tif');
    regGreenProjPath = strcat(pathTemplate,'_regProj_green.tif'); % strcat(pathTemplate,'_regProj_green.tif');
    regRedProjPath = strcat(pathTemplate,'_regProj_red.tif'); %sprintf('%s%s.tif', pathTemplate, regNameStr);
else
    sbxRegPath = '';
    regProjPath = '';
    regGreenProjPath = '';
    regRedProjPath = '';
end

if rem(sbxInfo.totScan, chunkSize) ~= 0
    chunkSize = max(factor(sbxInfo.totScan));
    fprintf('Non-integer number of chunks, changed chunk size to %i', chunkSize);
end
Nchunk = sbxInfo.totScan/chunkSize;
if isnan(binT) && sbxInfo.totScan > 5000
    binT = 155;
else
    binT = 1;
end
chanName = {'red','green'};
fprintf('\nProcessing %s: %i frames, %i planes, %i scans. %i channel. temporal binning = %i, chunkSize = %i (%2.1f chunks)\n', sbxInfo.exptName, sbxInfo.nframes, sbxInfo.Nplane, sbxInfo.totScan, sbxInfo.nchan, binT, chunkSize, Nchunk)
% Run pre-registration processing: correct optotune warping, perform rigid registration in 3D, and interpolate Z planes (optional)
if preaffToggle
    if sbxInfo.Nplane > 1  && fixSbx
        if (~exist(sbxFixPath,'file') || overwrite)
            fprintf('\n   Correcting sbx z order... ');
            sbxInfo = FixSBX(sbxInputPath, sbxInfo, flipZ, overwrite);
        end
        sbxInputPath = sbxFixPath;
    end
    
    % Write tifs from the raw data
    if ((~exist(rawProjPath,'file') && ~exist(rawGreenProjPath,'file') && ~exist(rawRedProjPath,'file')) || overwrite)
        fprintf('\n   Writing raw projection stack(s)... ');
        if sbxInfo.Nplane == 1
            [~, rawChan] = WriteSbxPlaneTif(sbxInputPath, sbxInfo, 1, 'dir',ZtifDir, 'name',fName, 'type','raw', 'binT',binT, 'verbose',true, 'chan','both', 'overwrite',overwrite, 'rescale',true );
            for chan = find(~cellfun(@isempty, rawChan))
                rawProjPath = sprintf('%s_rawProj_%s.tif', pathTemplate, chanName{chan} );
                WriteTiff(uint16(mean(rawChan{chan}, 3)), rawProjPath); %pipe.io.writeTiff( uint16(mean(rawChan{chan}, 3)), rawProjPath);
            end
        else
            WriteSbxProjection(sbxInputPath, sbxInfo, rawProjPath, 'binT',binT, 'verbose',true, 'chan',writeChan, 'dir',[ZtifDir,'Raw\'], 'name',sbxInfo.exptName);
        end
    end
    
    if sbxInfo.Nplane > 1
        % lensing correction for neurolabware data
        if ~exist(sbxOptPath,'file') || overwrite
            fprintf('\n   Correcting Neurolabware data... ');
            GetOptotuneWarp(sbxInputPath, sbxInfo, 'chan',refChan, 'type','rigid', 'edges',setEdges, 'save',true);  % , 'show',true , 'scale',scaleFactor reg , 'firstRefScan',500
        end
        if (~exist(optProjPath,'file') || overwrite)
            fprintf('\n   Writing sbxopt projection stack... ');
            WriteSbxProjection(sbxOptPath, sbxInfo, optProjPath, 'monochrome',false, 'RGB',false); % , 'dir',[ZtifDir,'Opt\'], 'name',[sbxInfo.exptName,'_Opt']
        end
        
        % 3D DFT shifts (rigid)
        optMean = loadtiff( optProjPath );
        if isempty(setEdges)
            if sbxInfo.nchan == 1
                optEdges = GetEdges( optMean(:,:,end), 'minInt',minInt, 'show',true ); % fprintf('edges = [%i, %i, %i, %i]\n', interpEdges )
            else
                optEdges = GetEdges( optMean(:,:,refChanInd,end), 'minInt',minInt/2^8, 'show',true  );
            end
        else
            optEdges = setEdges;
            if sbxInfo.nchan == 1
                ShowEdges(edges, optMean(:,:,end))
            else
                ShowEdges(setEdges, optMean(:,:,refChanInd,end))
            end
        end
        
        % Calculate rigid corrections
        if ~exist(shiftPath,'file') || overwrite
            fprintf('\n   Calculating  3D DFT shifts... ');
            CorrectData3D(sbxOptPath, sbxInfo, shiftPath, refChan, 'chunkSize',chunkSize, 'edges',optEdges, 'scale',scaleFactor);  
        end
        % make registered SBX file
        if ~exist(sbxDftPath,'file') || overwrite
            fprintf('\n   Writing sbxdft... ');
            MakeSbxDFT(sbxOptPath, sbxInfo, shiftPath, refChan, 'edges',optEdges, 'proj',true); % , 'zprojPath',zprojPath zproj_mean =
        end
        % write mean projection of rigid-corrected data
        if ~exist(dftProjPath,'file') || overwrite
            fprintf('\n   Writing sbxdft projection stack... ');
            WriteSbxProjection(sbxDftPath, sbxInfo, dftProjPath, 'dir',[ZtifDir,'DFT\'], 'name',[sbxInfo.exptName,'_dft'], 'monochrome',false, 'RGB',false); % , 'dir',[ZtifDir,'zInterp\'], 'name',[sbxInfo.exptName,'_zinterp']
        end
        
        % z interpolation
        dftMean = loadtiff( dftProjPath );
        if isempty(setEdges)
            if sbxInfo.nchan == 1
                dftEdges = GetEdges( dftMean(:,:,end), 'minInt',minInt, 'show',true ); % fprintf('edges = [%i, %i, %i, %i]\n', interpEdges )
            else
                dftEdges = GetEdges( dftMean(:,:,refChanInd,end), 'minInt',minInt/2^8, 'show',true  );
            end
        else
            dftEdges = setEdges;
            if sbxInfo.nchan == 1
                ShowEdges(edges, dftMean(:,:,end))
            else
                ShowEdges(setEdges, dftMean(:,:,refChanInd,end))
            end
        end
        
        if ~exist(sbxZpath,'file') || overwrite
            if ~exist(interpMatPath,'file')
                InterpZ(sbxDftPath, sbxInfo, interpMatPath, 'scale',scaleFactor, 'edges',dftEdges, 'chunkSize',chunkSize); %DFT_reg_z_interp(sbxDftPath, interpMatPath, refChan, scaleFactor, Nchunk, 'optotune',false, 'edges',dftEdges);
            end
            MakeSbxZ(sbxDftPath, sbxInfo, interpMatPath); % SBX_z_interp(sbxDftPath, interpMatPath);
        end
        if ~exist(interpProjPath, 'file') || overwrite
            fprintf('\n   Writing interpolated projection stack... ');
            interpMean = WriteSbxProjection(sbxZpath, sbxInfo, interpProjPath, 'dir',[ZtifDir,'zInterp\'], 'name',[sbxInfo.exptName,'_zinterp'], 'monochrome',false, 'RGB',false);
        end
    end
else
    sbxZpath = sbxInputPath;
end

% registration transform on each plane
if (regToggle && ~exist(sbxRegPath,'file')) || overwrite %true %
    if isfield(sbxInfo, 'scanLim')
        if isempty( regParams.refScan )
            if sbxInfo.Nrun <= 2
                regParams.refRun = 1;
            else
                regParams.refRun = 2;
            end
            regParams.refScan = sbxInfo.scanLim(regParams.refRun)+2:sbxInfo.scanLim(regParams.refRun+1)-2;
        end
    else
        regParams.refRun = NaN;
        regParams.refScan = [];
    end
    fprintf('\n   Performing planar registration... '); % (reference averaged over scans %i - %i)
    AlignPlanes( sbxInfo.path, sbxInfo, regParams, 'overwrite',overwrite, 'outPath',sbxRegPath );
end
if exist(sbxRegPath, 'file')
    if (~exist(regProjPath, 'file') || overwrite)
        WriteSbxProjection(sbxRegPath, sbxInfo, regProjPath, 'chan','both', 'overwrite',overwrite, 'dir',ZtifDir, 'name',fName, 'type','reg', 'monochrome',true, 'RGB',false);
    end
    if (~exist(regGreenProjPath, 'file') || overwrite)
        WriteSbxProjection(sbxRegPath, sbxInfo, regGreenProjPath, 'chan','green', 'overwrite',overwrite, 'monochrome',false, 'RGB',false);
    end
    if (~exist(regRedProjPath, 'file') || overwrite)
        WriteSbxProjection(sbxRegPath, sbxInfo, regRedProjPath, 'chan','red', 'overwrite',overwrite, 'monochrome',false, 'RGB',false);
    end
end
%{
if exist(sbxRegPath,'file') && ((~exist(regGreenProjPath,'file') && ~exist(regRedProjPath,'file')) || overwrite) %(~exist(regProjPath, 'file') || overwrite)
    if sbxInfo.Nplane == 1
        [~,regChan] = WriteSbxPlaneTif(sbxRegPath, sbxInfo, 1, 'dir',ZtifDir, 'name',sprintf('%s%s', fName, regNameStr), 'type','reg', 'binT',binT, 'verbose',true, 'chan','both', 'overwrite',overwrite );
        for chan = find(~cellfun(@isempty, regChan))
            regProjPath = sprintf('%s_regProj_%s.tif', pathTemplate, chanName{chan} );
            pipe.io.writeTiff( uint16(mean(regChan{chan}, 3)), regProjPath);
        end
    end
end
%}

% Write single-plane tifs from different stages of the registration
if writeZ && sbxInfo.Nplane > 1
    fprintf('\n   Writing single plane tifs... ');
    if sbxInfo.Nrun <= 2
        scaleFactor = 1;
    elseif sbxInfo.Nrun > 2 && sbxInfo.Nrun <= 4
        scaleFactor = 2;
    else
        scaleFactor = 4;
    end
    if isempty(Zwrite), Zwrite = unique([1:Zint:sbxInfo.Nplane, sbxInfo.Nplane]);  end
    % Write raw data tifs
    % {
    for z = 1:numel(Zwrite)
        WriteSbxPlaneTif(sbxInputPath, sbxInfo, Zwrite(z), 'dir',ZtifDir, 'name',sbxInfo.exptName, 'type','raw', ...
            'edge',[0,0,0,0], 'scale',scaleFactor, 'verbose',true, 'chan','both', 'binT',29 ); % , 'binT',4
        %WriteSbxPlaneTif(sbxCatPath, Zwrite(z), 'dir',tifDir, 'name',sbxInfo.exptName, 'type','raw', 'edge',[0,0,0,0], 'scale',scaleFactor, 'verbose',true, 'chan',1 );
    end
    % DFT tifs
    if exist(dftProjPath, 'file')
        dftMean = loadtiff( dftProjPath );
        if sbxInfo.nchan == 1
            dftEdges = GetEdges( dftMean(:,:,end), 'minInt',minInt, 'show',true );
        else
            dftEdges = GetEdges( dftMean(:,:,1,end), 'minInt',minInt, 'show',true );
        end
        close all;
        for z = 1:numel(Zwrite) %Zwrite %[1,30] %[6,14,19,29]%
            WriteSbxPlaneTif(sbxDftPath, Zwrite(z), 'dir',ZtifDir, 'name',sbxInfo.exptName, 'type','dft', 'edge',dftEdges, 'scale',scaleFactor, 'verbose',true, 'chan',1 );
        end
    end
    
    % Interp Tifs
    if exist(interpProjPath, 'file')
        interpMean = loadtiff( interpProjPath );
        if sbxInfo.nchan == 1
            interpEdges = GetEdges( interpMean(:,:,end), 'minInt',minInt, 'show',true );
        else
            interpEdges = GetEdges( interpMean(:,:,1,end), 'minInt',minInt/2^8, 'show',true );
        end
        close all;
        for z = 1:numel(Zwrite) %Zwrite %[1,30] %[6,14,19,29]%
            WriteSbxPlaneTif(sbxZpath, sbxInfo, Zwrite(z), 'dir',ZtifDir, 'name',sbxInfo.exptName, 'type','interp', 'edge',interpEdges, 'scale',scaleFactor, 'verbose',true, 'chan','both');
        end
    end
    %}
    % Write reg tifs
    if exist(sbxRegPath, 'file') && regToggle
        for z = 1:numel(Zwrite)
            WriteSbxPlaneTif(sbxRegPath, sbxInfo, Zwrite(z), 'dir',ZtifDir, 'name',fName, 'type','reg', 'edge',regParams.edges, 'scale',scaleFactor, 'binT',binT, 'verbose',true, 'chan','both' );
        end
    end
end
end