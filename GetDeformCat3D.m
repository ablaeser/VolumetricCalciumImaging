function [deformResultStruct, deformSubStruct, regParams, badInd] = GetDeformCat3D( infoStruct, varargin) % mouse, exptDate, runs
% Get data from registered 3D calcium imaging experiments
IP = inputParser;
addRequired( IP, 'infoStruct', @isstruct)
addOptional( IP, 'deformLim', [], @isstruct)
addParameter( IP, 'zoom', 2, @isnumeric ) 
addParameter( IP, 'edge', nan(1,4), @isnumeric ) %
addParameter( IP, 'basePrct', 10, @isnumeric) %rolling percentile value
addParameter( IP, 'window', 101, @isnumeric) %rolling percentile window
addParameter( IP, 'show', false, @islogical )
addParameter( IP, 'overwrite', false, @islogical )
parse( IP, infoStruct, varargin{:} ); %expt, 
digiZoom = IP.Results.zoom;
edge = IP.Results.edge;
deformLim = IP.Results.deformLim;
basePrct = IP.Results.basePrct; 
windowSize = IP.Results.window;
if mod(windowSize,2) ==0,  windowSize = windowSize+1; end
show = IP.Results.show;
overwrite = IP.Results.overwrite;
savePath =sprintf('%s%s_deformation.mat', infoStruct.dir, infoStruct.exptName);
regParams = struct();
if ~exist(savePath, 'file') || overwrite
    % Determine limits of real deformation (vs TurboReg failed)
    if isempty(deformLim)
        deformLim = struct('trans',[-Inf, Inf], 'scale',[-Inf, Inf], 'stretch',[-Inf,Inf], 'shear',[-Inf, Inf], 'shift',[-Inf, Inf]); %
    end
    
    % Initialize structure to store deformation data
    deformResultStruct = struct('RS_final',nan(infoStruct.totScan, infoStruct.Nplane), 'CS_final',nan(infoStruct.totScan, infoStruct.Nplane),'ZS_final',nan(infoStruct.totScan, infoStruct.Nplane),...
        'transAP',nan(infoStruct.totScan, infoStruct.Nplane), 'transML',nan(infoStruct.totScan, infoStruct.Nplane),...
        'scaleAP',nan(infoStruct.totScan, infoStruct.Nplane), 'scaleML',nan(infoStruct.totScan, infoStruct.Nplane),...
        'shearAP',nan(infoStruct.totScan, infoStruct.Nplane), 'shearML',nan(infoStruct.totScan, infoStruct.Nplane),...
        'shiftZ',nan(infoStruct.totScan, infoStruct.Nplane) ); % 'stretchAP',nan(infoStruct.totScan, infoStruct.Nplane), 'stretchML',nan(infoStruct.totScan, infoStruct.Nplane),

    % Get the results of rigid DFT shifting and z interpolation
    if infoStruct.Nplane > 1
        fprintf('\nGetting results of DFT registration');
        for r = flip(1:infoStruct.Nrun)
            [~,shiftPath] = FileFinder( infoStruct.runDir{r}, 'type','mat', 'contains','dftshifts'); % infoStruct(r).dir
            shiftData(r) = load(shiftPath{1}); %#ok<*AGROW>
        end
        deformResultStruct.RS_final = [shiftData.RS_final]'; 
        deformResultStruct.CS_final = [shiftData.CS_final]';
        deformResultStruct.ZS_final = [shiftData.ZS_final]';

        fprintf('\nGetting results of z interpolation');
        for r = flip(1:infoStruct.Nrun)
            [~,interpPath] = FileFinder( infoStruct.runDir{r}, 'type','mat', 'contains','zinterp'); % infoStruct(r).dir
            %interpPath = FileFind( infoStruct.runDir{r}, 'mat', false, @(x)(contains( x, 'zinterp' )) );
            interpData(r) = load(interpPath{1});
        end
        deformResultStruct.shiftZ = [interpData.ZS_chunk]' + [interpData.ZS1]';
    end

    % Get results of planar affine registration
    fprintf('\nGetting results of planar affine registration');
    regTforms = []; affineParams = [];   
    [~, tformPath] = FileFinder(infoStruct.dir, 'type','mat', 'contains','form'); %FileFind( infoStruct.dir, 'mat', false, @(x)(contains( x, '_affine_tforms' )) );
    if ~isempty(tformPath)
        % Load registration results
        regData = load( tformPath{1} ); %strcat(infoStruct.dir,infoStruct.name,'_affine_tforms.mat')
        if isfield(regData, 'regTform')
            regTforms = regData.regTform;
        elseif isfield(regData, 'tforms_all')
            regTforms = regData.tforms_all;
        end
        % Determine edges used for affine registration
        if isfield(regData, 'params')
            regParams = regData.params;
        end
        if any(isnan(edge))
            if isfield(regData, 'edges')
                edge = regData.edges;
            elseif isfield(regData, 'params') && isfield(regData.params, 'edges')
                edge = regData.params.edges;
            else
                error('Edges unknown');
            end
        end
        fprintf('\nEdges used for affine registration: [L, R, T, B] = [%i, %i, %i, %i]', edge )
        % Get conversion factors
        umPerPixel = (1/0.53)/digiZoom;
        MLpix = infoStruct(1).sz(1)-edge(3)-edge(4);
        APpix = infoStruct(1).sz(2)-edge(1)-edge(2);
        dT = infoStruct.Nplane/infoStruct.framerate;

        zGood = find(~any( cellfun( @isempty, regTforms ), 2 ))';
        if numel(zGood) < infoStruct.Nplane, warning('affine transform data missing for some planes'); end
        for z = zGood 
            for s = 1:infoStruct.totScan
                deformResultStruct.transAP(s,z) = regTforms{z,s}.T(3,1); % Trans AP (right/left +/-) trans_x
                deformResultStruct.transML(s,z) = regTforms{z,s}.T(3,2); % Trans ML (down/up +/-)  trans_y
                deformResultStruct.scaleAP(s,z) = regTforms{z,s}.T(1,1); % Scale AP (inflate/deflate >/< 1)  scale_x
                deformResultStruct.scaleML(s,z) = regTforms{z,s}.T(2,2); % Scale ML (inflate/deflate >/< 1) scale_y
                deformResultStruct.shearAP(s,z) = regTforms{z,s}.T(1,2); % Shear AP (tilt left/right +/-)  shear_x
                deformResultStruct.shearML(s,z) = regTforms{z,s}.T(2,1); % Shear ML (tilt down/right +/-) shear_y
            end
        end
        % Calculate stretch here so it can be used to censor bad datapoints
        deformResultStruct.stretchAP = (1/dT)*[nan(1,infoStruct.Nplane); diff(umPerPixel*APpix./deformResultStruct.scaleAP, 1, 1)];
        deformResultStruct.stretchML = (1/dT)*[nan(1,infoStruct.Nplane); diff(umPerPixel*MLpix./deformResultStruct.scaleML, 1, 1)];
    else
        fprintf('\nNo affine transform data found');
    end
    %{
    % Show the results
    if show
        Nsp_mid = 7; Nsp_right = 9;
        if infoStruct.Nplane == 1
            planeTicks = 1;
            Nsp_mid = 6; Nsp_right = 8;
        elseif infoStruct.Nplane == 15
            planeTicks = [1:3:15, 15];
        elseif infoStruct.Nplane == 30
            planeTicks = [1:5:30, 30];
        end
        runTicks = cumsum(infoStruct.Nscan);
        TL = [0.005,0];
        leftOpt = {[0.01,0.09], [0.06,0.04], [0.05,0.05]};  % {[vert, horz], [bottom, top], [left, right]}
        rightOpt = {[0.01,0.1], [0.06,0.04], [0.05,0.05]};  % {[vert, horz], [bottom, top], [left, right]}
        close all; clearvars sp
        figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'Color','w', 'PaperOrientation','landscape');
        % Results of preliminary registration
        SP(3) = subtightplot(3,3,7,leftOpt{:});
        MakeDeformPlot( deformCatStruct.ZS_final, 'DFT Z shift', TL, runTicks, planeTicks, deformLim.shift )
        % Row
        SP(1) = subtightplot(3,3,1,leftOpt{:});
        MakeDeformPlot( deformCatStruct.RS_final, 'DFT Row shift', TL, runTicks, planeTicks, deformLim.shift )
        title('Rigid DFT Pre-registration');
        % Column
        SP(2) = subtightplot(3,3,4,leftOpt{:});
        MakeDeformPlot( deformCatStruct.CS_final, 'DFT Column shift', TL, runTicks, planeTicks, deformLim.shift )
        linkaxes(SP,'x');
        axis tight;

        % Uncensored Results of affine registration
        % Z Shift
        if infoStruct.Nplane > 1
            sp(7) = subtightplot(Nsp_mid,3,20,rightOpt{:});
            MakeDeformPlot( deformCatStruct.shiftZ, 'Z Shift', TL, runTicks, planeTicks, deformLim.shift )
            xlabel('Volume Scan');
        end
        % Translation
        sp(1) = subtightplot(Nsp_mid,3,2,rightOpt{:});
        MakeDeformPlot( deformCatStruct.transAP, 'AP Translation', TL, runTicks, planeTicks, deformLim.trans )
        title('Turboreg Registration Results');
        sp(2) = subtightplot(Nsp_mid,3,5,rightOpt{:});
        MakeDeformPlot( deformCatStruct.transML, 'ML Translation', TL, runTicks, planeTicks, deformLim.trans )
        % Scaling
        sp(3) = subtightplot(Nsp_mid,3,8,rightOpt{:});
        MakeDeformPlot( deformCatStruct.scaleAP, 'AP Scale', TL, runTicks, planeTicks, deformLim.scale )
        sp(4) = subtightplot(Nsp_mid,3,11,rightOpt{:});
        MakeDeformPlot( deformCatStruct.scaleML, 'ML Scale', TL, runTicks, planeTicks, deformLim.scale )  
        % Shearing
        sp(5) = subtightplot(Nsp_mid,3,14,rightOpt{:});
        MakeDeformPlot( deformCatStruct.shearAP, 'AP Shear', TL, runTicks, planeTicks, deformLim.shear )
        sp(6) = subtightplot(Nsp_mid,3,17,rightOpt{:});
        MakeDeformPlot( deformCatStruct.shearML, 'ML Shear', TL, runTicks, planeTicks, deformLim.shear )
        if infoStruct.Nplane == 1, xlabel('Frame'); end
    end
    %}
    % Find and suppress frames with at least one bad deformation value
    badMat = ...
    (deformResultStruct.transAP < deformLim.trans(1) | deformResultStruct.transAP > deformLim.trans(2)) + (deformResultStruct.transML < deformLim.trans(1) | deformResultStruct.transML > deformLim.trans(2)) + ... 
    (deformResultStruct.scaleAP < deformLim.scale(1) | deformResultStruct.scaleAP > deformLim.scale(2)) + (deformResultStruct.scaleML < deformLim.scale(1) | deformResultStruct.scaleML > deformLim.scale(2)) + ... 
    (deformResultStruct.shearAP < deformLim.shear(1) | deformResultStruct.shearAP > deformLim.shear(2)) + (deformResultStruct.shearML < deformLim.shear(1) | deformResultStruct.shearML > deformLim.shear(2)) + ... 
    (deformResultStruct.stretchAP < deformLim.stretch(1) | deformResultStruct.stretchAP > deformLim.stretch(2)) + (deformResultStruct.stretchML < deformLim.stretch(1) | deformResultStruct.stretchML > deformLim.stretch(2)) + ...    
    (deformResultStruct.shiftZ < deformLim.shift(1) | deformResultStruct.shiftZ > deformLim.trans(2));
    badInd = find(badMat);
    fprintf('\nFound %i bad deformation frames, set to NaN\n', numel(badInd));
    deformCensStruct = deformResultStruct;
    deformCensStruct.transAP(badInd) = NaN;
    deformCensStruct.transML(badInd) = NaN;
    deformCensStruct.scaleAP(badInd) = NaN;
    deformCensStruct.scaleML(badInd) = NaN;
    deformCensStruct.shearAP(badInd) = NaN;
    deformCensStruct.shearML(badInd) = NaN;
    deformCensStruct.shiftZ(badInd) = NaN;

    % Break the results into submovie level structures, and convert to micron units
    deformSubStruct = repmat( struct('dft_RS',[], 'dft_CS',[], 'dft_ZS',[],...
        'transAP',[], 'transML',[], 'transMag',[], 'transAngle',[], 'DtransAP',[], 'DtransML',[], 'DtransMag',[], 'DtransAngle',[],...
        'scaleAP',[], 'scaleML',[], 'scaleMag',[], 'scaleAngle',[], 'stretchAP',[], 'stretchML',[], 'stretchMag',[], 'stretchAngle',[],...
        'shearAP',[], 'shearML',[], 'shearMag',[], 'shearAngle',[], 'DshearAP',[], 'DshearML',[], 'DshearMag',[], 'DshearAngle',[],... 
        'shiftZ',[], 'compressZ',[], 'DshiftZ',[], 'DcompressZ',[]), 1, infoStruct.Nrun );
    scanLims = [0, cumsum(infoStruct.Nscan)];

    for r = 1:infoStruct.Nrun
        subScans = scanLims(r)+1:scanLims(r+1);
        % DFT REGISTRATION RESULTS
        deformSubStruct(r).dft_RS = umPerPixel*deformCensStruct.RS_final(subScans,:);
        deformSubStruct(r).dft_CS = umPerPixel*deformCensStruct.CS_final(subScans,:);
        deformSubStruct(r).dft_ZS = deformCensStruct.ZS_final(subScans);
        
        % Z INTERPOLATION RESULTS
        deformSubStruct(r).shiftZ = deformCensStruct.shiftZ(subScans,:);
        
        % AFFINE REGISTRATION RESULTS (converting to um in the process)
        deformSubStruct(r).transAP = umPerPixel*deformCensStruct.transAP(subScans,:);
        deformSubStruct(r).transML = umPerPixel*deformCensStruct.transML(subScans,:);
        deformSubStruct(r).scaleAP = umPerPixel*APpix./deformCensStruct.scaleAP(subScans,:); %umPerPixel*APpix*deformCatStruct.scaleAP(subScans,:);
        deformSubStruct(r).scaleML = umPerPixel*MLpix./deformCensStruct.scaleML(subScans,:); %umPerPixel*MLpix*deformCatStruct.scaleML(subScans,:);
        deformSubStruct(r).shearAP = deformCensStruct.shearAP(subScans,:);
        deformSubStruct(r).shearML = deformCensStruct.shearML(subScans,:);

        % High-pass filter (optional)
        if basePrct > 0
            fprintf('\nHP Filtering: window = %i scans, %ith percentile subtraction', windowSize, basePrct);
            deformSubStruct(r).transAP = deformSubStruct(r).transAP - MovingPercentile(deformSubStruct(r).transAP, basePrct, windowSize, 'pre');
            deformSubStruct(r).transML = deformSubStruct(r).transML - MovingPercentile(deformSubStruct(r).transML, basePrct, windowSize, 'pre');
            deformSubStruct(r).scaleAP = deformSubStruct(r).scaleAP - MovingPercentile(deformSubStruct(r).scaleAP, basePrct, windowSize, 'pre');
            deformSubStruct(r).scaleML = deformSubStruct(r).scaleML - MovingPercentile(deformSubStruct(r).scaleML, basePrct, windowSize, 'pre');
            deformSubStruct(r).shearAP = deformSubStruct(r).shearAP - MovingPercentile(deformSubStruct(r).shearAP, basePrct, windowSize, 'pre');
            deformSubStruct(r).shearML = deformSubStruct(r).shearML - MovingPercentile(deformSubStruct(r).shearML, basePrct, windowSize, 'pre');
            deformSubStruct(r).shiftZ = deformSubStruct(r).shiftZ - MovingPercentile(deformSubStruct(r).shiftZ, basePrct, windowSize, 'pre');
        end
        % Calculate derivatives
        deformSubStruct(r).compressZ = [nan(numel(subScans),1), diff( deformCensStruct.shiftZ(subScans,:), 1, 2 )]; % spatial derivative
        deformSubStruct(r).DtransAP = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(r).transAP, 1, 1)];
        deformSubStruct(r).DtransML = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(r).transML, 1, 1)];
        deformSubStruct(r).stretchAP = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(r).scaleAP, 1, 1)];
        deformSubStruct(r).stretchML = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(r).scaleML, 1, 1)];
        deformSubStruct(r).DshearAP = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(r).shearAP, 1, 1)];
        deformSubStruct(r).DshearML = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(r).shearML, 1, 1)];
        deformSubStruct(r).DshiftZ = (1/dT)*[nan(1,infoStruct.Nplane); diff( deformSubStruct(r).shiftZ, 1, 1)];
        deformSubStruct(r).DcompressZ = (1/dT)*[nan(1,infoStruct.Nplane); diff( deformSubStruct(r).compressZ, 1, 1)];
        % Convert xy variables to polar coordinates
        [deformSubStruct(r).transAngle, deformSubStruct(r).transMag] = cart2pol( deformSubStruct(r).transAP, deformSubStruct(r).transML ); 
        [deformSubStruct(r).scaleAngle, deformSubStruct(r).scaleMag] = cart2pol( deformSubStruct(r).scaleAP, deformSubStruct(r).scaleML ); 
        [deformSubStruct(r).shearAngle, deformSubStruct(r).shearMag] = cart2pol( deformSubStruct(r).shearAP, deformSubStruct(r).shearML ); 
        [deformSubStruct(r).DtransAngle, deformSubStruct(r).DtransMag] = cart2pol( deformSubStruct(r).DtransAP, deformSubStruct(r).DtransML );
        [deformSubStruct(r).stretchAngle, deformSubStruct(r).stretchMag] = cart2pol( deformSubStruct(r).stretchAP, deformSubStruct(r).stretchML );
        [deformSubStruct(r).DshearAngle, deformSubStruct(r).DshearMag] = cart2pol( deformSubStruct(r).DshearAP, deformSubStruct(r).DshearML);
        % Break up expansive, compressive, and mixed stretch
        expInd = deformSubStruct(r).stretchAngle >= 0 & deformSubStruct(r).stretchAngle <= pi/2;
        compInd = deformSubStruct(r).stretchAngle <= -pi/2 & deformSubStruct(r).stretchAngle >= -pi;
        mixedInd = (deformSubStruct(r).stretchAngle < pi & deformSubStruct(r).stretchAngle > pi/2) | (deformSubStruct(r).stretchAngle < 0 & deformSubStruct(r).stretchAngle > -pi/2);
        deformSubStruct(r).stretchExp = zeros( size(deformSubStruct(r).scaleAP) );
        deformSubStruct(r).stretchExp( expInd ) = deformSubStruct(r).stretchMag( expInd );
        deformSubStruct(r).stretchComp = zeros( size(deformSubStruct(r).scaleAP) );
        deformSubStruct(r).stretchComp( compInd ) = deformSubStruct(r).stretchMag( compInd );
        deformSubStruct(r).stretchMixed = zeros( size(deformSubStruct(r).scaleAP) );
        deformSubStruct(r).stretchMixed( mixedInd ) = deformSubStruct(r).stretchMag( mixedInd );
        %{
        for z = flip(1:infoStruct.Nplane)
            % Primary variables
            [deformSubStruct(r).transAngle(:,z), deformSubStruct(r).transMag(:,z)] = cart2pol( deformSubStruct(r).transAP(:,z), deformSubStruct(r).transML(:,z) ); 
            [deformSubStruct(r).scaleAngle(:,z), deformSubStruct(r).scaleMag(:,z)] = cart2pol( deformSubStruct(r).scaleAP(:,z), deformSubStruct(r).scaleML(:,z) ); 
            [deformSubStruct(r).shearAngle(:,z), deformSubStruct(r).shearMag(:,z)] = cart2pol( deformSubStruct(r).shearAP(:,z), deformSubStruct(r).shearML(:,z) ); 

        end
        
        % Z compression
        deformSubStruct(r).compressZ = [nan(numel(subScans),1), diff( deformCensStruct.shiftZ(subScans,:), 1, 2 )];
        
        % TIME DERIVATIVES
        % Translation
        deformSubStruct(r).DtransAP = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(r).transAP, 1, 1)];
        deformSubStruct(r).DtransML = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(r).transML, 1, 1)];
        deformSubStruct(r).DtransMag = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(r).transMag, 1, 1)];
        deformSubStruct(r).DtransAngle = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(r).transAngle, 1, 1)];
        % Scaling -> stretch
        deformSubStruct(r).stretchAP = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(r).scaleAP, 1, 1)];
        deformSubStruct(r).stretchML = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(r).scaleML, 1, 1)];
        deformSubStruct(r).stretchMag = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(r).scaleMag, 1, 1)];
        deformSubStruct(r).stretchAngle = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(r).scaleAngle, 1, 1)];
        % Shearing
        deformSubStruct(r).DshearAP = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(r).shearAP, 1, 1)];
        deformSubStruct(r).DshearML = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(r).shearML, 1, 1)];
        deformSubStruct(r).DshearMag = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(r).shearMag, 1, 1)];
        deformSubStruct(r).DshearAngle = (1/dT)*[nan(1,infoStruct.Nplane); diff(deformSubStruct(r).shearAngle, 1, 1)];
        %{
        % Scaling
        deformSubStruct(r).stretchAP = nan( infoStruct.Nscan(r), infoStruct.Nplane );
        deformSubStruct(r).stretchAP(2:end,:) = (1/dT)*diff( deformSubStruct(r).scaleAP, 1, 1 );
        deformSubStruct(r).stretchML = nan( infoStruct.Nscan(r), infoStruct.Nplane );
        deformSubStruct(r).stretchML(2:end,:) = (1/dT)*diff( deformSubStruct(r).scaleML, 1, 1 );
        % Shearing
        deformSubStruct(r).DshearAP = nan( infoStruct.Nscan(r), infoStruct.Nplane );
        deformSubStruct(r).DshearAP(2:end,:) = (1/dT)*diff( deformSubStruct(r).shearAP, 1, 1 );
        deformSubStruct(r).DshearML = nan( infoStruct.Nscan(r), infoStruct.Nplane );
        deformSubStruct(r).DshearML(2:end,:) = (1/dT)*diff( deformSubStruct(r).shearML, 1, 1 );
        % Z Shift
        deformSubStruct(r).DshiftZ = nan( infoStruct.Nscan(r), infoStruct.Nplane );
        deformSubStruct(r).DshiftZ(2:end,:) = (1/dT)*diff( deformSubStruct(r).shiftZ, 1, 1 );
        %}
        % Z Shift
        deformSubStruct(r).DshiftZ = (1/dT)*[nan(1,infoStruct.Nplane); diff( deformSubStruct(r).shiftZ, 1, 1)];
        deformSubStruct(r).DcompressZ = (1/dT)*[nan(1,infoStruct.Nplane); diff( deformSubStruct(r).compressZ, 1, 1)];
        %}
    end
    %{
    if show
        % Censored final results
        % Z shift
        if infoStruct.Nplane > 1
            sp(16) = subtightplot(Nsp_right,3,27,rightOpt{:});
            MakeDeformPlot( vertcat(deformSubStruct.shiftZ), 'Z Shift (planes)', TL, runTicks, planeTicks, deformLim.shift )
            xlabel('Volume Scan');
        end
        % Translation
        sp(8) = subtightplot(Nsp_right,3,3,rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.transAP), 'AP Trans (um)', TL, runTicks, planeTicks, [-Inf,Inf] ) % deformLim.trans
        title('Final, Censored Deformation');
        sp(9) = subtightplot(Nsp_right,3,6,rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.transML), 'ML Trans (um)', TL, runTicks, planeTicks, [-Inf,Inf] )
        % Scaling
        sp(10) = subtightplot(Nsp_right,3,9,rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.scaleAP), 'AP Scale (um)', TL, runTicks, planeTicks, [-Inf,Inf] ) % deformLim.scale
        sp(11) = subtightplot(Nsp_right,3,12,rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.scaleML), 'ML Scale (um)', TL, runTicks, planeTicks, [-Inf,Inf] )  
        % Stretch
        sp(12) = subtightplot(Nsp_right,3,15,rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.stretchAP), 'AP Stretch (um/s)', TL, runTicks, planeTicks, [-Inf,Inf] ) % deformLim.stretch
        sp(13) = subtightplot(Nsp_right,3,18,rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.stretchML), 'ML Stretch (um/s)', TL, runTicks, planeTicks, [-Inf,Inf] )
        % Shearing
        sp(14) = subtightplot(Nsp_right,3,21,rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.shearAP), 'AP Shear', TL, runTicks, planeTicks, [-Inf,Inf] )
        sp(15) = subtightplot(Nsp_right,3,24,rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.shearML), 'ML Shear', TL, runTicks, planeTicks, [-Inf,Inf] )
        if infoStruct.Nplane == 1
            linkaxes(sp,'x');
            xlabel('Frame');
        else
            linkaxes(sp,'xy');
            impixelinfo;
        end
        axis tight;
    end
    %}
    fprintf('\nSaving %s', savePath);
    save(savePath, 'infoStruct','deformLim','deformSubStruct','deformResultStruct','deformCensStruct','badInd','expInd','compInd','mixedInd',...
        'affineParams','tformPath','dT','APpix','MLpix','edge','basePrct','windowSize','digiZoom','umPerPixel','zGood','scanLims');
else
    fprintf('\nLoading %s', savePath);
    load(savePath);
end

if show
    runTicks = cumsum(infoStruct.Nscan);
    timeTicks = 1:5000:infoStruct.totScan;
    TL = [0.005,0];
    leftOpt = {[0.01,0.09], [0.06,0.04], [0.05,0.05]};  % {[vert, horz], [bottom, top], [left, right]}
    rightOpt = {[0.01,0.1], [0.06,0.04], [0.05,0.05]};  % {[vert, horz], [bottom, top], [left, right]}
    close all; clearvars sp
    figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'Color','w', 'PaperOrientation','landscape');
    if infoStruct.Nplane == 1
        planeTicks = 1;
        % Translation 
        sp(1) = subtightplot(8, 2, 1, rightOpt{:});
        MakeDeformPlot( deformResultStruct.transAP, 'AP Translation', TL, runTicks, planeTicks, deformLim.trans )
        title('Turboreg Registration Results');
        sp(2) = subtightplot(8, 2, 3, rightOpt{:});
        MakeDeformPlot( deformResultStruct.transML, 'ML Translation', TL, runTicks, planeTicks, deformLim.trans )
        % Scaling
        sp(3) = subtightplot(8, 2, 5, rightOpt{:});
        MakeDeformPlot( deformResultStruct.scaleAP, 'AP Scale', TL, runTicks, planeTicks, deformLim.scale )
        sp(4) = subtightplot(8, 2, 7, rightOpt{:});
        MakeDeformPlot( deformResultStruct.scaleML, 'ML Scale', TL, runTicks, planeTicks, deformLim.scale )  
        % Stretch
        sp(5) = subtightplot(8, 2, 9, rightOpt{:});
        MakeDeformPlot( deformResultStruct.stretchAP, 'AP Stretch', TL, runTicks, planeTicks, deformLim.stretch )
        sp(6) = subtightplot(8, 2, 11, rightOpt{:});
        MakeDeformPlot( deformResultStruct.stretchML, 'ML Stretch', TL, timeTicks, planeTicks, deformLim.stretch )
        % Shearing
        sp(7) = subtightplot(8, 2, 13, rightOpt{:});
        MakeDeformPlot( deformResultStruct.shearAP, 'AP Shear', TL, runTicks, planeTicks, deformLim.shear )
        sp(8) = subtightplot(8, 2, 15, rightOpt{:});
        MakeDeformPlot( deformResultStruct.shearML, 'ML Shear', TL, timeTicks, planeTicks, deformLim.shear ) % runTicks
        xlabel('Frame');
        
        % Censored final results
        % Translation
        sp(9) = subtightplot(8, 2, 2, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.transAP), 'AP Trans (um)', TL, runTicks, planeTicks, [-Inf,Inf] ) % deformLim.trans
        title('Final, Censored Deformation');
        sp(10) = subtightplot(8, 2, 4, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.transML), 'ML Trans (um)', TL, runTicks, planeTicks, [-Inf,Inf] )
        % Scaling
        sp(11) = subtightplot(8, 2, 6, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.scaleAP), 'AP Scale (um)', TL, runTicks, planeTicks, [-Inf,Inf] ) % deformLim.scale
        sp(12) = subtightplot(8, 2, 8, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.scaleML), 'ML Scale (um)', TL, runTicks, planeTicks, [-Inf,Inf] )  
        % Stretch
        sp(13) = subtightplot(8, 2, 10, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.stretchAP), 'AP Stretch (um/s)', TL, runTicks, planeTicks, [-Inf,Inf] ) % deformLim.stretch
        sp(14) = subtightplot(8, 2, 12, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.stretchML), 'ML Stretch (um/s)', TL, runTicks, planeTicks, [-Inf,Inf] )
        % Shearing
        sp(15) = subtightplot(8, 2, 14, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.shearAP), 'AP Shear', TL, runTicks, planeTicks, [-Inf,Inf] )
        sp(16) = subtightplot(8, 2, 16, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.shearML), 'ML Shear', TL, runTicks, planeTicks, [-Inf,Inf] )
        
        linkaxes(sp,'x');
        xlabel('Frame');
        axis tight;
    else
        if infoStruct.Nplane == 15
            planeTicks = [1:3:15, 15];
        elseif infoStruct.Nplane == 30
            planeTicks = [1:5:30, 30];
        else
            planeTicks = round(linspace(1, infoStruct.Nplane, infoStruct.Nplane/4));
        end
        % Results of preliminary registration
        SP(3) = subtightplot(3,3,7,leftOpt{:});
        MakeDeformPlot( deformResultStruct.ZS_final, 'DFT Z shift', TL, runTicks, planeTicks, deformLim.shift )
        % Row
        SP(1) = subtightplot(3,3,1,leftOpt{:});
        MakeDeformPlot( deformResultStruct.RS_final, 'DFT Row shift', TL, runTicks, planeTicks, [-Inf,Inf] ) % deformLim.shift
        title('Rigid DFT Pre-registration');
        % Column
        SP(2) = subtightplot(3,3,4,leftOpt{:});
        MakeDeformPlot( deformResultStruct.CS_final, 'DFT Column shift', TL, runTicks, planeTicks, [-Inf,Inf] ) % deformLim.shift
        linkaxes(SP,'x');
        axis tight;
        
        % Uncensored Results of affine registration
        % Z Shift
        sp(9) = subtightplot(7, 3, 20, rightOpt{:});
        MakeDeformPlot( deformResultStruct.shiftZ, 'Z Shift', TL, runTicks, planeTicks, deformLim.shift )
        xlabel('Volume Scan');
        % Translation 
        sp(1) = subtightplot(7, 3, 2, rightOpt{:});
        MakeDeformPlot( deformResultStruct.transAP, 'AP Translation', TL, runTicks, planeTicks, deformLim.trans )
        title('Turboreg Registration Results');
        sp(2) = subtightplot(7, 3, 5, rightOpt{:});
        MakeDeformPlot( deformResultStruct.transML, 'ML Translation', TL, runTicks, planeTicks, deformLim.trans )
        % Scaling
        sp(3) = subtightplot(7, 3, 8, rightOpt{:});
        MakeDeformPlot( deformResultStruct.scaleAP, 'AP Scale', TL, runTicks, planeTicks, deformLim.scale )
        sp(4) = subtightplot(7, 3, 11, rightOpt{:});
        MakeDeformPlot( deformResultStruct.scaleML, 'ML Scale', TL, runTicks, planeTicks, deformLim.scale ) 
        % Stretch
        sp(5) = subtightplot(7, 3, 9, rightOpt{:});
        MakeDeformPlot( deformResultStruct.stretchAP, 'AP Stretch', TL, runTicks, planeTicks, deformLim.stretch )
        sp(6) = subtightplot(7, 3, 12, rightOpt{:});
        MakeDeformPlot( deformResultStruct.stretchML, 'ML Stretch', TL, timeTicks, planeTicks, deformLim.stretch )
        % Shearing
        sp(7) = subtightplot(7, 3, 14, rightOpt{:});
        MakeDeformPlot( deformResultStruct.shearAP, 'AP Shear', TL, runTicks, planeTicks, deformLim.shear )
        sp(8) = subtightplot(7, 3, 17, rightOpt{:});
        MakeDeformPlot( deformResultStruct.shearML, 'ML Shear', TL, runTicks, planeTicks, deformLim.shear )

        % Censored final results
        % Z shift
        sp(18) = subtightplot(9, 3, 27, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.shiftZ), 'Z Shift (planes)', TL, runTicks, planeTicks, deformLim.shift )
        xlabel('Volume Scan');
        % Translation
        sp(10) = subtightplot(9, 3, 3, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.transAP), 'AP Trans (um)', TL, runTicks, planeTicks, [-Inf,Inf] ) % deformLim.trans
        title('Final, Censored Deformation');
        sp(11) = subtightplot(9, 3, 6, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.transML), 'ML Trans (um)', TL, runTicks, planeTicks, [-Inf,Inf] )
        % Scaling
        sp(12) = subtightplot(9, 3, 9, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.scaleAP), 'AP Scale (um)', TL, runTicks, planeTicks, [-Inf,Inf] ) % deformLim.scale
        sp(13) = subtightplot(9, 3, 12, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.scaleML), 'ML Scale (um)', TL, runTicks, planeTicks, [-Inf,Inf] )  
        % Stretch
        sp(14) = subtightplot(9, 3, 15, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.stretchAP), 'AP Stretch (um/s)', TL, runTicks, planeTicks, [-Inf,Inf] ) % deformLim.stretch
        sp(15) = subtightplot(9, 3, 18, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.stretchML), 'ML Stretch (um/s)', TL, runTicks, planeTicks, [-Inf,Inf] )
        % Shearing
        sp(16) = subtightplot(9, 3, 21, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.shearAP), 'AP Shear', TL, runTicks, planeTicks, [-Inf,Inf] )
        sp(17) = subtightplot(9, 3, 24, rightOpt{:});
        MakeDeformPlot( vertcat(deformSubStruct.shearML), 'ML Shear', TL, runTicks, planeTicks, [-Inf,Inf] )

        linkaxes(sp,'xy');
        impixelinfo;
        axis tight;
    end
end
end