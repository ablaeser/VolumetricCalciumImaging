clear; clc; close all;
dataTablePath = 'D:\MATLAB\Movie Processing\ExperimentIndex.xlsx'; %'C:\Users\ablaeser\OneDrive\OneDrive - Beth Israel Lahey Health\Andrew Blaeser\MATLAB\Dura\3D\3DduraData.xlsx'; % 
dataDir = 'D:\2photon\';  %'C:\2photon';
% Parse data table
dataCol = struct('mouse',1, 'date',2, 'FOV',3, 'volume',4, 'run',5, 'done',6);
dataTable = readcell( dataTablePath );  dataTable(1,:) = [];
Nexpt = size(dataTable, 1);
dataTable(:,dataCol.date) = cellfun(@num2str, dataTable(:,dataCol.date), 'UniformOutput',false);
missingData = cellfun(@all, cellfun(@ismissing, dataTable, 'UniformOutput',false), 'UniformOutput',true );
% Initialize variables
expt = repmat( struct('mouse','', 'date',[], 'dir','', 'name','', 'runs',[], 'Nruns',NaN), Nexpt, 1); % , 'proc',[]
deformLim = struct('trans',[-20, 20], 'scale',[0.95, 1.05], 'stretch',100*[-1,1], 'shear',[-0.03, 0.03], 'shift',[-3.5, 3.5]);
runInfo = cell(1,Nexpt); Tscan = cell(1,Nexpt); segParams = cell(1,Nexpt); ROI = cell(1,Nexpt); 
% Determine which experiments to analyze
xPresent = 1; % 
Npresent = numel(xPresent);
overwriteROI = false;
for x = xPresent % 51
    % Parse data table
    expt(x).mouse = dataTable{x,dataCol.mouse};
    expt(x).date = dataTable{x,dataCol.date};  if isnumeric(expt(x).date), expt(x).date = num2str(expt(x).date); end
    if missingData(x,dataCol.FOV)
        expt(x).fov = 1;
        expt(x).dir = sprintf('%s%s\\%s\\', dataDir, expt(x).mouse, expt(x).date); %sprintf('%s%s\\%s_%s\\', dataDir, expt(x).mouse, expt(x).date, expt(x).mouse);
        expt(x).name = sprintf('%s_%s', expt(x).mouse, expt(x).date);
    else
        expt(x).fov = dataTable{x,dataCol.FOV};
        expt(x).dir = sprintf('%s%s\\%s_FOV%i\\', dataDir, expt(x).mouse, expt(x).date, expt(x).fov); % sprintf('%s%s\\%s_FOV%i_%s\\', dataDir, expt(x).mouse, expt(x).date, expt(x).fov, expt(x).mouse);
        expt(x).name = sprintf('%s_%s_FOV%i', expt(x).mouse, expt(x).date, expt(x).fov);
    end
    fprintf('\n\n%s', expt(x).name )
    if isnumeric(dataTable{x,dataCol.run}), expt(x).runs = dataTable{x,dataCol.run}; else, expt(x).runs = str2num( dataTable{x,dataCol.run} ); end
    expt(x).Nruns = numel(expt(x).runs); expt(x).Nroi = NaN; expt(x).Naxon = NaN;
    
    % Get basic run-level metadata and data, and maybe perform registration on individual runs
    for r = flip(expt(x).runs)
        runInfo{x}(r) = MakeInfoStruct( dataDir, expt(x).mouse, expt(x).date, r, expt(x).fov );
        RegisterCat3D( runInfo{x}(r), 'overwrite',false, 'writeZ',false, 'Zint',3, 'chunk',20, 'minInt',1500, 'refChan','green', 'preaff',true, 'fix',false); 
    end
    [Tscan{x}, runInfo{x}] = GetTime(runInfo{x}); 

    % Concatenate metdata and single plane movies   
    expt(x).Nrow = runInfo{x}(1).sz(1); expt(x).Ncol = runInfo{x}(1).sz(2); expt(x).Nplane = runInfo{x}(1).otlevels; expt(x).Nchan = runInfo{x}(1).nchan;
    expt(x).Nscan = floor([runInfo{x}.nframes]/expt(x).Nplane); expt(x).totScan = sum(expt(x).Nscan); expt(x).totFrame = sum([runInfo{x}.totFrame]); 
    expt(x).scanLims = [0, cumsum(expt(x).Nscan)];
    expt(x).frameRate = runInfo{x}(1).framerate; expt(x).scanRate = runInfo{x}(1).framerate/expt(x).Nplane; 
    if expt(x).Nplane > 1
        expt(x).sbx = strcat(expt(x).dir, expt(x).name, '.sbx_interp');
    else
        expt(x).sbx = strcat(expt(x).dir, expt(x).name, '.sbx_affine');
    end
    catInfo(x) = ConcatenateRunInfo(expt(x), runInfo{x}, 'suffix','sbxcat', 'overwrite',false); % Get concatenated metadata
    expt(x).zoom = str2double(catInfo(x).config.magnification_list(catInfo(x).config.magnification,:));
    %expt(x).umPerPixel = (1/0.53)/expt(x).zoom; 
   
    % Concatenate runs into a single sbx file
    interRunShift = ConcatenateExptRuns(expt(x), runInfo{x}, catInfo(x), 'refRun',2 );
    %ConcatenateSinglePlane(expt(x), runInfo{x}, catInfo(x), 2, [100,100,40,40]); %interRunShift =
    
    % Run registration on concatenated data
    % {
    regParams.refChan = 'green';  
    regParams.refRun = 1; 
    regParams.refScan = []; 
    regParams.chunkSize = 1;
    regParams.histmatch = false; %true; %  
    regParams.avgT = 0;
    regParams.binXY = 2;
    regParams.binT = 1;
    regParams.prereg = false;
    regParams.highpass = 0;
    regParams.lowpass = 0;
    regParams.medFilter = [0,0,0];
    regParams.minInt = 1500;
    regParams.edges = [80,80,20,20]; % []; %  [80,90,20,20]; %[60,30,20,20]; % [80,90,20,20]; % [20,60,20,20]; 
    regParams.name = '';
    regParams.method = 'rigid'; % 'affine';
    RegisterCat3D( catInfo(x), regParams, 'overwrite',false, 'writeZ',true, 'Zint',3, 'chunk',30, 'preaff',false); % expt(x).mouse, expt(x).date, expt(x).runs , 'cat',false
    
    % Get deformation data
    [deformCat{x}, deform{x}, regParams, badInd] = GetDeformCat3D( catInfo(x), deformLim, 'show',true, 'overwrite',false, 'window',find(Tscan{x}{1}<=32,1,'last') );  %  deformCat, affParams
    
    % Segment movies
    %{
    zProj = 3:14;    
    badMat = false( size(deformCat{x}.scaleAP) ); badMat(badInd) = true; 
    censScans = find(any(badMat(:,zProj), 2))';
    segEdges = [104,90,30,40];
    SegmentCat3D(catInfo(x), 'chan','green', 'zProj',zProj, 'censScans',censScans, 'chunkSize',872, 'overwrite',false, 'minFoot',50, 'edges',segEdges); % 'edges', segParams{x}.edges,  affParams(x).edges , 'censScans',segParams{x}.censScans segParams{x}.zProj
    segParams{x} = GetSegParams( catInfo(x) );
    % Generate ROI
    ROI{x} = MakeROI3D(expt(x), 'overwrite',overwriteROI, 'corrPrct',75);  % , 'corrPrct',90, [ROI{x}, preROI{x}]
    expt(x).Nroi = numel(ROI{x});
    expt(x).roiProj = VisualizeSegmentation(expt(x), ROI{x}, 'overwrite',overwriteROI);   %  
    %ROI{x} = WriteROIproj(expt(x), catInfo(x), ROI{x}, 'edges',segParams{x}.edges, 'overwrite',overwriteROI); % ROI = , 'rSet',32  overwriteROI
    %}

end
