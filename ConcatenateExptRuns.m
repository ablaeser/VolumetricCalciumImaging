function interRunShift = ConcatenateExptRuns(expt, runInfo, catInfo, varargin )
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'runInfo', @isstruct )
addRequired( IP, 'catInfo', @isstruct )
%addParameter( IP, 'overwrite', false, @islogical )
addParameter( IP, 'minInt', 1000, @isnumeric )
addParameter( IP, 'maxEdge', [100, 100, 120, 120], @isnumeric ) % [60,60,40,40]
addParameter( IP, 'setEdge', [], @isnumeric )
addParameter( IP, 'refRun', 1, @isnumeric )
addParameter( IP, 'refChan', 'green', @ischar )
addParameter( IP, 'sbx', 'sbxz', @ischar )
addParameter( IP, 'proj', 'interp', @ischar )
addParameter( IP, 'ext', 'sbxcat', @ischar )
addParameter( IP, 'overwrite', false, @islogical ) % for scanbox, 1 = green, 2 = red. -1 = both
parse( IP, expt, runInfo, catInfo, varargin{:} ); % mouse, exptDate,
minIntOrig = IP.Results.minInt;
if expt.Nchan == 2 && expt.Nplane > 1
    minInt = minIntOrig/2^8;
else
    minInt = minIntOrig;
end
edgeMax = IP.Results.maxEdge;
edgeSet = IP.Results.setEdge;
sbxType = IP.Results.sbx;
projType = IP.Results.proj;
refRun = IP.Results.refRun;
refChan = IP.Results.refChan;
refChanInd = find(strcmpi(refChan, {'red','green'} ));
if isempty(refChanInd),  error('Invalid reference channel'); end
catExt = IP.Results.ext; %'.sbxcat'; % '.sbxcat'
overwrite = IP.Results.overwrite;

catDir = strcat(expt.dir,'Concat\');  mkdir(catDir)
catName = expt.name; % sprintf('%s_FOV%i', , expt.fov);
catPathRoot = sprintf('%s%s', expt.dir, catName);
%catInfoPath = strcat(catPathRoot, '.mat' ); %[expt.dir, expt.name, expt.fov, '.mat'];
catSbxPath = strcat(catPathRoot, '.', catExt);
catProjPath = sprintf('%s%s_catProj.tif', expt.dir, catName );
interRunShift = zeros(expt.Nruns, 3);
if expt.Nruns > 1
    if overwrite || ~exist(catSbxPath,'file')
        Nchan = runInfo(1).nchan; Nscan = [runInfo.Nscan];  %totScan = sum(Nscan);
        sbxPath = cell(1,expt.Nruns);  runProjPath = cell(1,expt.Nruns);  runProj = cell(expt.Nruns, 1);
        cropProj = cell(expt.Nruns, 1); shiftProj = cell(1,expt.Nruns); nullProj = cell(1,expt.Nruns); %hpProj = cell(expt.Nruns, 1);
        if isempty(edgeSet), projEdges = zeros( expt.Nruns, 4 ); else, projEdges = repmat(edgeSet, expt.Nruns, 1 ); end
        if expt.Nplane > 1
            Nplane = runInfo(1).Nplane;   % Nrow = runInfo(1).sz(1); Ncol = runInfo(1).sz(2);
            for r = 1:expt.Nruns % expt.runs
                sbxPath{r} = sprintf('%s%s.%s', runInfo(r).dir, runInfo(r).fileName, sbxType ); % sbx_affine  sbxz
                runProjPath{r} = sprintf('%s%s_%sProj.tif', runInfo(r).dir, runInfo(r).fileName, projType ); % affine
                runProj{r} = loadtiff( runProjPath{r} );
                if ndims(runProj{r}) == 4, runProj{r} = squeeze(runProj{r}(:,:,refChanInd,:)); end
                if isempty(edgeSet)
                    projEdges(r,:) = GetEdges( runProj{r}(:,:,end), 'minInt',minInt, 'show',true );
                else
                    if ndims(runProj{r}) == 4
                        ShowEdges(projEdges(r,:), runProj{r}(:,:,refChanInd,end));
                    else
                        ShowEdges(projEdges(r,:), runProj{r}(:,:,end));
                    end
                    %pause;
                end
                %saveastiff( runProj{r}, sprintf('%sTest\\%s_Proj_run%i.tif', expt.dir, expt.name, r ) );
            end
            close all;
            % keep all edges within edgeMax
            if isempty(edgeSet)
                edgeSet = max(projEdges, [], 1); %[90,100,60,100]; % %edgeSet(2) =  edgeSet(2) + (683-539)
            end
            aboveEdge = edgeMax - edgeSet < 0;
            edgeSet(aboveEdge) = edgeMax(aboveEdge);


            %filtSigma = 5;
            for r = 1:expt.Nruns % expt.runs
                cropProj{r} = runProj{r}(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),:);
                %hpProj{r} = cropProj{r} - imgaussfilt(cropProj{r}, filtSigma);
                saveastiff( cropProj{r}, sprintf('%s%s_Proj_run%i.tif', catDir, expt.name, r ) );
                %saveastiff( hpProj{r}, sprintf('%s%s_HP_run%i.tif', catDir, expt.name, r ) );
            end

            % Calculate shifts between runs
            zUse = 1:Nplane;
            %refRun = 2;  %10:Nplane;
            interRunDFTshift = zeros(expt.Nruns,3);  interRunZshift = zeros(numel(zUse), expt.Nruns);
            for r = 1:expt.Nruns
                if r ~= refRun
                    if Nplane > 1
                        interRunDFTshift(r,:) = dftregistration3D(fftn(cropProj{refRun}(:,:,zUse)),  fftn(cropProj{r}(:,:,zUse)), 4);
                        interRunZshift(:,r) = InterpolateZshift(cropProj{refRun}(:,:,zUse),  cropProj{r}(:,:,zUse), 2);
                    else
                        zUse = 1;
                        tempOut = dftregistration(fft2(cropProj{refRun}),  fft2(cropProj{r})); % [error,diffphase,net_row_shift,net_col_shift]
                        interRunDFTshift(r,1:2) = tempOut(3:4);
                    end
                end
            end
            %interRunShift = zeros(expt.Nruns,3);
            interRunShift(:,[1,2,3]) = interRunDFTshift(:,[2,1,3]); % imtranslate expects [x,y,z], not [y,x,z]
            interRunShift = round( interRunShift );

            % Apply shifts to each run's projection
            shiftProj = cell(expt.Nruns,expt.Nchan);
            for r = 1:expt.Nruns
                if Nplane > 1
                    shiftProj{r} = imtranslate( runProj{r}, interRunShift(r,:) );
                else
                    shiftProj{r} = imtranslate( runProj{r}, interRunShift(r,1:2) );
                end
                saveastiff( [cropProj{refRun}, shiftProj{r}(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),:)], sprintf('%s%s_shiftProj_run%i.tif', catDir, expt.name, r ) );
            end
            shiftProjCat = cat( 4, shiftProj{:} ); % projCat = cat( 4, nullProj{:} );
            for z = 1:Nplane
                %saveastiff( squeeze( projCat(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),z,:) ), sprintf('%sTest\\%s_ind_z%i.tif', expt.dir, expt.name, z ) )
                saveastiff( squeeze( shiftProjCat(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),z,:) ), sprintf('%s%s_cat_z%i.tif', catDir, expt.name, z ) )
            end

            % Make metadata file
            %catInfo = ConcatenateRunInfo(expt, runInfo, 'suffix','sbxcat'); %SpoofSBXinfo3D(Nx, Ny, Nplane, totScan, Nchan);
            %save( catInfoPath, 'info' );
        else
            for r = expt.runs
                sbxPath{r} = sprintf('%s%s.sbx', runInfo(r).dir, runInfo(r).fileName ); % sbx_affine  sbxz
                runProjPath{r} = sprintf('%s%s_%s.tif', runInfo(r).dir, runInfo(r).fileName, refChan ); % affine
                runProj{r} = loadtiff( runProjPath{r} );
                if isempty(edgeSet)
                    projEdges(r,:) = GetEdges( runProj{r}(:,:,end), 'minInt',minInt, 'show',true );
                else
                    ShowEdges(projEdges(r,:), runProj{r});
                end
                close all;
            end

            % keep all edges within edgeMax
            if isempty(edgeSet)
                edgeSet = max(projEdges, [], 1); %[90,100,60,100]; % %edgeSet(2) =  edgeSet(2) + (683-539)
            end
            aboveEdge = edgeMax - edgeSet < 0;
            edgeSet(aboveEdge) = edgeMax(aboveEdge);
            % Crop the reference projections
            for r = expt.runs
                cropProj{r} = runProj{r}(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),:);
                saveastiff( cropProj{r}, sprintf('%s%s_Proj_run%i.tif', catDir, expt.name, r ) );
            end

            close all;

            % Calculate shifts between runs
            interRunDFTshift = zeros(expt.Nruns,2);
            for r = 1:expt.Nruns
                if r ~= refRun
                    tempOut = dftregistration(fft2(cropProj{refRun}),  fft2(cropProj{r})); % [error,diffphase,net_row_shift,net_col_shift]  hpProj
                    interRunDFTshift(r,1:2) = tempOut(3:4);
                end
            end
            interRunShift = round( flip(interRunDFTshift, 2) ); % flip shifts to be compatible with imtranslate

            % Apply shifts to each run's projection
            for r = 1:expt.Nruns
                nullProj{r} = imtranslate( runProj{r}, interRunShift(r,1:2) );
                saveastiff( nullProj{r}(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),:), sprintf('%s%s_nullProj_run%i.tif', catDir, expt.name, r ) );
                shiftProj{r} = nullProj{r};
                saveastiff( [cropProj{refRun}, shiftProj{r}(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),:)], sprintf('%s%s_shiftProj_run%i.tif', catDir, expt.name, r ) );
            end
            projCat = cat( 4, nullProj{:} ); shiftProjCat = cat( 4, shiftProj{:} );
            for z = 1:expt.Nplane
                saveastiff( squeeze( projCat(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),z,:) ), sprintf('%s%s_ind_z%i.tif', catDir, expt.name, z ) )
                saveastiff( squeeze( shiftProjCat(edgeSet(3):end-edgeSet(4),edgeSet(1):end-edgeSet(2),z,:) ), sprintf('%s%s_cat_z%i.tif', catDir, expt.name, z ) )
            end
        end

        % Write the sbxcat file
        fprintf('\n     Writing %s\n', catSbxPath); tic
        rw = pipe.io.RegWriter(catSbxPath, catInfo, catExt, true);
        w = waitbar(0, sprintf('writing %s',catExt));
        for r = 1:expt.Nruns %expt.runs
            % Load run stack
            fprintf('\n   Loading %s... ', sbxPath{r}); tic
            
            if Nplane > 1
                runStack = readSBX(sbxPath{r}, runInfo(r), 1, Nscan(r), -1, []); % [c,x,y,z,t]
                if Nchan == 2
                    runStack = permute(runStack, [2,3,4,5,1]); % [x,y,z,t,c]
                    % Apply shifts to each scan
                    if ~all(interRunShift(r,:) == 0)
                        for s = 1:Nscan(r)
                            for c = 1:Nchan
                                runStack(:,:,:,s,c) = imtranslate( runStack(:,:,:,s,c), interRunShift(r,:)  );
                            end
                        end
                    end
                    runStack = permute(runStack, [5,1,2,3,4]); % [c,x,y,z,t]
                    runStack = reshape(runStack, [size(runStack,[1,2,3]), prod(size(runStack,[4,5]))]);
                    rw.write( runStack ); %
                elseif Nchan == 1
                    % Apply shifts to each scan
                    if ~all(interRunShift(r,:) == 0)
                        for s = 1:Nscan(r)
                            runStack(:,:,:,s) = imtranslate( runStack(:,:,:,s), interRunShift(r,:)  );
                        end
                    end
                    % Write the data to sbxcat
                    runStack = reshape(runStack, [size(runStack,[1,2]), prod(size(runStack,[3,4]))] );
                    rw.write( runStack ); % rw.write(squeeze(uint16(tempScan)));
                end
            else
                if Nchan == 2
                    runStack = readSBX(sbxPath{r}, runInfo(r), 1, Nscan(r), -1, []); % [c,x,y,t]
                    runStack = permute(runStack, [2,3,4,1]); % [x,y,t,c]
                    % Apply shifts to each scan
                    if ~all(interRunShift(r,:) == 0)
                        for s = 1:Nscan(r)
                            for c = 1:Nchan
                                runStack(:,:,s,c) = imtranslate( runStack(:,:,s,c), interRunShift(r,:)  );
                            end
                        end
                    end
                    runStack = permute(runStack, [4,1,2,3]); % [c,x,y,t]
                    rw.write( runStack ); %
                elseif Nchan == 1
                    runStack = readSBX(sbxPath{r}, runInfo(r), 1, Nscan(r), -1, []); % [x,y,t]
                    % Apply shifts to each scan
                    if ~all(interRunShift(r,:) == 0)
                        for s = 1:Nscan(r)
                            runStack(:,:,s) = imtranslate( runStack(:,:,s), interRunShift(r,:)  );
                        end
                    end
                    rw.write( runStack ); % Write the data to sbxcat
                end
            end
            waitbar( r/expt.Nruns, w );
            toc
        end
        rw.delete;
        delete(w);
    else
        fprintf('\n%s already exists!', catSbxPath);
    end
    %{
    for z = [2,4,6,8,10,14]
        WriteSbxPlaneTif(catSbxPath, catInfo, z, 'dir',catDir, 'name',[catName,'_cat'], 'overwrite',true  );
    end
    %}
    WriteSbxProjection(catSbxPath, catInfo, catProjPath, 'dir',catDir, 'name',[catName,'_cat'], 'overwrite',overwrite); % , 'firstScan',100, 'Nscan',500
elseif ~exist(catSbxPath, 'file') || overwrite
    sbxSourcePath = strcat(runInfo.path,'z');
    fprintf('\nSingle run experiment: copying %s to %s', sbxSourcePath, catSbxPath);
    copyfile(sbxSourcePath, catSbxPath);
    projSourcePath = strcat(runInfo.dir, runInfo.fileName, '_interpProj.tif');
    copyfile(projSourcePath, catProjPath);
end
end