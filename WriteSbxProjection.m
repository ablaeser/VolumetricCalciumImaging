function [projStack, tifStack] = WriteSbxProjection(sbxPath, sbxInfo, projPath, varargin)
% Extract data from an SBX file, and then perform a projection (mean or max) across frames
IP = inputParser;
addRequired( IP, 'sbxPath', @ischar )
addRequired( IP, 'sbxInfo', @isstruct)
addOptional( IP, 'projPath', '', @ischar)
addParameter( IP, 'chan', 'both', @ischar ) % 'green', 'red', 'both'. for scanbox, 1 = green, 2 = red. -1 = both
%addParameter( IP, 'chan', -1, @isnumeric ) % for scanbox, 1 = green, 2 = red. -1 = both
addParameter( IP, 'z', [], @isnumeric ) 
addParameter( IP, 'projType', 'mean', @ischar)
addParameter( IP, 'firstScan', 1, @isnumeric )
addParameter( IP, 'Nscan', -1, @isnumeric )
addParameter( IP, 'edges', [0,0,0,0], @isnumeric ) % [left, right, top, bottom]
addParameter( IP, 'scale', 1, @isnumeric ) %
addParameter( IP, 'binT', 1, @isnumeric ) 
addParameter( IP, 'dir', '', @ischar ) 
addParameter( IP, 'name', '', @ischar )
addParameter( IP, 'type', '', @ischar ) 
addParameter( IP, 'monochrome', false, @islogical ) % write a monochrome movie for each plane?
addParameter( IP, 'RGB', false, @islogical ) % write an RGB monochrome movie for each plane?
addParameter( IP, 'rescale', false, @islogical )
addParameter( IP, 'verbose', true, @islogical );
addParameter( IP, 'overwrite', false, @islogical );
parse( IP, sbxPath, sbxInfo, projPath, varargin{:} ); 
writeChan = IP.Results.chan;  
projType = IP.Results.projType;
firstScan = IP.Results.firstScan;
edges = IP.Results.edges;
scaleFactor = IP.Results.scale;
binT = IP.Results.binT;
verbose = IP.Results.verbose;
Nscan = IP.Results.Nscan;
zSet = IP.Results.z;
saveDir = IP.Results.dir;  mkdir(saveDir);
saveName = IP.Results.name;
sbxType = IP.Results.type;
monochrome = IP.Results.monochrome;
RGB = IP.Results.RGB;
rescaleIntToggle = IP.Results.rescale;
overwrite = IP.Results.overwrite;
if isempty(zSet), zSet = 1:sbxInfo.Nplane; end
%zSet = flip(zSet);
maxZ = @(x)(max(x,[],3));
if ~isempty(sbxType), sbxType = strcat('_',sbxType); end
fileExists = exist(projPath,'file');
outputToggle = nargout > 0;
tic;
% Determine which channel(s) to use

[usePMT, usePMTname] = DeterminePMT(writeChan, sbxInfo); % usePMT
if ~(~outputToggle && ~overwrite && fileExists)
    % Get each plane, crop, resize, and mean project (can't necessarily get the full data at once due to memory constraints)
    %projStack = zeros(sbxInfo.width, sbxInfo.height, sbxInfo.Nplane, 2); % color order: red, green
    for Z = flip(1:numel(zSet)) % flip(zSet) %flip(1:sbxInfo.otlevels)
        if verbose, fprintf('\nZ = %i', Z); end
        %if ~isempty(saveName), tempName = saveName; else, tempName = ''; end  % sprintf('%s_plane%01i',saveName, Z);
        [~, tempChan] = WriteSbxPlaneTif(sbxPath, sbxInfo, zSet(Z), 'verbose',verbose, 'dir',saveDir, 'name',saveName, 'overwrite',overwrite, ...
            'edges',edges, 'scale',scaleFactor, 'firstScan',firstScan, 'Nscan',Nscan, 'chan',usePMTname, 'zeros',true, 'binT',binT, 'RGB',RGB, 'monochrome',monochrome ); % 
        if strcmpi(projType,'mean')
            tempProj = cellfun(@mean, tempChan, {3,3}, 'UniformOutput',false);
        else
            tempProj = cellfun(maxZ, tempChan, 'UniformOutput',false);
        end
        for pmt = usePMT %find(~cellfun(@isempty, tempProj)) %intersect(writeChanInd, )
            projStack(:,:,zSet(Z),pmt) = tempProj{pmt};
        end
    end
    projStack = projStack(:,:,:,usePMT);

    if verbose, toc; end
    % Save tif (optional)
    if ~isempty( projPath ) || overwrite
        if verbose, fprintf('\nWriting %s', projPath); end
        if numel(usePMT) > 1 %strcmpi(writeChan, 'both')
            chanName = {'red','green'}; % {'green','red'}; %PMTname = 
            chanInd = [1,2]; %[2,1]; 
            [projDir, ~] = fileparts(projPath);
            tifStack = zeros(size(projStack,1), size(projStack,2), size(projStack,3), 3);
            for chan = usePMT
                saveastiff(uint16(projStack(:,:,:,chan)), sprintf('%s\\%s_%s.tif', projDir, sbxInfo.fileName, chanName{chan}));
                if ~rescaleIntToggle
                    tifStack(:,:,:,chanInd(chan)) = projStack(:,:,:,chan)/256;
                else
                    tempStack = projStack(:,:,:,chan);
                    chanLower = prctile(tempStack(:), 1);
                    chanUpper = max(tempStack(:)); %prctile(stackChan{chan}(:), 1);
                    if verbose, fprintf('\nRescaling %s channel: [%i, %i] -> [0, 255]', chanName{chan}, chanLower, chanUpper); end
                    tifStack(:,:,:,chanInd(chan)) = rescale(tempStack, 0, 2^8-1, 'inputMin',chanLower); % min(stackChan{chan}(:))
                end
            end
            tifStack = uint8(tifStack);
        else
            tifStack = uint16(projStack); % (:,:,:,writeChanInd)
        end
        WriteTiff(tifStack, projPath); %pipe.io.writeTiff(tifStack, projPath);
        if verbose, fprintf('... done!\n'); end
    end
elseif fileExists && outputToggle
    tifStack = loadtiff(projPath);
    tifStack = permute(tifStack, [1,2,4,3]);
    projStack = flip(tifStack(:,:,:,[1,2]),4);
else
    if verbose, fprintf('\n%s already exists and no output was requested', projPath); end
end
end