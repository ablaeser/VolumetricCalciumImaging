function [roiRGB, axonRGB] = VisualizeSegmentation(expt, ROI, varargin) % sbxPath, 

IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'ROI', @isstruct )
addOptional( IP, 'axon', [], @isstruct )
addParameter( IP, 'overwrite', false, @islogical)
parse( IP, expt, ROI, varargin{:} );  
axon = IP.Results.axon;
overwrite = IP.Results.overwrite;
RGBOpt = struct('overwrite',true, 'message',false, 'append',false, 'big',true, 'color',true );
roiDir = strcat(expt.dir, 'ROI\'); mkdir(roiDir);
Nroi = numel(ROI);
Naxon = numel(axon);

% ROI
roiCmap = distinguishable_colors(Nroi);
roiLabelMat = zeros(expt.Nrow, expt.Ncol, expt.Nplane); % size(maxProj)
for r = 1:Nroi,  roiLabelMat(ROI(r).ind) = r;  end
if expt.Nplane > 1
    roiRGB = label2rgb3d(roiLabelMat,roiCmap,'w');
    roiRGB = permute(roiRGB, [1,2,4,3]);
else
    roiRGB = label2rgb(roiLabelMat,roiCmap,'w');
end
roiProjPath = strcat(expt.dir, expt.name,'_roiProj.tif');
if ~exist(roiProjPath, 'file') || overwrite % && Nroi > 0
    fprintf('\nWriting %s\n', roiProjPath);
    saveastiff( uint16(roiRGB), roiProjPath, RGBOpt );
end

% Axons
if Naxon > 0
    axonCmap = distinguishable_colors(Naxon);      
    axonLabelMat = max(cat(4, axon.labelVol ),[], 4);
    if expt.Nplane > 1
        axonRGB = label2rgb3d(axonLabelMat,axonCmap,'w');
        axonRGB = permute(axonRGB, [1,2,4,3]);
    else
        axonRGB = label2rgb(axonLabelMat,axonCmap,'w');
    end
    axonProjPath = strcat(expt.dir, expt.name,'_axonProj.tif');
    if ~exist(axonProjPath, 'file') || overwrite
        fprintf('\nWriting %s\n', axonProjPath);
        saveastiff( uint16(axonRGB), axonProjPath, RGBOpt );
    end
else
    axonRGB = zeros(size(roiRGB));
end
end