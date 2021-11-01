function out = WriteSbxZproj(sbxPath, sbxInfo, projPath, firstScan, Nscan, pmt, zProj, varargin)
p = inputParser;
addOptional(p, 'proj_type', 'mean');
addOptional(p, 'write_tiff', true); %check whether to save the tproj as a tiff, or just return it
parse(p, varargin{:});
p = p.Results;
if pmt == -1 && sbxInfo.nchan == 2
    Nchan = 2;
    projDim = 4;
else
    Nchan = 1;
    projDim = 3;
end
%load each volume scan (parallelized) and z-project
tic
A = cell(firstScan+Nscan-1,1);
h = parfor_progressbar(Nscan,'z-projecting...');
parfor scan = firstScan:firstScan+Nscan-1
    vol = readSBX(sbxPath, sbxInfo, scan, 1, pmt, []);
    if Nchan == 1
        vol = vol(:,:,zProj);
    else
        vol = vol(:,:,:,zProj);
    end
    vol(vol==0) = NaN;
    if strcmpi(p.proj_type,'mean')
        slice = squeeze( mean(vol,projDim,'omitnan') );
    elseif strcmpi(p.proj_type,'max')
        slice = squeeze( max(vol,[],projDim,'omitnan') );
    elseif strcmpi(p.proj_type,'median')
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
out = cat(projDim, A{:});
clearvars A;

%write the file
if p.write_tiff
    if ndims(out) == 4
        out = permute(out, [2,3,1,4]);
    end
    %out = reshape(out,size(out,1),size(out,2),[]);
    fprintf('\nWriting %s', projPath); 
    pipe.io.writeTiff(uint16(out), projPath);
end
end