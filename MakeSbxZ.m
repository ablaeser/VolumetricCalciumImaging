function MakeSbxZ(sbxPath, sbxInfo, shiftPath, varargin)

p = inputParser;
addOptional(p,'edges',[0,0,0,0]);
parse(p,varargin{:});
p = p.Results;

load(shiftPath,'-mat');
RS = RS1 + RS2 + RS3 + RS_chunk;
CS = CS1 + CS2 + CS3 + CS_chunk;
ZS = ZS1 + ZS_chunk;
%{
figure;
imagesc( ZS );
set(gca, 'Xtick',[20,50,60,90,100,105,150,190,195]) % changed Test scans
impixelinfo;
%}

tic
fprintf('\nLoading %s', sbxPath); tic;
sbx_data = readSBX(sbxPath, sbxInfo, 1, sbxInfo.Nscan, -1, []); toc

fprintf('\nApplying z interpolation shifts');  tic;
if sbxInfo.nchan == 1
    sbx_data = ApplyZShiftInterpolateFBS(sbx_data, ZS, CS, RS); toc % apply shifts, this step is slow
    sbx_data = reshape(sbx_data, [size(sbx_data,[1,2]), prod(size(sbx_data,[3,4]))] ); 
else
    sbx_data = permute(sbx_data, [2,3,4,5,1]); % ApplyZShiftInterpolateFBS expects spatial dimensions first
    % Interpolate each color, serial
    for c = 1:2
        sbx_data(:,:,:,:,c) = ApplyZShiftInterpolateFBS(sbx_data(:,:,:,:,c), ZS, CS, RS);
    end
    % Interpolate each color in parallel - seems to be slower than serial
    %{
    sbx_data = num2cell(sbx_data, [1,2,3,4]);
    sbx_data = squeeze(sbx_data);
    tic
    parfor c = 1:2
        sbx_data{c} = ApplyZShiftInterpolateFBS(sbx_data{c}, ZS, CS, RS);
    end
    sbx_data = cat(5, sbx_data{:}); toc
    sbx_data = permute(sbx_data, [5,1,2,3,4]);  toc % regwriter expects color dim first
    %}
    sbx_data = permute(sbx_data, [5,1,2,3,4]); toc % regwriter expects color dim first
    sbx_data = reshape(sbx_data, [size(sbx_data,[1,2,3]), prod(size(sbx_data,[4,5]))] ); 
end
toc

fprintf('\nWriting .sbxz'); tic
rw = pipe.io.RegWriter(sbxPath, sbxInfo, '.sbxz', true); 
rw.write(sbx_data);
rw.delete;
toc

%{
%open regwriter object for writing z-interpolated figure
rw = pipe.io.RegWriter(sbxPath, sbxInfo, '.sbxz', true);
w = waitbar(0,'writing .sbxz');

%work one volume at a time
for s = 1:p.Nt
    %read the volume
    raw_vol = readSBX(sbxPath, sbxInfo, s, 1, -1, []); % pipe.imread(path,Nz*(i-1)+1, Nz,-1,[]);
    %crop it based on edges
    raw_vol = raw_vol(:,p.edges(3)+1:end-p.edges(4),p.edges(1)+1:end-p.edges(2),:);
    %if using optotune correction, perform correction
    if sbxInfo.nchan == 1
        warpvol(1,:,:,:) = raw_vol;
    else
        warpvol = raw_vol;
    end
    %get the DFT shifts for this volume
    CS_vol = CS(:,s);
    RS_vol = RS(:,s);
    ZS_vol = ZS(:,s);

    reg_vol = zeros(size(warpvol));
    for c = 1:sbxInfo.nchan
        reg_vol(c,:,:,:) = ApplyZShiftInterpolateFBS(squeeze(warpvol(c,:,:,:)),ZS_vol,CS_vol,RS_vol);
    end

    %crop the registered image based on edge-sizes
    reg_vol_masked = zeros(size(reg_vol,1),size(reg_vol,2)+p.edges(3)+p.edges(4),size(reg_vol,3)+p.edges(1)+p.edges(2),size(reg_vol,4));
    reg_vol_masked(:,p.edges(3)+1:p.edges(3)+size(reg_vol,2),p.edges(1)+1:p.edges(1)+size(reg_vol,3),:) = reg_vol;
    reg_vol_masked = uint16(reg_vol_masked);

    %write volume to .sbxz file
    rw.write(squeeze(reg_vol_masked));
    waitbar(s/p.Nt);
end
rw.delete;
delete(w);
%}
end