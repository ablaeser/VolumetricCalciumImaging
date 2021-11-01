% Define test data
sbxPath = 'D:\2photon\DL89\171122_DL89\DL89_171122.sbxz'; expt(x).sbx;
refPath = strcat(expt(x).dir, expt(x).name,'_interpProj.tif');
testScans = 1:Nscan; %100:100:expt(x).totScan;   
NtestScans = numel(testScans);
testPlane = 15;
testDir = 'D:\2photon\DL89\171122_DL89\Test\Affine\';

% Load the test data, load metadata and the z-interpolated movie
Nchan = 1;
Nplane = expt(x).Nplane;
Nscan = expt(x).totScan;
Nx = expt(x).Nx; 
Ny = expt(x).Ny; 

rawData = WriteSbxPlaneTif(sbxPath, testPlane);

% Use the z-interpolated projections as the reference
fprintf('\n     Loading %s', refPath); 
refMean = loadtiff( refPath ); 
% GetEdges( refMean(:,:,end), 'minInt',1100, 'show', true );
edges = [67, 174, 47, 70];  % [L, R, T, B]

% Define the reference image
refChan = 1;  %1; %1 = red, 2 = green
refRun = 2;
scanLims = [0, cumsum(expt(x).Nscan)];
refScan = scanLims(refRun)+1:scanLims(refRun+1)-1;
if isempty(refScan), refScan = ceil(Nscan/2)-25:ceil(Nscan/2)+25; end
NrefScan = numel(refScan);
fprintf('\nUsing reference scans %i - %i', refScan(1), refScan(end));
refStack = rawData(:,:,refScan);
refStack = refStack(:,:,find( sum( sum(refStack, 1), 2) )); % exclude blank frames from averaging
refIm = mean( refStack, 3); %uint16();

% Perform affine registrations and extract deformation
fix(clock)
fprintf('\n');
tic
%{
turboParam(1).binxy = 1; turboParam(1).binframe = 1; turboParam(1).hp = false; turboParam(1).pre = false; turboParam(1).sigma = 5; turboParam(1).name = sprintf('%s_nobin_nohp', expt(x).name);
tforms{1} = AffineTurboReg(refIm, 'mov_path',sbxPath, 'optotune_level',testPlane, 'startframe',testScans(1), 'nframes',NtestScans, 'edges',edges, 'dir',testDir, 'name',turboParam(1).name, ...
    'binxy',turboParam(1).binxy, 'highpass',turboParam(1).hp, 'binframes',turboParam(1).binframe, 'pre_register',turboParam(1).pre);
%tforms{1} = pipe.reg.turboreg(refIm, 'mov_path',sbxPath, 'optotune_level',testPlane, 'startframe',testScans(1), 'nframes',NtestScans, 'edges',edges, 'pmt',refChan, ...
%    'binxy',turboParam(1).binxy, 'highpass',turboParam(1).hp, 'binframes',turboParam(1).binframe, 'pre_register',turboParam(1).pre);
deformStruct(1) = ExtractTurboDeform(tforms{1});
toc

turboParam(2).binxy = 2; turboParam(2).binframe = 1; turboParam(2).hp = false; turboParam(2).pre = false; turboParam(2).sigma = 5; turboParam(2).name = sprintf('%s_bin_nohp', expt(x).name);
tforms{2} = AffineTurboReg(refIm, 'mov_path',sbxPath, 'optotune_level',testPlane, 'startframe',testScans(1), 'nframes',NtestScans, 'edges',edges, 'dir',testDir, 'name',turboParam(2).name, ...
    'binxy',turboParam(2).binxy, 'highpass',turboParam(2).hp, 'binframes',turboParam(2).binframe, 'pre_register',turboParam(2).pre);
deformStruct(2) = ExtractTurboDeform(tforms{2});
toc

turboParam(3).binxy = 1; turboParam(3).binframe = 1; turboParam(3).hp = true; turboParam(3).pre = false; turboParam(3).sigma = 5; turboParam(3).name = sprintf('%s_nobin_hp', expt(x).name);
tforms{3} = AffineTurboReg(refIm, 'mov_path',sbxPath, 'optotune_level',testPlane, 'startframe',testScans(1), 'nframes',NtestScans, 'edges',edges, 'dir',testDir, 'name',turboParam(3).name, ...
    'binxy',turboParam(3).binxy, 'highpass',turboParam(3).hp, 'binframes',turboParam(3).binframe, 'pre_register',turboParam(3).pre);
deformStruct(3) = ExtractTurboDeform(tforms{3});
toc
%}
turboParam(4).binxy = 2; turboParam(4).binframe = 1; turboParam(4).hp = true; turboParam(4).pre = false; turboParam(4).sigma = 5; turboParam(4).name = sprintf('%s_z%02d_bin_hp', expt(x).name, testPlane);
tforms{4} = AffineTurboReg(sbxPath, refIm, 'optotune_level',testPlane, 'startframe',testScans(1), 'nframes',NtestScans, 'edges',edges, 'dir',testDir, 'name',turboParam(4).name, ...
    'binxy',turboParam(4).binxy, 'highpass',turboParam(4).hp, 'binframes',turboParam(4).binframe, 'pre_register',turboParam(4).pre);
deformStruct(4) = ExtractTurboDeform(tforms{4});
toc
%{
turboParam(5).binxy = 2; turboParam(5).binframe = 1; turboParam(5).hp = true; turboParam(5).pre = true; turboParam(5).sigma = 5; turboParam(5).name = sprintf('%s_bin_hp_dft', expt(x).name);
tforms{5} = AffineTurboReg(refIm, 'mov_path',sbxPath, 'optotune_level',testPlane, 'startframe',testScans(1), 'nframes',NtestScans, 'edges',edges, 'dir',testDir, 'name',turboParam(5).name, ...
    'binxy',turboParam(5).binxy, 'highpass',turboParam(5).hp, 'binframes',turboParam(5).binframe, 'pre_register',turboParam(5).pre);
deformStruct(5) = ExtractTurboDeform(tforms{5});
%}


%% Compare results - which method works best?
deformVars = fieldnames(deformStruct); Nvars = numel(deformVars);
opt = {[0.03,0.08], [0.06,0.04], [0.1,0.1]};  % {[vert, horz], [bottom, top], [left, right]}
scaleLim = [0.9, 1.1];
colorMat = distinguishable_colors(5);
close all; clearvars sp h;
figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'Color','w', 'PaperOrientation','landscape');
for v = 1:Nvars
    sp(v) = subtightplot(Nvars, 1, v, opt{:});
    if ismember(v, [3,4])
        line( testScans([1,end]), [1,1], 'color','k' ); hold on;
        ylim(scaleLim);
    else
        line( testScans([1,end]), [0,0], 'color','k' ); hold on;
    end
    for r = 4:5
        h(r) = plot( testScans, deformStruct(r).(deformVars{v}), 'Color',colorMat(r,:) ); hold on;
        %pause;
    end
    %if v == 3, legend('Location','northwest'); end % legend(h, deformVars, 'Location','northwest');
    ylabel( deformVars{v}, 'Interpreter','none' );
end
xlabel('Scans');
linkaxes(sp,'x');

%%
% No preprocessing
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','no_pre', 'Nscan',Nscan, ...
        'edges',edges, 'binxy',1, 'binframes',1, 'highpass',0, 'lowpass',0, 'pre_register',false );
    
% edge only
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','edge_only', 'Nscan',Nscan, ...
        'edges',edges+5, 'binxy',1, 'binframes',1, 'highpass',0, 'lowpass',0, 'pre_register',false );    

% DFT only
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','dft_only', 'Nscan',Nscan, ...
        'edges',edges, 'binxy',1, 'binframes',1, 'highpass',0, 'lowpass',0, 'pre_register',true );
    
% HP only
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','hp_only', 'Nscan',Nscan, ...
        'edges',edges, 'binxy',1, 'binframes',1, 'highpass',5, 'lowpass',0, 'pre_register',false );
% LP only
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','lp_only', 'Nscan',Nscan, ...
    'edges',edges, 'binxy',1, 'binframes',1, 'highpass',0, 'lowpass',1, 'pre_register',false );

% BinXY only
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','binxy_only', 'Nscan',Nscan, ...
        'edges',edges, 'binxy',2, 'binframes',1, 'highpass',0, 'lowpass',0, 'pre_register',false );
    
% bin t only
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','bint_only', 'Nscan',Nscan, ...
        'edges',edges, 'binxy',1, 'binframes',2, 'highpass',0, 'lowpass',0, 'pre_register',false );
    
% median only
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','med_only', 'Nscan',Nscan, ...
    'edges',edges, 'binxy',1, 'binframes',1, 'highpass',0, 'lowpass',0, 'pre_register',false, 'medFilter',[1,1,1] );

%% Edge tests
% edge only
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','edge_only', 'Nscan',Nscan, ...
        'edges',edges+5, 'binxy',2, 'binframes',1, 'highpass',5, 'lowpass',1, 'pre_register',false );    
    
% edge_hp
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','edge_hp', 'Nscan',Nscan, ...
        'edges',edges+5, 'binxy',1, 'binframes',1, 'highpass',5, 'lowpass',0, 'pre_register',false );
    
% edge_lp
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','edge_lp', 'Nscan',Nscan, ...
        'edges',edges+5, 'binxy',1, 'binframes',1, 'highpass',0, 'lowpass',1, 'pre_register',false );
        
% edge_bp
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','edge_bp', 'Nscan',Nscan, ...
        'edges',edges+5, 'binxy',1, 'binframes',1, 'highpass',5, 'lowpass',1, 'pre_register',false );

% edge_pre
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','edge_pre', 'Nscan',Nscan, ...
        'edges',edges+5, 'binxy',1, 'binframes',1, 'highpass',0, 'lowpass',0, 'pre_register',true );
    
% edge_pre_bp
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','edge_pre', 'Nscan',Nscan, ...
        'edges',edges+5, 'binxy',1, 'binframes',1, 'highpass',5, 'lowpass',1, 'pre_register',true );
    
% edge_med
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','edge_med', 'Nscan',Nscan, ...
        'edges',edges+5, 'binxy',1, 'binframes',1, 'highpass',0, 'lowpass',0, 'pre_register',false, 'medFilter',[1,1,1] );

% edge_binxy_bp_dft
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','edge_binxy_bp_dft', 'Nscan',Nscan, ...
    'edges',edges+5, 'binxy',2, 'binframes',1, 'highpass',5, 'lowpass',1, 'pre_register',true );

% edge_binxy_hp_med_dft
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','edge_binxy_hp_med_dft', 'Nscan',Nscan, ...
    'edges',edges+5, 'binxy',2, 'binframes',1, 'highpass',5, 'lowpass',0, 'medFilter',[1,1,1], 'pre_register',true );

% edge_binxy_med_dft
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','edge_binxy_med_dft', 'Nscan',Nscan, ...
    'edges',edges+5, 'binxy',2, 'binframes',1, 'highpass',0, 'lowpass',0, 'medFilter',[1,1,1], 'pre_register',true );

% edge_med_dft !!
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','edge_med_dft', 'Nscan',Nscan, ...
    'edges',edges+5, 'binxy',1, 'binframes',1, 'highpass',0, 'lowpass',0, 'medFilter',[3,3,3], 'pre_register',true );

% edge_med_binxy_dft 
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','edge_med_binxy_dft', 'Nscan',Nscan, ...
    'edges',edges+5, 'binxy',2, 'binframes',1, 'highpass',0, 'lowpass',0, 'medFilter',[3,3,3], 'pre_register',true );

% binxy_med_dft
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','binxy_med_dft', 'Nscan',Nscan, ...
    'edges',edges, 'binxy',2, 'binframes',1, 'highpass',0, 'lowpass',0, 'medFilter',[1,1,1], 'pre_register',true );

% edge_binxy_bp
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','edge_binxy_bp', 'Nscan',Nscan, ...
    'edges',edges+5, 'binxy',2, 'binframes',1, 'highpass',5, 'lowpass',1, 'pre_register',false );  

% edge10_binxy_pre_bp
AffineTurboReg(sbxPath, sbxInfo, refVol(:,:,z), 'pmt',refChan, 'optotune_level',z, 'dir',tifDir, 'name','edge10_binxy_bp_dft', 'Nscan',Nscan, ...
    'edges',edges+10, 'binxy',2, 'binframes',1, 'highpass',5, 'lowpass',1, 'pre_register',true );

%% Test effects of different isolated affine transforms
inputImage = loadtiff('D:\2photon\DL67\170601\DL67_170601_affineProj.tif');
%{
deformCatStruct.transAP(s,z) = affineData.tforms_all{z,s}.T(3,1); % Trans AP (right/left +/-) trans_x
deformCatStruct.transML(s,z) = affineData.tforms_all{z,s}.T(3,2); % Trans ML (down/up +/-)  trans_y
deformCatStruct.scaleAP(s,z) = affineData.tforms_all{z,s}.T(1,1); % Scale AP (inflate/deflate >/< 1)  scale_x
deformCatStruct.scaleML(s,z) = affineData.tforms_all{z,s}.T(2,2); % Scale ML (inflate/deflate >/< 1) scale_y
deformCatStruct.shearAP(s,z) = affineData.tforms_all{z,s}.T(1,2); % Shear AP (tilt left/right +/-)  shear_x
deformCatStruct.shearML(s,z) = affineData.tforms_all{z,s}.T(2,1); % Shear ML (tilt down/right +/-) shear_y
%}
transAP = 0;
transML = 0;
scaleAP = 1;
scaleML = 1;
shearAP = 0.02;
shearML = 0;


tform = affine2d();
tform.T = [scaleAP, shearAP, 0; shearML, scaleML, 0; transAP, transML, 1]; %[1, 0, 0; 0, 1, 0; 0, 0, 1];
outputImage = imwarp(inputImage, tform, 'OutputView',imref2d(size(inputImage)) );

close all; clearvars sp;
figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
sp(1) = subplot(1,2,1); 
axis image;
imshow( inputImage, [] )
sp(2) = subplot(1,2,2);
imshow( outputImage, [] )
axis image;
linkaxes(sp,'xy');
impixelinfo;
