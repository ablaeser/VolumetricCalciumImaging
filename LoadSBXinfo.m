function info = LoadSBXinfo(inputPath) % input path can be .sbx or .mat

% Make sure we're opening the info .mat file
[parent, fname, ext] = fileparts(inputPath);

% Find the name of the info file and load it
if ~isempty(parent), parent = [parent '\']; end
basePath = [parent fname];
infoPath = [basePath '.mat'];
load(infoPath); % creates the info structure 

% Determine path of the corresponding sbx file
if ~contains(ext, 'sbx')
    sbxPath = [basePath '.sbx'];
else
    sbxPath = [basePath, ext];
end
%info.fid = fopen(sbxPath);

% If bidirectional scanning, double the records per buffer
if info.scanmode == 0,  info.recordsPerBuffer = info.recordsPerBuffer*2;  end
% How many channels were acquired?
switch info.channels
    case 1
        info.nchan = 2;      % both PMT0 & 1
    case 2
        info.nchan = 1;      % PMT 0
    case 3
        info.nchan = 1;      % PMT 1
    case -1  % added for David Tingley 10/29/21 - may not apply in general
        info.nchan = info.chan.nchan; %2;
end
if info.nchan == 1
    factor = 2;
elseif info.nchan == 2
    factor = 1;
else
    error('Expect nchan = 1 or 2 only');
end


info.nsamples = (info.sz(2)*info.recordsPerBuffer*2*info.nchan);   % bytes per record

info.Nrow = info.sz(1); 
info.Ncol = info.sz(2);

d = dir(sbxPath);
if isfield(info, 'scanbox_version') && info.scanbox_version >= 2
    info.max_idx = d.bytes/info.recordsPerBuffer/info.sz(2)*factor/4 - 1;
else
    info.max_idx =  d.bytes/info.bytesPerBuffer*factor - 1;
end

% Appended useful information
info.nframes = info.max_idx + 1;
info.optotune_used = false;
info.otlevels = 1;
if isfield(info, 'volscan') && info.volscan > 0, info.optotune_used = true; end
if ~isfield(info, 'volscan') && ~isempty(info.otwave), info.optotune_used = true; end
if info.optotune_used
    if isfield(info, 'otwave')
        info.otlevels = length(info.otwave);
    elseif isfield(info, 'otparam')
        info.otlevels = 30;
    else
        error('Unable to determine number of optotune planes used');
    end
end

info.framerate = 15.49; %
if info.scanmode == 0,  info.framerate = 30.98;  end

info.height = info.sz(2);
info.width = info.recordsPerBuffer;

end