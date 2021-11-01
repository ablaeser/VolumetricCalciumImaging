function [chunkLims, Nchunk, chunkLength] = MakeChunkLims(firstScan, lastScan, totScan, varargin)
IP = inputParser;
addRequired( IP, 'firstScan', @isnumeric )
addRequired( IP, 'lastScan', @isnumeric )
addOptional( IP, 'totScan', lastScan, @isnumeric )
addParameter( IP, 'N',NaN, @isnumeric )
addParameter( IP, 'size',NaN, @isnumeric )
parse( IP, firstScan, lastScan, totScan, varargin{:} );
totScan = IP.Results.totScan;
%if isinf(totScan),totScan = lastScan;
Nchunk = IP.Results.N;
chunkSize = IP.Results.size;
%clearvars chunkLims
if isnan(Nchunk)
    chunkLims(:,1) = (firstScan:chunkSize:lastScan)'; 
    chunkLims(:,2) = chunkLims(:,1) + chunkSize - 1;
    chunkLims( chunkLims(:,1) > totScan, : ) = [];
    if chunkLims(end,2) > totScan, chunkLims(end,2) = totScan; end
    if chunkLims(end,2) > lastScan, chunkLims(end,2) = lastScan; end
    Nchunk = size(chunkLims, 1);
elseif ~isnan(Nchunk)
    chunkLims = round(linspace(firstScan,lastScan,Nchunk+1))';
    chunkLims(:,2) = circshift(chunkLims(:,1), -1) - 1;
    chunkLims(end,:) = []; 
    chunkLims(end) = lastScan;
end
chunkLength = diff(chunkLims,1,2)+1;
%sum(chunkLength)
end

