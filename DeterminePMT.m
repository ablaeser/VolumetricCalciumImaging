function [usePMTind, usePMTname] = DeterminePMT(requestColor, sbxInfo)
% Which PMTs are available (PMT gain > 0)?
PMTname = {'green','red'};
pmtGain = [sbxInfo.config.pmt0_gain, sbxInfo.config.pmt1_gain];
pmtOn = pmtGain > 0; 
pmtAvail = PMTname(pmtOn);
% Which PMTs were requested?
if ~isempty(pmtAvail)
    if strcmpi(requestColor, 'both')
        pmtRequest = {'green','red'};
        if sbxInfo.nchan == 1, pmtRequest = pmtRequest{1}; end
    else
        pmtRequest = requestColor;
    end
    usePMTind = find(strcmpi(pmtRequest, pmtAvail));
    if isempty(usePMTind) 
        warning('Requested PMT (%s) not available, using %s instead!', requestColor, pmtAvail);
        pmtRequest = pmtAvail;
        usePMTind = find(strcmpi(pmtRequest, pmtAvail));
    end
else
    error('No PMTs available!');
end
if numel(usePMTind) == 1
    usePMTname = PMTname{usePMTind};
else
    usePMTname = 'both';
end


%{
[chan, pmtOn] = DetermineChanSBX(sbxInfo)

%find();
Non = sum(pmtOn); %numel(pmtOn);
pmtName = {'green','red'};
if sbxInfo.nchan == 2
    chan = 'both';
elseif sbxInfo.nchan == 1 && Non == 1
    chan = pmtName{pmtOn};
else
    chan = 'green'; % if boht PMTs were on, but only one channel was saved, assume it was green
end
%}
end

