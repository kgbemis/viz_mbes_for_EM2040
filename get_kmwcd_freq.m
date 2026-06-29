function get_kmwcd_freq(datadir,filename,filename2)
% function y=read_kmall(datadir,filename,filename2)
%    reads an EM2040 file, collects metadata, saves metadata and data into
%       structure and file
%    input:
%       datadir = directory in which file(s) reside
%       filename = name of file containing the kmwcd data
%       filename2 = name of kmall file with other data
%

% set kmwcd file name
fname = fullfile(datadir,filename);

% get basic information about 
%KMALLfileinfo = CFF_kmall_file_info(fname);

% read metadata from kmwcd file
KMALLdata = CFF_read_kmall(fname);
wcdat=KMALLdata.EMdgmMWC;  % water column metadata

% read metadata from kmall file if needed
if ~isempty(filename2)
    fname2 = fullfile(datadir,filename2);
    KMALLdata2 = CFF_read_kmall(fname2);
end % isempty 

% check and print basic file information
fprintf('reading KMWCD file: %s \n',fname)
fprintf('start date = %s \n',...
    datestr(datetime(wcdat(1).header.time_sec,'ConvertFrom','posixtime')))
fprintf('end date = %s \n',...
    datestr(datetime(wcdat(end).header.time_sec,'ConvertFrom','posixtime')))
fprintf('central frequency of center sector = %4.0f kHz \n',...
    wcdat(1).sectorData(1).centreFreq_Hz./1000)

% pull most of setup data
totalNdgm=length(wcdat);
ckNrx=arrayfun(@(x) x.rxInfo.numBeams,wcdat);
if abs(mean(diff(ckNrx)))>0
    fprintf('WARNING: inconsistent number of beams across datagrams\n')
end
posixtimes=arrayfun(@(x) x.header.time_sec,wcdat);
tnanosec=arrayfun(@(x) x.header.time_nanosec,wcdat);
numOfDgms= arrayfun(@(x) x.partition.numOfDgms,wcdat);
dgmNum= arrayfun(@(x) x.partition.dgmNum,wcdat);
pingcount=arrayfun(@(x) x.cmnPart.pingCnt,wcdat);
numbeams=arrayfun(@(x) x.rxInfo.numBeams,wcdat);
fprintf('ping labels run from %d to %d \n',pingcount(1),pingcount(end))

% set which datagrams will pull
% for now do full file (ok for UNH exp data)
Ndgm=totalNdgm; 
startDgm=1;
endDgm=Ndgm;

% need to know buffer size in order to do all setup
numSamps1=wcdat(1).beamData_p.numSampleData; % not sure this is correct
maxWCSampIdx1=numSamps1(1);
    fprintf('number samples (ping 1, beam 1) = %d\n',maxWCSampIdx1)
maxSamps=zeros(Ndgm,1);
keepSamps=zeros(ckNrx(1),totalNdgm);
for i=1:totalNdgm
    numSamps=wcdat(i).beamData_p.numSampleData;
    maxSamps(i)=max(numSamps);
    keepSamps(:,i)=numSamps;
end
maxWCSampIdx=max(maxSamps(startDgm:endDgm));
    fprintf('maximum number samples for all pings = %d\n',max(maxSamps))
    fprintf('maximum number samples for pings to read = %d\n',maxWCSampIdx)


% check if any datagrams are split
DatagramNum=zeros(totalNdgm,1);
NumDatagrams=zeros(totalNdgm,1);
part_count=0;
split_count=0;
for idgm=1:totalNdgm 
    % partition - has the datagram number and split info
    DatagramNum(idgm)=wcdat(idgm).partition.dgmNum;
    NumDatagrams(idgm)=wcdat(idgm).partition.numOfDgms;
    if DatagramNum ~= NumDatagrams
        fprintf('datagram %d is %dth partition of a datagram \n',idgm,DatagramNum)
        part_count=part_count+1;
    end
    if NumDatagrams>1
        fprintf('WARNING: may be a split datagram at %d \n',idgm)
        split_count=split_count+1;
    end
end
if split_count~=0 & part_count~=0
    fprintf('Potentially %d split datagrams\n',max([split_count part_count]))
else
    fprintf('No split datagrams detected\n')
end

% nothing more as this file just to spit out basic information like
% frequency


