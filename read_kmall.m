function y=read_kmall(datadir,filename,filename2,startDgm,endDgm,outdir)
% function y=read_kmall(datadir,filename,filename2)
%    reads an EM2040 file, collects metadata, saves metadata and data into
%       structure and file
%    input:
%       datadir = directory in which file(s) reside
%       filename = name of file containing the kmwcd data
%       filename2 = name of kmall file with other data
%       startDgm = datagram (ping) count to start at (can be [])
%       endDgm = datagram (ping) count to end at (can be [])
%

% set kmwcd file name
fname = fullfile(datadir,filename);

% get basic information about file but may not need this
%KMALLfileinfo = CFF_kmall_file_info(fname);

% read metadata from kmwcd file
KMALLdata = CFF_read_kmall(fname);
wcdat=KMALLdata.EMdgmMWC;  % water column metadata

% read metadata from kmall file if needed
if ~isempty(filename2)
    fname2 = fullfile(datadir,filename2);
    KMALLdata2 = CFF_read_kmall(fname2);
    KMALLdata2.EMdgmMRZ
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
Nrx=arrayfun(@(x) x.rxInfo.numBeams,wcdat);
if abs(mean(diff(Nrx)))>0
    fprintf('WARNING: inconsistent number of beams across datagrams\n')
end
posixtimes=arrayfun(@(x) x.header.time_sec,wcdat);
tnanosec=arrayfun(@(x) x.header.time_nanosec,wcdat);
numOfDgms= arrayfun(@(x) x.partition.numOfDgms,wcdat);
dgmNum= arrayfun(@(x) x.partition.dgmNum,wcdat);
pingcount=arrayfun(@(x) x.cmnPart.pingCnt,wcdat);
numbeams=arrayfun(@(x) x.rxInfo.numBeams,wcdat);
numsectors=arrayfun(@(x) x.txInfo.numTxSectors,wcdat);
    maxNumSectors=max(numsectors);
fprintf('ping labels run from %d to %d \n',pingcount(1),pingcount(end))

% check that which datagrams to pull is already set (and set if it is not)
% default if not pre-set, is to pull all datagrams
if isempty(startDgm)
    startDgm=1;
end
if isempty(endDgm)
    Ndgm=totalNdgm-startDgm+1; 
    endDgm=Ndgm;
else
    Ndgm=endDgm-startDgm+1;
end

% need to know buffer size in order to do all setup
numSamps1=wcdat(1).beamData_p.numSampleData; % not sure this is correct
maxWCSampIdx1=numSamps1(1);
    fprintf('number samples (ping 1, beam 1) = %d\n',maxWCSampIdx1)
maxSamps=zeros(Ndgm,1);
minSamps=zeros(Ndgm,1);
keepSamps=zeros(Nrx(1),totalNdgm);
for i=1:totalNdgm
    numSamps=wcdat(i).beamData_p.numSampleData;
    maxSamps(i)=max(numSamps);
    minSamps(i)=min(numSamps);
    keepSamps(:,i)=numSamps;
end
maxWCSampIdx=max(maxSamps(startDgm:endDgm));
    fprintf('maximum number samples for all pings = %d\n',max(maxSamps))
    fprintf('maximum number samples for pings to read = %d\n',maxWCSampIdx)
    fprintf('maximum variation in number of samples within a single ping = %d \n',max(maxSamps-minSamps))


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
% this section will need more info output if any split datagrams are
% detected 


% structure for rest 
% loop over all datagrams (or all that want to read)
%     read amplitude data
%     store data in wcdat
%     store any metadata pulled from kmall rather than kmwcd in wcdat
% maybe end loop
% save data

% setup vectors for storage
SoundSpeed=zeros(Ndgm,1);
SampFreq=zeros(Ndgm,1);
TVGFuncApplied=zeros(Ndgm,1);
TVGOffset=zeros(Ndgm,1);
TxBeamWidth=zeros(Ndgm,maxNumSectors);
startRangeSampNum=zeros(Ndgm,max(numbeams));

% loop over set datagrams
for idgm=startDgm:endDgm

    % first get all the 
    % txInfo - get the number of sectors and set the center sector
    % sectorData
        % EM2040 tranmits in multiple sectors -- need to check how actually
        % set up - could be a 3-sector or 4-sector or other    
        NumSectors=wcdat(idgm).txInfo.numTxSectors;
        switch NumSectors
            case 1
                cenSec=1;
            case 2
                cenSec=1;
            case 3
                cenSec=2;
            otherwise
                fprintf('unexpected number of sectors\n')
                cenSec=1;  %arbitrary choice just so set
        end
        % also need the transmit beamwidth
        for isec=1:NumSectors
            TxBeamWidth(idgm,isec)=wcdat(idgm).sectorData(isec).txBeamWidthAlong_deg;
        end
 
% rxInfo 
%  already have Nrx = numbeams from above
%  could probably move all of this up there too (out of loop!)
    SoundSpeed(idgm)=wcdat(idgm).rxInfo.soundVelocity_mPerSec;
    SampFreq(idgm)=wcdat(idgm).rxInfo.sampleFreq_Hz;
    TVGFuncApplied(idgm)=wcdat(idgm).rxInfo.TVGfunctionApplied;
    TVGOffset(idgm)=wcdat(idgm).rxInfo.TVGoffset_dB;

% beamData_p
    startRangeSampNum(idgm,:)=wcdat(idgm).beamData_p.startRangeSampleNum;

 %   numSamps=wcdat(idgm).beamData_p.numSampleData; % not sure this is correct
%    xmitSectNum=wcdat(idgm).beamData_p.beamTxSectorNum;

% read beamAmp from binary file
%    beamAmpdata=read_bin_kmall(fname,wcdat(idgm));
%    beamAngle(idgm,:)=wcdat(idgm).beamData_p.beamPointAngReVertical_deg';
%    DR=wcdat(idgm).beamData_p.detectedRangeInSamples; % no idea if this makes sense
%    beamAmp=zeros(Nrx,1);
%    for ibeam=1:Nrx
%        tempBeamAmp=beamAmpdata(ibeam).sampleAmplitude05dB_p;
%        beamAmp(ibeam,1:numSamps(ibeam))=tempBeamAmp;
%    end
%    wcdat(idgm).beamData_p.beamAmp=beamAmp;

% probably need to pad out the beamAmp records here? or not? let's try not
% for now
% some thoughts from google/MatlabCentral on better padding operation:
%   >> D = {[1,2,3,4],[5,6],[7,8,9]};
%   >> M = zeros(5,5);
%   >> for k = 1:numel(D), M(k,1:numel(D{k})) = D{k}; end
%   >> M
%   M =
%     1  2  3  4  0
%     5  6  0  0  0
%     7  8  9  0  0
%     0  0  0  0  0
%     0  0  0  0  0
% this requires knowing the largest vector length but I already get and
% need that 
% I think it makes more sense than 
%   beamAmp = [beamAmp zeros(Nrx,maxWCSampIdx-length(beamAmp(1,:)))-999];
% which I think assumes that all beams have same number of samples
%
% should be able to implement new method above when pull beamAmp data
% originally

end

% store all data in structure so can manipulate elsewhere
% extract filecode for building output filenames
[~,filecode,~]=fileparts(fname);
% save number of sample information
sampinfofile=fullfile(outdir,['sampinfo_' filecode '.mat']);
save(sampinfofile,'maxSamps','minSamps','keepSamps')
% save key metadata to see if changes ever
keymetafile=fullfile(outdir,['keymeta_' filecode '.mat']);
save(keymetafile,'SoundSpeed','SampFreq','Nrx','TVGFuncApplied',...
    'TVGOffset','TxBeamWidth','startRangeSampNum')
% store necessary metadata for gridding
gridmeta.SoundSpeed=SoundSpeed; % 1-D vector
gridmeta.SampFreq=SampFreq; % 1-D vector
%gridmeta.beamAngle=beamAngle; % 2-D grid
gridmeta.numBeams=Nrx; % 1-D vector
% RxBeamWidth   degrees
% TxBeamWidth   degrees
% TVGFuncApplied
% TVGOffset
% beamAmp

