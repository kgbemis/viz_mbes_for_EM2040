% play with reading the kmall files

% used Liz's parsing codes as a start but have pulled in files from the
% CoFFee code to read the new format and better grasp structure of files

%% manual entry of parameters
% sonar parameters
%maxWCSampIdx = 1060; % setting a fixed watercolumn buffer in meters 
%           (this is 600 in Liz's code)
%       maxWCSampIdx needs to be set for the particular data set
%       outer beams should have maximum number of samples based on plotted
%       patterns for number of samples so can use the following after
%       metadata is read:
%           numSamps=wcdat(idgm).beamData_p.numSampleData;
%           maxWCSampIdx=numSamps(1);
RxBeamWidth=1; % setting receive beamwidth arbitrarily (need to find where to read)
cartspeed=0.05; % setting velocity for along track direction (0.34 max at UNH)

% operations choices
dotest=1;
describe_datagrams=1;
% always reads kmall file info and kmwcd data
check_dgmtimes=1; % always computes times, this just controls plots
vizping=1;
do3Dviz=0;
savefiles=1;

outdir='../MBES_mat_files/';
%expcode='geoA_Apr15_single_15Hz_slow_hot';
expcode='';

%% open an EM data file 

filelocation = ['../MBES_raw_data/' expcode '/'];
%filelocation = 'e:UNH_tank_experiment\EM2040\UNH_EM2040_40_Apr_13_2023\';
%filelocation = 'e:UNH_tank_experiment\EM2040\UNH_EM2040_40_Apr2023_geoA\';
%filename='0009_20200918_094230.kmall';
%filename='0010_20200918_095915.kmall';
%filename='0004_20230406_123411.kmall'; % apr 6 - flow setup issues
%filename='0004_20230406_123411.kmwcd';  %   too long (3407 pings) to read 
                                        %   in one gulp
%filename='0017_20230405_163521.kmwcd'; % apr 5 - hot flow good
%filename='0005_20230413_141335.kmwcd';
%filename='0006_20230413_143146.kmwcd';
%filename='0007_20230413_143346.kmwcd';
%filename='0008_20230413_143848.kmwcd';
%filename='0009_20230413_143947.kmwcd';
%filename='0010_20230413_144205.kmwcd';
%filename='0011_20230413_144251.kmwcd';
%filename='0000_20230413_160721.kmwcd';
%filename='0001_20230413_160826.kmwcd';
%filename='0002_20230413_161101.kmwcd';
%filename='0003_20230413_161138.kmwcd';
%filename='0004_20230413_161651.kmwcd';


%filename='0017_20230405_163521.kmwcd'; % early run


%filename='0020_20230413_171151.kmwcd'; % apr13-575Hz
%filename='0021_20230413_171231.kmwcd'; % apr13-B-hot flow-300 kHz
%filename='0022_20230413_172311.kmwcd'; % apr13-B-hot flow-200 kHz-gas 
%filename='0023_20230413_172348.kmwcd'; % apr13-B-hot flow-200 kHz-gas 
%filename='0024_20230413_172544.kmwcd'; % apr13-320Hz 
%filename='0024_20230413_172544.kmwcd'; % apr13-B-hot flow-600 kHz
%filename='0025_20230413_172618.kmwcd'; % apr13-B-
%filename='0026_20230413_174909.kmwcd'; % apr13-? flow-? kHz
%filename='0027_20230413_175000.kmwcd'; titlestr='apr13-B-? flow-200kHz';
% try 46 (geoC 400kHz looks like good plume)
%filename='0046_20230413_180331.kmwcd';
% try 48
%filename='0048_20230413_180417.kmwcd';
%filename='0053_20230413_181816.kmwcd';
filename='0075_20230414_193631.kmwcd'; % Liz got plume but I don't :-(
filename2='0075_20230414_193631.kmall'; % Liz got plume but I don't :-(
%filename='0077_20230413_190346.kmwcd';
%filename='0078_20230413_190400.kmwcd';
%filename='0079_20230413_190422.kmwcd';
%filename='0081_20230413_190450.kmwcd';
%filename='0065_20230413_184653.kmwcd';

%filename='0024_20230414_163451.kmwcd'; % geoA 400hz
%filename='0023_20230414_163427.kmwcd'; 
%filename='0046_20230414_175425.kmwcd';
%filename='0048_20230414_175530.kmwcd';

%filename='0000_20230415_123519.kmwcd'; % geoA hot slow 200Hz
%filename='0001_20230415_123600.kmwcd'; % geoA hot slow 200Hz
% no 0002 file
%filename='0003_20230415_124320.kmwcd'; % geoA hot slow 200Hz
%filename='0004_20230415_124403.kmwcd'; % geoA hot slow 200Hz
%filename='0005_20230415_124611.kmwcd'; % geoA hot slow 300Hz
%filename='0006_20230415_124653.kmwcd'; % geoA hot slow 300Hz
%filename='0007_20230415_124807.kmwcd'; % geoA hot slow 400Hz
%filename='0008_20230415_124857.kmwcd'; % geoA hot slow 400Hz

%filename='0000_20230415_193500.kmwcd'; % calibration 
%filename='0003_20230415_200249.kmwcd';

%filename='0004_20230415_204544.kmwcd'; % test file for student work
%filename2='0004_20230415_204544.kmall'; % absolute calibration speed test - different setup

%filename='0017_20230405_163521.kmwcd'; % early run
%filename2='0017_20230405_163521.kmall'; 

fname = fullfile(filelocation,filename);
fprintf('reading file: %s \n',fname)

% open selected file
fid = fopen(fname,'r');

%% test reading the data file
% this part of code copied from CoFFee
%   Authors: Alex Schimel (NGU, alexandre.schimel@ngu.no) and Yoann
%   Ladroit (NIWA, yoann.ladroit@niwa.co.nz)
%   2017-2021; Last revision: 20-08-2021

if dotest
        fprintf('Initial header info in file:\n')
    % Datagram length in bytes. The length field at the start (4 bytes) 
    % and end of the datagram (4 bytes) are included in the length count.
    out_struct.numBytesDgm = fread(fid,1,'uint32');
        fprintf('size of next datagram: %d\n',out_struct.numBytesDgm)

    % Multi beam datagram type definition, e.g. #AAA
    out_struct.dgmType = fscanf(fid,'%c',4);
        fprintf('datagram type: %s\n',out_struct.dgmType)

    % Datagram version.
    out_struct.dgmVersion = fread(fid,1,'uint8');
        fprintf('datagram version: %d\n',out_struct.dgmVersion)

    % System ID. Parameter used for separating datagrams from different
    % echosounders if more than one system is connected to SIS/K-Controller.
    out_struct.systemID = fread(fid,1,'uint8');
        fprintf('system ID: %d\n',out_struct.systemID)

    % Echo sounder identity, e.g. 124, 304, 712, 2040, 2045 (EM 2040C)
    out_struct.echoSounderID = fread(fid,1,'uint16');
        fprintf('echo sounder identity: %d\n',out_struct.echoSounderID)

    % UTC time in seconds. Epoch 1970-01-01. time_nanosec part to be added for
    % more exact time.
    out_struct.time_sec = fread(fid,1,'uint32');
        fprintf('time seconds: %d\n',out_struct.dgmVersion)

    % Nano seconds remainder. time_nanosec part to be added to time_sec for
    % more exact time.
    out_struct.time_nanosec = fread(fid,1,'uint32');
        fprintf('time nanoseconds: %d\n',out_struct.time_nanosec)

        fprintf('\n')
    % close file as next steps will reopen
    fclose(fid);

end % dotest

%% Call CFF function to get datagram structure in kmall file
KMALLfileinfo = CFF_kmall_file_info(fname);
%   * |KMALLfileinfo|: structure about datagrams in KMALLfilename, with fields:
%     * |file_name|: input file name
%     * |fileSize|: file size in bytes
%     * |dgm_num|: number of datagram in file
%     * |dgm_type_code|: datagram type as string, e.g. '#IIP' (Kongsberg .kmall format)
%     * |dgm_type_text|: datagram type description (Kongsberg .kmall format) 
%     * |dgm_type_version|: version for this type of datagram, as int (Kongsberg .kmall format)
%     * |dgm_counter|: counter for this type and version of datagram in the
%           file. There should not be multiple versions of a same type in
%           a same file, but we never know...
%     * |dgm_start_pif|: position of beginning of datagram in file 
%     * |dgm_size|: datagram size in bytes
%     * |dgm_sys_ID|: System ID. Parameter used for separating datagrams from different echosounders.
%     * |dgm_EM_ID|: Echo sounder identity, e.g. 124, 304, 712, 2040, 2045 (EM 2040C)
%     * |sync_counter|: number of bytes found between this datagram and the
%           previous one (any number different than zero indicates a sync error)
%     * |date_time|: datagram date in datetime format
%     * |parsed|: flag for whether the datagram has been parsed. Initiated
%               at 0 at this stage. To be later turned to 1 for parsing.
% structure also has
%     * |list_dgm_type|: this has the unique dgm_type_codes
%     * |list_dgm_counter|: not sure if counts of instances?

%% explore what is in the kmall file
if describe_datagrams
    figure(3)
    plot(KMALLfileinfo.dgm_num,'.')
    fprintf('min dgm_num= %d; max dgm_num=%d\n',...
        min(KMALLfileinfo.dgm_num),max(KMALLfileinfo.dgm_num))
    fprintf('min sync_counter= %d; max sync_counter=%d\n',...
        min(KMALLfileinfo.sync_counter),max(KMALLfileinfo.sync_counter))
    fprintf('min parsed= %d; max parsed=%d\n',...
        min(KMALLfileinfo.parsed),max(KMALLfileinfo.parsed))

    figure(4)
    b1=barh(KMALLfileinfo.list_dgm_counter,...
        'FaceColor',[0.7 0.9 0.7],'EdgeColor',[0.7 0.9 0.7]);
    yticklabels(KMALLfileinfo.list_dgm_type)
    ylabel('Datagram Types')
    xlabel('Counts')
    datalabels=string(b1(1).YData);
    datalabelpos=max(KMALLfileinfo.list_dgm_counter)+2;
    potdatalabels={'IIP','IOP','SVP','SKM','MWC','SVT','CHE'};
    fullinfolabellist={'IIP=Installation parameters and sensor setup',...
        'IOP= Runtime parameters as chosen by operator',...
        'SVP= sensor data from sound velocity profile or CTD',...
        'SKM= sensor KM binar sensor format',...
        'MWC= Multibeam water column datagram', ...
        'SVT= sensor data for sound velocity at transducer',...
        'CHE= Compatibility data for heave',...
        'SPO= sensor data for position',...
        'CPO= Compatibility data for position',...
        'MRZ= Multibeam raw range and depth datagram',...
        'SCL= sensor data from clock'};
    %infolabels

    if sscanf(version('-release'),'%d')>2019
    text(datalabelpos*ones(size(yticklabels)),b1(1).XEndPoints,datalabels,'VerticalAlignment','middle')
    %text(2000*ones(size(yticklabels)),b1(1).XEndPoints,infolabels,'VerticalAlignment','middle')
    else
    text(datalabelpos*ones(size(yticklabels)),b1(1).XData,datalabels,'VerticalAlignment','middle')
    %text(2000*ones(size(yticklabels)),b1(1).XData,infolabels,'VerticalAlignment','middle')
    end
end % describe_datagrams

%% now we want to read some datagrams and figure out what to do with the data!
% this should read all the data
KMALLdata = CFF_read_kmall(fname);
% 5/31/23 - try making a cell with both .kmwcd and .kmall names
%fnamepair{1,1}=fname;
%fnamepair{2,1}=fname2;
%KMALLdata = CFF_read_kmall(fnamepair);
% the second version will only work if I upgrade the version to higher than
% 2020a

% now have a set of stuctures in KMALLdata that contain the datagrams for
% each type plus an info strucure

%% get the non-water column information

if ~isempty(filename2)
    fname2 = fullfile(filelocation,filename2);
    KMALLdata2 = CFF_read_kmall(fname2);

end % isempty 

%% start parsing out the water column data

wcdat=KMALLdata.EMdgmMWC;  % water column data

% next step - convert datagram time to matlab datenumber
fprintf('start date = %s \n',...
    datestr(datetime(wcdat(1).header.time_sec,'ConvertFrom','posixtime')))
fprintf('end date = %s \n',...
    datestr(datetime(wcdat(end).header.time_sec,'ConvertFrom','posixtime')))
fprintf('central frequency of center sector = %f \n',wcdat(1).sectorData(1).centreFreq_Hz)

% Note: might be able to replace this section with with four lines that
% just index into the single number, e.g. "wcdat.header.time_sec"   
% Second Note: pingcount is really an incremeted label that may or may not
% start at 1 ==> subtract to get number of pings between two times
Ndgm=length(wcdat);
%%dgmtimes=NaT(Ndgm,1);
%numOfDgms=zeros(Ndgm,1);
%dgmNum=zeros(Ndgm,1);
%posixtimes=zeros(Ndgm,1);
%pingcount=zeros(Ndgm,1);
%tnanosec=zeros(Ndgm,1);
%ckNrx=zeros(Ndgm,1);
%for i=1:Ndgm
%    dgmtimes(i)=datetime(wcdat(i).header.time_sec,'ConvertFrom','posixtime');
%    posixtimes(i)=wcdat(i).header.time_sec;
%    tnanosec(i)=wcdat(i).header.time_nanosec;
%    numOfDgms(i)= wcdat(i).partition.numOfDgms;
%    dgmNum(i)= wcdat(i).partition.dgmNum;
%    pingcount(i)=wcdat(i).cmnPart.pingCnt;
%    ckNrx(i)=wcdat(i).rxInfo.numBeams;
%end
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
%
ckNumSectors=arrayfun(@(x) x.txInfo.numTxSectors,wcdat);
dgmtimes=arrayfun(@(x) datetime(x.header.time_sec,'ConvertFrom','posixtime'),wcdat);
if check_dgmtimes
    figure(5)
    subplot(211)
    plot(dgmtimes,numOfDgms,'p',dgmtimes,dgmNum)
    subplot(212)
    plot(dgmtimes,pingcount,'.-')
    
    fprintf('first ping count = %d\n',pingcount(1))
    fprintf('last ping count = %d\n',pingcount(end))
    fprintf('total counts = %d\n',pingcount(end)-pingcount(1)+1)
end

figure(11)
p = arrayfun(@(a) plot(a.beamData_p.beamTxSectorNum),wcdat);
xlabel('beam number')
ylabel('sector number')


%% section loop over pings out
% if numOfDgms always 1, can just proceed
% assuming this for now

% setting these above now
% setting a fixed watercolumn buffer
%   maxWCSampIdx = 600;
% setting receive beamwidth arbitrarily (need to find where to read)
%   RxBeamWidth=1;
% setting velocity for along track direction
%   cartspeed=0.34;

% need to know buffer size in order to do all setup
numSamps1=wcdat(1).beamData_p.numSampleData; % not sure this is correct
maxWCSampIdx1=numSamps1(1);
    fprintf('number samples (ping 1, beam 1) = %d\n',maxWCSampIdx1)
maxSamps=zeros(Ndgm,1);
for i=1:Ndgm
    numSamps=wcdat(i).beamData_p.numSampleData;
    maxSamps(i)=max(numSamps);
end
maxWCSampIdx=max(maxSamps);
    fprintf('maximum number samples = %d\n',maxWCSampIdx)


%dgmtime=NaT(num,1); would create empty datetime array
YY=zeros(ckNrx(1),maxWCSampIdx,Ndgm);
ZZ=zeros(ckNrx(1),maxWCSampIdx,Ndgm);
SV=zeros(ckNrx(1),maxWCSampIdx,Ndgm);
TT=zeros(ckNrx(1),maxWCSampIdx,Ndgm);
XX=zeros(ckNrx(1),maxWCSampIdx,Ndgm);
xBottom=zeros(Ndgm,ckNrx(1));
dTime=zeros(Ndgm,1);

for idgm=1:Ndgm  % just do the first ping for now
    %fprintf('reading %d-th ping now\n',idgm)
    % header
        % time coming from datagram header - already pulled
        %dgmtime(idgm)=datetime(wcdat(idgm).header.dgtime,'ConvertFrom','posixtime');
        %thistime=second(dgmtimes(idgm))+tnanosec(idgm)/1e9;
        thistime=posixtimes(idgm)+tnanosec(idgm)/1e9;
    % partition - has the datagram number and split info
    DatagramNum=wcdat(idgm).partition.dgmNum;
    NumDatagrams=wcdat(idgm).partition.numOfDgms;
    % cmnPart - not sure if will need this info
    % txInfo - not sure if will need this info
        NumSectors=wcdat(idgm).txInfo.numTxSectors;
        switch NumSectors
            case 1
                cenSec=1;
            case 2
                cenSec=1;
            case 3
                cenSec=2;
        end
    % sectorData - not sure if will need this info
        % EM2040 tranmits in multiple sectors -- need to check how actually
        % set up - could be a 3-sector or 4-sector or other
        TxBeamWidth=zeros(NumSectors,1);
        for isec=1:NumSectors
            TxBeamWidth(isec)=wcdat(idgm).sectorData(isec).txBeamWidthAlong_deg;
        end
    % rxInfo - not sure if will need this info
    SoundSpeed=wcdat(idgm).rxInfo.soundVelocity_mPerSec;
    SampFreq=wcdat(idgm).rxInfo.sampleFreq_Hz;
    TVGFuncApplied=wcdat(idgm).rxInfo.TVGfunctionApplied;
    TVGOffset=wcdat(idgm).rxInfo.TVGoffset_dB;
        %RxBeamWidth=wcdat(iping).rxInfo.nothere    beamAngle=wcdat(i).beamData_p.beamPointAngReVertical_deg;
    startRangeSampNum=wcdat(idgm).beamData_p.startRangeSampleNum;
    numSamps=wcdat(idgm).beamData_p.numSampleData; % not sure this is correct
    % maxWCSampIdx=numSamps(1); % actually need this to be the max, max 
    %       so set it outside this loop
    xmitSectNum=wcdat(idgm).beamData_p.beamTxSectorNum;
    %beamNum=wcdat(i).beamData_p.?  ask Liz what this is
    %beamAmp - need to read from binary file still
    beamAmpdata=read_bin_kmall(fname,wcdat(idgm));
    %DR huh?
    Nrx=wcdat(idgm).rxInfo.numBeams;
    beamAngle=wcdat(idgm).beamData_p.beamPointAngReVertical_deg';
    DR=wcdat(idgm).beamData_p.detectedRangeInSamples; % no idea if this makes sense
    beamAmp=zeros(Nrx,1);
    for ibeam=1:Nrx
        tempBeamAmp=beamAmpdata(ibeam).sampleAmplitude05dB_p;
        beamAmp(ibeam,1:numSamps(ibeam))=tempBeamAmp;
    end
    wcdat(idgm).beamData_p.beamAmp=beamAmp;
    %wcdat(idgm).beamData_p.beamPhase=

        % this part is very much Liz's code
    if DatagramNum == NumDatagrams
        pingidx = idgm;  % since we know ping numbers, just set rather than counting

        pingTime(pingidx) = thistime;    % datetime
        %elapsedsec=second(pingTime(idgm)-pingTime(1));
        elapsedsec=pingTime(idgm)-pingTime(1);
        %elapsedsec=posixtimes(idgm)-posixtimes(1);

        %size(beamAmp)
        %size(zeros(Nrx,maxWCSampIdx-length(beamAmp(1,:)))-999)
        beamAmp = [beamAmp zeros(Nrx,maxWCSampIdx-length(beamAmp(1,:)))-999];

        %disp({'SoundSpeed=' SoundSpeed});
        %disp({'SampFreq=' SampFreq});

        % range in meters
        %range =  (1:(length(beamAmp(1,:))))*SoundSpeed/10/2/(SampFreq/100);
        range =  (1:(length(beamAmp(1,:))))*SoundSpeed/2/SampFreq;

        % pad with zeros out to the maximum range
        range = [range zeros(1,maxWCSampIdx-length(range))];   

        % this is depth
        %z = cos(beamAngle'/100*pi/180)*range;  
        z = cos(beamAngle*pi/180)*range;  % my angle in straight degrees

        % this is across-trackd distance
        %y = sin(beamAngle'/100*pi/180)*range;
        y = sin(beamAngle*pi/180)*range;

        Awc = beamAmp/2;

        X = TVGFuncApplied;
        C = TVGOffset;
        clear TS
        %RTval=10*log10((RxBeamWidth/10)*pi/1800*(TxBeamWidth/10)*pi/1800);
        RTval=10*log10(RxBeamWidth*pi/180*TxBeamWidth(cenSec)*pi/180);
        TS = Awc + RTval.*ones(size(Awc)) - double(X).*log10(ones(length(Awc(:,1)),1)*range) + 40*log10(ones(length(Awc(:,1)),1)*range) - double(C);
        tsBuf1(pingidx,1:length(TS(:,1)),1:length(TS(1,:))) = TS;

        %steeringangle = beamAngle/100;
        steeringangle = beamAngle;
        recieveAngle = 1./cos(steeringangle*pi/180);

        %RxRad = (RxBeamWidth/10)*pi/1800;
        RxRad = RxBeamWidth*pi/180;
        Length = 2*range*sin(RxRad/2);

        %TxRad = (TxBeamWidth/10)*pi/1800*recieveAngle;
        TxRad = TxBeamWidth(cenSec)*pi/180*recieveAngle;
        %Width = 2*ones(length(TxRad),length(range)).*range.*(sin(TxRad./2)).';
        Width = 2*ones(length(TxRad),length(range)).*range.*(sin(TxRad./2));
        BeamArea = Length.*Width;

        Tau = 3500/1e5;
        Sv=zeros(size(TS));
        for ii = 1:size(TS,2)

            Vol_log = 10*log10(BeamArea(:,ii)*Tau*SoundSpeed/2);
            Sv(:,ii) = TS(:,ii) - Vol_log;
        end

        for aa = 1:Nrx
            if DR(aa) ~= 0
                yBottom(pingidx,aa) = y(aa,DR(aa));
                zBottom(pingidx,aa) = z(aa,DR(aa));
            else
                yBottom(pingidx,aa) = 0;
                zBottom(pingidx,aa) = 0;
            end
        end


        [~,nadir_idx] = min(abs(yBottom(pingidx,:)-35));

        if zBottom(nadir_idx) ~= 0
            nadir_idx = nadir_idx;
        elseif zBottom(nadir_idx + 1) ~= 0
            nadir_idx = nadir_idx + 1;
        else
            nadir_idx = nadir_idx - 1;
        end
        
        if vizping
             figure(1)
             subplot(211)
             pcolor(y,z,TS); 
                shading flat; axis equal; 
                set(gca,'ydir','reverse');
             set(gca,'fontname','Times'); 
                caxis([-140 -60]); colorbar
             title(['Ping ' num2str(pingidx) ': TS'])
             ylim([0 6])
             hold on
             plot(yBottom(pingidx,:),zBottom(pingidx,:),'r.')
             hold off 

             subplot(212)
             pcolor(y,z,Sv); 
                shading flat; axis equal; 
                set(gca,'ydir','reverse');
             set(gca,'fontname','Times'); 
                caxis([-140 -60]); colorbar
             title(['Ping ' num2str(pingidx) ' :Sv'])
             hold on
             plot(yBottom(pingidx,:),zBottom(pingidx,:),'r.')
             ylim([0 6])
             hold off
        end % vizping

    %    cut_Sv  = Sv(:,1:nadir_idx);
    %    cut_y   = y(:,1:nadir_idx);
    %    cut_z   = z(:,1:nadir_idx);

   %     bin_edges = -150:10:150;

    %                 b_x = beam_x(pingidx,:);
    %                 b_y = beam_y(pingidx,:);
    %                 b_z = beam_z(pingidx,:);

    %    for iping = 1:(length(bin_edges)-1)
    %       bin_idx = cut_y < bin_edges(iping+1) & cut_y > bin_edges(iping);
    %          test_mean(i)    = mean(cut_Sv(bin_idx));
    %          test_max(i)     = max(cut_Sv(bin_idx));
    %        test1 = movmean(cut_Sv(bin_idx),10);
    %        test_mmax(iping) = max(test1);
    %        y_mid(iping) =  mean(cut_y(bin_idx));
    %    end 

    %    y_mid_bin(pingidx,:) = y_mid;
    %    z_nadir(pingidx) = zBottom(nadir_idx);
    %    Sv_mmax(pingidx,:) = test_mmax;

    YY(:,:,idgm)=y(:,1:maxWCSampIdx);
    ZZ(:,:,idgm)=z(:,1:maxWCSampIdx);
    SV(:,:,idgm)=Sv(:,1:maxWCSampIdx);
    TT(:,:,idgm)=TS(:,1:maxWCSampIdx);
    XX(:,:,idgm)=cartspeed*elapsedsec*ones(256,maxWCSampIdx);
    xBottom(idgm,:)=cartspeed*elapsedsec*ones(1,256);
    dTime(idgm)=elapsedsec;

    end % inside datagram split check

    %pause
end

%% summary visualizations
if do3Dviz
    pingtoplot=fix(Ndgm/2); % 180
    thisbeam=128; % 150
   
%plot along track profile
figure(6)
pcolor(squeeze(XX(thisbeam,:,:)),squeeze(ZZ(thisbeam,:,:)),squeeze(TT(thisbeam,:,:)))
shading flat
set(gca,'ydir','reverse');
set(gca,'fontname','Times'); 
if sscanf(version('-release'),'%d')<2022
    caxis([-140 -40]); colorbar
else
    clim([-140 -40]); colorbar    
end
title(['end Ping ' num2str(pingidx) ' :TS'])
hold on
 plot(xBottom(:,thisbeam),zBottom(:,thisbeam),'r.')
 plot(xBottom(pingtoplot,thisbeam).*ones(2,1),[0 6],'k')
 ylim([0 6])
hold off

figure(7)
pcolor(squeeze(XX(thisbeam,:,:)),squeeze(ZZ(thisbeam,:,:)),squeeze(SV(thisbeam,:,:)))
shading flat
set(gca,'ydir','reverse');
set(gca,'fontname','Times'); 
if sscanf(version('-release'),'%d')<2022
    caxis([-140 -40]); colorbar
else
    clim([-140 -40]); colorbar    
end
title(['end Ping ' num2str(pingidx) ' :SV'])
hold on
 plot(xBottom(:,thisbeam),zBottom(:,thisbeam),'r.')
 ylim([0 6])
hold off

% plot a summary of pings
figure(8)
pcolor(mean(YY(:,:,fix(Ndgm/2):ceil(3*Ndgm/4)),3),...
    mean(ZZ(:,:,fix(Ndgm/2):ceil(3*Ndgm/4)),3),...
    mean(TT(:,:,fix(Ndgm/2):ceil(3*Ndgm/4)),3))
shading flat
set(gca,'ydir','reverse');
set(gca,'fontname','Times'); 
if sscanf(version('-release'),'%d')<2022
    caxis([-140 -40]); colorbar
else
    clim([-140 -40]); colorbar    
end
title(['end Ping ' num2str(pingidx) ' :Sv'])
hold on
 plot(mean(yBottom,1),mean(zBottom,1),'r.')
 ylim([0 6])
hold off

% plot a specific ping
figure(9)
pcolor(YY(:,:,pingtoplot),ZZ(:,:,pingtoplot),TT(:,:,pingtoplot))
shading flat
set(gca,'ydir','reverse');
set(gca,'fontname','Times'); 
if sscanf(version('-release'),'%d')<2022
    caxis([-140 -40]); colorbar
else
    clim([-140 -40]); colorbar    
end
title(['end Ping ' num2str(pingidx) ' :Sv'])
hold on
 plot(mean(yBottom,1),mean(zBottom,1),'r.')
 plot(YY(thisbeam,:,pingtoplot),ZZ(thisbeam,:,pingtoplot),'k')
 ylim([0 6])
hold off

%figure(10)
%p = patch(isosurface(XX,YY,ZZ,SV,-80));
%isonormals(XX,YY,ZZ,SV,p)
%p.FaceColor = 'red';
%p.EdgeColor = 'none';
%daspect([1 1 1])
%view(3); 
%axis tight
%camlight 
%lighting 

end % do3Dviz


%% save files

if savefiles

% outdir='../MBES_mat_data/'; % moving to top
datafile=filename(1:end-6);
source=filename(end-4:end);
outmatfile=fullfile(outdir,[datafile '_' source '.mat']);
outvizfile=fullfile(outdir,[datafile '_' source '_viz.mat']);
outallstruct=fullfile(outdir,[datafile '_bothstructs.mat']);
fprintf('saving to file %s \n',outmatfile)
save(outmatfile,'wcdat')
save(outvizfile,'XX','YY','ZZ','SV','TT','xBottom','yBottom','zBottom','dTime');
save(outallstruct,'KMALLdata','KMALLdata2')
end % savefiles

fclose(fid);

return
%%
% every thing beyond here is Liz's - mostly using as guide to which
% datagrams are needed and other steps
%********************************************************************
%********************************************************************
%********************************************************************
%********************************************************************
%********************************************************************
%********************************************************************

%% start reading data file

%tic
pingidx = 0;
rangeSB = [];
BufSize = 1000;
tsBuf1 = zeros(BufSize,256,2500);
stopp = 0;
bufCt1 = 0;

while ~stopp

    sz = fread(fid,1,'uint32');     % the size of the next datagram. note: the size does not include this field
        fprintf('size of next datagram: %d\n',sz)
    startid = fread(fid,1,'uint8'); % the start id of a datagram should be 2
        fprintf('start id: %d\n',startid)
    
    if startid == 2
       
        datatype = fread(fid,1,'char'); % this is the type of datagram that we are reading
         fprintf('data type id is %d',datatype)
       
        % char(83) = S, seabead image data
        if datatype == 83
            [DateVal,TimeVal,PingCount,MeanAbsorp,PulseLength,NormIncRange,TVGstartRange,...
                TVGstopRange,NormincBSN,OblqincBSO,TxBeamWidth,TVGcrossoverAng,NumBeamsHere,...
                beamIndexNum,sortDir,NumSampsPerBeam,BeamCenterSampNum,SimradSnippets,BeamDetAmp] = readEM_S(fid);
        pause
        % char(82) = R, runtime parameters data    
        elseif datatype == 82
            [ModelNum,DateVal,TimeVal,PingCount,SystemSerNum,Mode,MinDepth,MaxDepth,AbsCoef,TxPulseLength,...
                TxBeamWidth, TxPower,RxBeamWidth,RxBandWidth,RxFixedGain, TVGXover,MaxPortCoverage, MaxStbdCoverage,...
                MaxPortSwath,BeamSpacing,MaxStbdSwath] = readEM_R(fid);
            
        elseif datatype == 78
            [DateVal,TimeVal,PingCount,SSPD,SampFreq,TiltAng,FocRange,SignalLength,SectorXmitDelay,CentFreq,MeanAbsorpt,SigWavID,XmitSecNum,BW,BeamAngle,XmitSecNumByBeam,TwoWayTT,DetInfo,BS] = readEM_N(fid);
            %pause
        
        % char(107) = k, watercolumn data    
        elseif datatype == 107
            [DateVal,TimeVal,PingCount,NumDatagrams,DatagramNum,NumXmitSect,...
                NumRecBeams,NumBeamsHereTemp, SoundSpeed, SampFreq, TxTimeHeave, TVGFuncApplied, TVGOffset,...
                centFreqTemp,tiltAngle, xmitSectNumTemp,beamAngleTemp,startRangeSampNumTemp,numSampsTemp,xmitSectNumRecTemp,...
                beamNumTemp,beamAmpTemp,DRtemp] = readEM_k(fid);
          
            % convert datagram time to matlab datenumber
            YY = floor(DateVal/10000);
            MM = floor( (DateVal-YY*10000)/100);
            DD = floor( (DateVal-YY*10000-MM*100));
            dn = datenum(YY,MM,DD)+TimeVal/1000/3600/24;
            
            % there may be more than one water column datagram due to the
            % size limit on a single datagram, and we want to stitch all of
            % the data for a single ping together.  The offset calculated
            % below allows the data to be associated with its proper beam
            % number
            NumBeamsHere(DatagramNum) = NumBeamsHereTemp;
            if DatagramNum > 1
                offsetidx = sum(NumBeamsHere(1:(DatagramNum-1)));
            else
                offsetidx = 0;
            end
            
            indx = 1:length(beamAngleTemp);
            beamAngle(offsetidx + indx) = beamAngleTemp;
            startRangeSampNum(offsetidx + indx) = startRangeSampNumTemp;
            numSamps(offsetidx + indx) = numSampsTemp;
            xmitSectNum(offsetidx + indx) = xmitSectNumRecTemp;
            beamNum(offsetidx + indx) = beamNumTemp;
            beamAmp(offsetidx + indx,1:length(beamAmpTemp(1,:))) = beamAmpTemp;
            DR(offsetidx + indx) = DRtemp;
            
            % if we have read all of the data for a single ping (plot this
            % data)
            if DatagramNum == NumDatagrams
                %disp( [num2str(floor(TimeVal/1000/3600)) ':' num2str(floor(rem(TimeVal/1000/3600,1)*60),'%2d') ':' num2str(rem(rem(TimeVal/1000/3600,1)*60,1)*60,'%2.2f')])

                pingidx = pingidx + 1

                pingTime(pingidx) = TimeVal/1000;    % time in seconds since midnight
                
                beamAmp = [beamAmp zeros(NumRecBeams,maxWCSampIdx-length(beamAmp(1,:)))-999];

                %disp({'SoundSpeed=' SoundSpeed});
                %disp({'SampFreq=' SampFreq});
                
                % range in meters
                range =  (1:(length(beamAmp(1,:))))*SoundSpeed/10/2/(SampFreq/100);
                
                % pad with zeros out to the maximum range
                range = [range zeros(1,maxWCSampIdx-length(range))];   
                
                % this is depth
                z = cos(beamAngle'/100*pi/180)*range;  
                
                % this is across-trackd distance
                y = sin(beamAngle'/100*pi/180)*range;
                
                Awc = beamAmp/2;
                
                X = TVGFuncApplied;
                C = TVGOffset;
                clear TS
                TS = Awc + 10*log10((RxBeamWidth/10)*pi/1800*(TxBeamWidth/10)*pi/1800) - X*log10(ones(length(Awc(:,1)),1)*range) + 40*log10(ones(length(Awc(:,1)),1)*range) - C;
                tsBuf1(pingidx,1:length(TS(:,1)),1:length(TS(1,:))) = TS;
                
                steeringangle = beamAngle/100;
                recieveAngle = 1./cos(steeringangle*pi/180);
                
                RxRad = (RxBeamWidth/10)*pi/1800;
                Length = 2*range*sin(RxRad/2);
                
                TxRad = (TxBeamWidth/10)*pi/1800*recieveAngle;
                Width = 2*ones(length(TxRad),length(range)).*range.*(sin(TxRad./2)).';
                BeamArea = Length.*Width;
                
                Tau = 3500/1e5;
                for ii = 1:size(TS,2)
                    
                    Vol_log = 10*log10(BeamArea(:,ii)*Tau*SoundSpeed/2);
                    Sv(:,ii) = TS(:,ii) - Vol_log;
                end
                    
                for aa = 1:NumBeamsHere
                    if DR(aa) ~= 0
                        yBottom(pingidx,aa) = y(aa,DR(aa));
                        zBottom(pingidx,aa) = z(aa,DR(aa));
                    else
                        yBottom(pingidx,aa) = 0;
                        zBottom(pingidx,aa) = 0;
                    end
                end
                
                
                [~,nadir_idx] = min(abs(yBottom(pingidx,:)-35));
                
                if zBottom(nadir_idx) ~= 0
                    nadir_idx = nadir_idx;
                elseif zBottom(nadir_idx + 1) ~= 0
                    nadir_idx = nadir_idx + 1;
                else
                    nadir_idx = nadir_idx - 1;
                end
                
%                 subplot(211)
%                 pcolor(y,z,TS); shading flat; axis equal; set(gca,'ydir','reverse');
%                 set(gca,'fontname','Times'); caxis([-100 0]); colorbar
%                 title(['Ping ' num2str(pingidx) ': TS'])
%                 ylim([0 300])
%                 hold on
%                 plot(yBottom(pingidx,:),zBottom(pingidx,:),'r.')
%                 hold off 
%                 
%                 subplot(212)
%                 pcolor(y,z,Sv); shading flat; axis equal; set(gca,'ydir','reverse');
%                 set(gca,'fontname','Times'); caxis([-100 0]); colorbar
%                 title(['Ping ' num2str(pingidx) ' :Sv'])
%                 hold on
%                 plot(yBottom(pingidx,:),zBottom(pingidx,:),'r.')
%                 ylim([0 300])
%                 hold off
%                 drawnow
% %                 plot(y(bin_idx),z(bin_idx),'rs')
%                 pause(0.05)
                
                cut_Sv  = Sv(:,1:nadir_idx);
                cut_y   = y(:,1:nadir_idx);
                cut_z   = z(:,1:nadir_idx);
                
                bin_edges = -150:10:150;
                
%                 b_x = beam_x(pingidx,:);
%                 b_y = beam_y(pingidx,:);
%                 b_z = beam_z(pingidx,:);
                
                for i = 1:(length(bin_edges)-1)
                   bin_idx = cut_y < bin_edges(i+1) & cut_y > bin_edges(i);
          %          test_mean(i)    = mean(cut_Sv(bin_idx));
          %          test_max(i)     = max(cut_Sv(bin_idx));
                    test1 = movmean(cut_Sv(bin_idx),10);
                    test_mmax(i) = max(test1);
                    y_mid(i) =  mean(cut_y(bin_idx));
                end 
                
                y_mid_bin(pingidx,:) = y_mid;
                z_nadir(pingidx) = zBottom(nadir_idx);
                Sv_mmax(pingidx,:) = test_mmax;
                            
            end
        
        % all other datagrams, skip and move to next    
        else
            % char(73) = I, extra detects datagram
            % char(196) = A, attitude datagram
            
            fread(fid,sz-2,'uint8');
        end
    else
        % the start id was not equal to 2, so either this is not a file
        % containing simrad datagrams or we have made an error reading
        % it...
        disp('start id not equal to 2!')
        pause(.1)
        if feof(fid)
            stopp = 1;
        end
    end
end

fclose(fid)

%HeatMap_0030.nadir_z = z_nadir;
%HeatMap_0030.across_y = y_mid_bin;
%HeatMap_0030.Sv = Sv_mmax;

%save('HeatMap_0030.mat','HeatMap_0030')
   
