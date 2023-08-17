% read all kmwcd/kmall structures to parse out all important metadata

% used Liz's parsing codes as a start but have pulled in files from the
% CoFFee code to read the new format and better grasp structure of files

%% manual entry of parameters
% sonar parameters
maxWCSampIdx = 600; % setting a fixed watercolumn buffer in meters
RxBeamWidth=1; % setting receive beamwidth arbitrarily (need to find where to read)
cartspeed=0.34; % setting velocity for along track direction

% operations choices

indir='../MBES_mat_files/';


%% open mat file with previously read KMALLdata structures

load(fullfile(indir,'0075_20230414_193631_bothstructs.mat'))


%% datagram types
% this is what we might find
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

% pull actual field lists
fields1=fieldnames(KMALLdata);
fields2=fieldnames(KMALLdata2);
combinedfields=unique([fields1; fields2]);
display(combinedfields)

%% metadata wanted 
% focus on metadata that isn't naturally obtained in the water column 
% processing
% sonar system properties

pruntime=KMALLdata.EMdgmIIP;
fprintf('runtime date = %s \n',...
    datestr(datetime(pruntime.header.time_sec,'ConvertFrom','posixtime')))
fprintf('system & echosounder = %d & %d \n',...
    pruntime.header.systemID, pruntime.header.echoSounderID)

pruntime2=KMALLdata2.EMdgmIIP;
fprintf('runtime date = %s \n',...
    datestr(datetime(pruntime2.header.time_sec,'ConvertFrom','posixtime')))
fprintf('system & echosounder = %d & %d \n',...
    pruntime2.header.systemID, pruntime.header.echoSounderID)

% need to make above check automatically for both entries and report
% whether consistent and the consistent information

return

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
% setting a fixed watercolumn buffer in meters
%   maxWCSampIdx = 600;
% setting receive beamwidth arbitrarily (need to find where to read)
%   RxBeamWidth=1;
% setting velocity for along track direction
%   cartspeed=0.34;

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

    YY(:,:,idgm)=y(:,1:600);
    ZZ(:,:,idgm)=z(:,1:600);
    SV(:,:,idgm)=Sv(:,1:600);
    TT(:,:,idgm)=TS(:,1:600);
    XX(:,:,idgm)=cartspeed*elapsedsec*ones(256,600);
    xBottom(idgm,:)=cartspeed*elapsedsec*ones(1,256);
    dTime(idgm)=elapsedsec;

    end % inside datagram split check

    %pause
end

%% summary visualizations
if do3Dviz
    pingtoplot=Ndgm/2; % 180
    thisbeam=150; % 150
   
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
pcolor(mean(YY(:,:,Ndgm/2:3*Ndgm/4),3),...
    mean(ZZ(:,:,Ndgm/2:3*Ndgm/4),3),mean(TT(:,:,Ndgm/2:3*Ndgm/4),3))
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

