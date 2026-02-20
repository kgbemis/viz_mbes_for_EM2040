function fill_grid()
%

% input data needed
% variable      units       source
% SoundSpeed    m/s
% SampFreq      ?
% beamAmp       ?
% beamAngle     degrees
% RxBeamWidth   degrees
% TxBeamWidth   degrees
% TVGFuncApplied
% TVGOffset

% import combined metadata and water column data structure


% setup to pull in water column data
YY=zeros(ckNrx(1),maxWCSampIdx,Ndgm);
ZZ=zeros(ckNrx(1),maxWCSampIdx,Ndgm);
SV=zeros(ckNrx(1),maxWCSampIdx,Ndgm);
TT=zeros(ckNrx(1),maxWCSampIdx,Ndgm);
XX=zeros(ckNrx(1),maxWCSampIdx,Ndgm);
xBottom=zeros(Ndgm,ckNrx(1));
dTime=zeros(Ndgm,1);


% loop over datagrams
for idgm=startDgm:endDgm  % just do the first ping for now


    if DatagramNum == NumDatagrams
        pingidx = idgm-startDgm+1;  % changing to use this as the count

        % store both actual datetime and also elapsed time in seconds
        %   although unclear was actually saving this at all!!!
        pingTime(pingidx) = thistime;    % datetime
        elapsedsec=pingTime(pingidx)-pingTime(1);

        % what is this doing
        beamAmp = [beamAmp zeros(Nrx,maxWCSampIdx-length(beamAmp(1,:)))-999];

        % range in meters
        range = (1:(length(beamAmp(1,:))))*SoundSpeed/2/SampFreq;
        % pad with zeros out to the maximum range
        range = [range zeros(1,maxWCSampIdx-length(range))];   

        % this is depth
        z = cos(beamAngle*pi/180)*range;  % my angle in straight degrees

        % this is across-track distance
        y = sin(beamAngle*pi/180)*range;

        % pre calc terms of sonar equation
        Awc = beamAmp/2;
        X = TVGFuncApplied;
        C = TVGOffset;
        clear TS  % why this? is length of TS varying? this is inherited from Liz's code
        RTval=10*log10(RxBeamWidth*pi/180*TxBeamWidth(cenSec)*pi/180);
        TS = Awc + RTval.*ones(size(Awc)) - double(X).*log10(ones(length(Awc(:,1)),1)*range) + 40*log10(ones(length(Awc(:,1)),1)*range) - double(C);
        %tsBuf1(pingidx,1:length(TS(:,1)),1:length(TS(1,:))) = TS;

        steeringangle = beamAngle;
        recieveAngle = 1./cos(steeringangle*pi/180);

        RxRad = RxBeamWidth*pi/180;
        Length = 2*range*sin(RxRad/2);


        TxRad = TxBeamWidth(cenSec)*pi/180*recieveAngle;
        Width = 2*ones(length(TxRad),length(range)).*range.*(sin(TxRad./2));
        BeamArea = Length.*Width;

        Tau = 3500/1e5;
        Sv=zeros(size(TS));
        for ii = 1:size(TS,2)
            Vol_log = 10*log10(BeamArea(:,ii)*Tau*SoundSpeed/2);
            Sv(:,ii) = TS(:,ii) - Vol_log;
        end

    % store gridded data and bottom returns for this (pingidx) datagram
    YY(:,:,pingidx)=y(:,1:maxWCSampIdx);
    ZZ(:,:,pingidx)=z(:,1:maxWCSampIdx);
    SV(:,:,pingidx)=Sv(:,1:maxWCSampIdx);
    TT(:,:,pingidx)=TS(:,1:maxWCSampIdx);
    XX(:,:,pingidx)=cartspeed*elapsedsec*ones(256,maxWCSampIdx);
    xBottom(pingidx,:)=cartspeed*elapsedsec*ones(1,256);
    dTime(pingidx)=elapsedsec;

    end % inside datagram split check

% end of datagram loop
end
