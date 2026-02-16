function fill_grid()

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

% end of datagram loop
end
