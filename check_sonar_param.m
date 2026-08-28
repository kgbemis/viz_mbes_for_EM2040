function check_sonar_param(keymetafile)
%
%

[~,datalabel,~]=fileparts(keymetafile);

keymeta=load(keymetafile);

SoundSpeed=keymeta.SoundSpeed;
SampFreq=keymeta.SampFreq;
Nrx=keymeta.Nrx;
TVGFuncApplied=keymeta.TVGFuncApplied;
TVGOffset=keymeta.TVGOffset;
TxBeamWidth=keymeta.TxBeamWidth;
RxBeamWidth=keymeta.RxBeamWidth;
cenSec=keymeta.cenSec;
cenFreq=keymeta.cenFreq;
humFreq=keymeta.humFreq;
startRangeSampNum=keymeta.startRangeSampNum;
xmitSectNum=keymeta.xmitSectNum;
beamAngle=keymeta.beamAngle;

Npings=length(SoundSpeed);

figure(1)
tiledlayout('vertical')
% first tile
nexttile
plot(1:Npings,SoundSpeed,'+')
xlabel('ping number')
ylabel('sound speed (m/s)')
title(datalabel)
% next tile
nexttile
plot(1:Npings,SampFreq,'x')
xlabel('ping number')
ylabel('sampling frequency (Hz)')
title(datalabel)
% next tile
nexttile
plot(1:Npings,Nrx,'x')
xlabel('ping number')
ylabel('number of beams')
title(datalabel)
% next tile
nexttile
plot(1:Npings,TVGFuncApplied,'x')
hold on
plot(1:Npings,TVGOffset,'+')
hold off
xlabel('ping number')
ylabel('TVG')
title(datalabel)
legend('TVGFuncApplied','TVGOffset (dB)')


figure(2)
tiledlayout('horizontal')
% next tile
nexttile
pcolor(TxBeamWidth)
shading flat
ylabel('ping number')
xlabel('sector')
cb=colorbar;
cb.Label.String='Transmit Beam Width (degrees)';
title(datalabel)

% next tile
nexttile
pcolor(startRangeSampNum)
shading flat
cb2=colorbar;
cb2.Label.String='Starting Sample Number for Range';
not0=sum(find(startRangeSampNum(:)~=0));
text(10,75,['not 0 = ' num2str(not0)])
title(datalabel)
ylabel('ping number')
xlabel('beam number')

% next tile
nexttile
pcolor(xmitSectNum)
shading flat
cb3=colorbar;
cb3.Label.String='Xmit Sect Num';
mxn=mean(xmitSectNum(:));
sxn=std(xmitSectNum(:));
text(10,25,[num2str(mxn) '+/-' num2str(sxn)],'Color','w')
not1=sum(find(xmitSectNum(:)~=1));
is1=sum(find(xmitSectNum(:)==1));
text(10,50,['not 1 = ' num2str(100*not1/is1) '%'],'Color','w')
text(10,75,['not 1 = ' num2str(not1)],'Color','w')
title(datalabel)
ylabel('ping number')
xlabel('beam number')

figure(3)
tiledlayout('horizontal')
% next tile - TX beam width
nexttile
    plot(TxBeamWidth)
    ylim([0 1.5])
    xlabel('ping number')
    ylabel('Transmit Beam Width (degrees)')
    title(datalabel)
    [ntx,mtx]=size(TxBeamWidth);
    size(TxBeamWidth)
    nsec=min([ntx mtx]);
    for i=1:nsec
        text(1,1.0+i/20,['sec ' num2str(i) ' bw = ' num2str(mean(TxBeamWidth(:,i)))])
    end
% next tile - RX beamwidth
nexttile
    plot(RxBeamWidth)
    ylim([0 1.5])
    xlabel('ping number')
    ylabel('Recieve Beam Width (degrees)')
    title(datalabel)
% next tile - center frequency
nexttile
    plot(humFreq,'--')
    hold on
        plot(cenFreq)
        plot(cenFreq(:,fix(mean(cenSec))),':k')
        for i=1:nsec
            text(1,75+i*30,['sec ' num2str(i) ' Cfreq = ' num2str(mean(cenFreq(:,i)))])
        end    
        text(1,25,['nominal Cfreq = ' num2str(mean(humFreq))])
    hold off
    ylim([0 800])
    xlabel('ping number')
    ylabel('Central Frequency (Hz)')
    title(datalabel)
    legend('nominal','sectors')

figure(4)
colormap(hsv(25))
tiledlayout('horizontal')
%next tile
nexttile
pcolor(beamAngle)
shading flat
cb4=colorbar;
cb4.Label.String='beam angle (degrees)';
ylabel('ping number')
xlabel('beam number')


