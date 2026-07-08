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
startRangeSampNum=keymeta.startRangeSampNum;
xmitSectNum=keymeta.xmitSectNum;

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
title(datalabel)
ylabel('ping number')
xlabel('beam number')

% next tile
nexttile
pcolor(xmitSectNum)
shading flat
cb3=colorbar;
cb3.Label.String='Xmit Sect Num';
title(datalabel)
ylabel('ping number')
xlabel('beam number')

figure(3)
tiledlayout('horizontal')
% next tile
nexttile
plot(TxBeamWidth)

xlabel('ping number')
ylabel('Transmit Beam Width (degrees)')
title(datalabel)

