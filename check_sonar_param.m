function check_sonar_param(keymetafile)
%
%

[~,datalabel,~]=fileparts(keymetafile);

keymeta=load(keymetafile);

SoundSpeed=keymeta.SoundSpeed;
SampFreq=keymeta.SampFreq;
Nrx=keymeta.Nrx;

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
%nexttile
%pcolor(allSamps)
%shading flat
%colorbar
%title(datalabel)
%xlabel('ping number')
%ylabel('beam number')
