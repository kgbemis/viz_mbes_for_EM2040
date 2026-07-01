function check_sampinfo(sampinfofile)
%
%

[~,datalabel,~]=fileparts(sampinfofile);

sampinfo=load(sampinfofile);
minSamps=sampinfo.minSamps;
maxSamps=sampinfo.maxSamps;
allSamps=sampinfo.keepSamps;
Npings=length(minSamps);

figure(1)
tiledlayout('horizontal')
% first tile
nexttile
plot(1:Npings,minSamps,'+')
hold on
    plot(1:Npings,maxSamps,'x')
hold off
legend('min','max')
xlabel('ping number')
ylabel('number of samples along range')
title(datalabel)
% next tile
nexttile
pcolor(allSamps)
shading flat
colorbar
title(datalabel)
xlabel('ping number')
ylabel('beam number')

