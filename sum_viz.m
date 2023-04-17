function sum_viz(filecode,pingtoplot,thisbeam,doprint,printstr,svclims,tsclims)
% visualize MBES data that has already been read from the .kmwcd data files

% load data for summary visualization
datadir='..\MBES_mat_files\';
vizdat=load(fullfile(datadir,[filecode '_viz.mat ']));
    XX=vizdat.XX;
    YY=vizdat.YY;
    ZZ=vizdat.ZZ;
    SV=vizdat.SV; % backscattering volume
    TT=vizdat.TT; % target strength????
    xBottom=vizdat.xBottom; % automatic floor detection results
    yBottom=vizdat.yBottom;
    zBottom=vizdat.zBottom;
load(fullfile(datadir,[filecode '.mat ']),'wcdat');

% extract critical MBES and run information from mat file
Ndgm=length(wcdat); % this is number of WC datagrams
fprintf('number WC datagrams = %d\n',Ndgm)
startdate=datetime(wcdat(1).header.time_sec,'ConvertFrom','posixtime');
startstr=char(startdate);

% setup info for labeling
fileinfo = split(filecode,'_');

% set viz parameters
if isempty(pingtoplot)
    pingtoplot=fix(Ndgm/2); 
end
if isempty(thisbeam)
    thisbeam=150; 
end
fprintf('plotting ping %d and beam %d \n',pingtoplot,thisbeam)
if isempty(svclims)
    svclims=[-140 -40];
end
if isempty(tsclims)
    tsclims=[-140 -40];
end
   
%plot along track profile
f6=figure(6);
pcolor(squeeze(XX(thisbeam,:,:)),squeeze(ZZ(thisbeam,:,:)),squeeze(TT(thisbeam,:,:)))
ylim([0 6])
shading flat
set(gca,'ydir','reverse');
set(gca,'fontname','Times'); 
if sscanf(version('-release'),'%d')<2022
    caxis([-140 -40]); colorbar
else
    %clim([-140 -40]); colorbar    
    clim(tsclims); colorbar    
end
title(['run ' fileinfo{1} ': starts at ' startstr ': TS'])
hold on
 plot(xBottom(:,thisbeam),zBottom(:,thisbeam),'r.')
 plot(xBottom(pingtoplot,thisbeam).*ones(2,1),[0 6],'k')
hold off
xlabel('estimated distance along track (m)')
ylabel('depth below sonar (m)')

f7=figure(7);
pcolor(squeeze(XX(thisbeam,:,:)),squeeze(ZZ(thisbeam,:,:)),squeeze(SV(thisbeam,:,:)))
ylim([0 6])
shading flat
set(gca,'ydir','reverse');
set(gca,'fontname','Times'); 
if sscanf(version('-release'),'%d')<2022
    caxis([-140 -40]); colorbar
else
    %clim([-140 -40]); colorbar    
    clim(svclims); colorbar    
end
title(['run ' fileinfo{1} ': starts at ' startstr ': SV'])
hold on
 plot(xBottom(:,thisbeam),zBottom(:,thisbeam),'r.')
hold off
xlabel('estimated distance along track (m)')
ylabel('depth below sonar (m)')


% plot a summary of pings
f8=figure(8);
pcolor(mean(YY(:,:,fix(Ndgm/2):ceil(3*Ndgm/4)),3),...
    mean(ZZ(:,:,fix(Ndgm/2):ceil(3*Ndgm/4)),3),...
    mean(TT(:,:,fix(Ndgm/2):ceil(3*Ndgm/4)),3))
shading flat
set(gca,'ydir','reverse');
set(gca,'fontname','Times'); 
if sscanf(version('-release'),'%d')<2022
    caxis(tsclims); colorbar
else
    clim(tsclims); colorbar    
end
title(['run ' fileinfo{1} ': starts at ' startstr ': average Ping return: TS'])
hold on
 plot(mean(yBottom,1),mean(zBottom,1),'m.')
 plot(yBottom(pingtoplot,:),zBottom(pingtoplot,:),'r.')
 ylim([0 6])
hold off
xlabel('distance across swath (m)')
ylabel('depth below sonar (m)')

% plot a summary of pings
f9=figure(9);
pcolor(mean(YY(:,:,fix(Ndgm/2):ceil(3*Ndgm/4)),3),...
    mean(ZZ(:,:,fix(Ndgm/2):ceil(3*Ndgm/4)),3),...
    mean(SV(:,:,fix(Ndgm/2):ceil(3*Ndgm/4)),3))
shading flat
set(gca,'ydir','reverse');
set(gca,'fontname','Times'); 
if sscanf(version('-release'),'%d')<2022
    caxis(svclims); colorbar
else
    clim(svclims); colorbar    
end
title(['run ' fileinfo{1} ': starts at ' startstr ': average Ping return: SV'])
hold on
 plot(mean(yBottom,1),mean(zBottom,1),'m.')
 plot(yBottom(pingtoplot,:),zBottom(pingtoplot,:),'r.')
 ylim([0 6])
hold off
xlabel('distance across swath (m)')
ylabel('depth below sonar (m)')


% plot a specific ping
f10=figure(10);
pcolor(YY(:,:,pingtoplot),ZZ(:,:,pingtoplot),TT(:,:,pingtoplot))
ylim([0 6])
shading flat
set(gca,'ydir','reverse');
set(gca,'fontname','Times'); 
if sscanf(version('-release'),'%d')<2022
    caxis([-140 -40]); colorbar
else
    clim([-140 -40]); colorbar    
end
title(['run ' fileinfo{1} ': starts at ' startstr ': Ping ' num2str(pingtoplot) ': TS'])
hold on
 plot(yBottom(pingtoplot,:),zBottom(pingtoplot,:),'r.')
 plot(YY(thisbeam,:,pingtoplot),ZZ(thisbeam,:,pingtoplot),'k')
hold off
xlabel('distance across swath (m)')
ylabel('depth below sonar (m)')

% plot a specific ping
f11=figure(11);
pcolor(YY(:,:,pingtoplot),ZZ(:,:,pingtoplot),SV(:,:,pingtoplot))
ylim([0 6])
shading flat
set(gca,'ydir','reverse');
set(gca,'fontname','Times'); 
if sscanf(version('-release'),'%d')<2022
    caxis([-140 -40]); colorbar
else
    clim([-140 -40]); colorbar    
end
title(['run ' fileinfo{1} ': starts at ' startstr ': Ping ' num2str(pingtoplot) ': SV'])
hold on
 plot(yBottom(pingtoplot,:),zBottom(pingtoplot,:),'r.')
 %plot(YY(thisbeam,:,pingtoplot),ZZ(thisbeam,:,pingtoplot),'k')
hold off
xlabel('distance across swath (m)')
ylabel('depth below sonar (m)')


%figure(12)
%p = patch(isosurface(XX,YY,ZZ,SV,-80));
%isonormals(XX,YY,ZZ,SV,p)
%p.FaceColor = 'red';
%p.EdgeColor = 'none';
%daspect([1 1 1])
%view(3); 
%axis tight
%camlight 
%lighting 

%%
% print figures
if doprint
    f6name=[filecode '_alongtrack_TS_' printstr '.png'];
    print(f6,f6name,'-r300','-dpng')
    f7name=[filecode '_alongtrack_SV' printstr '.png'];
    print(f7,f7name,'-r300','-dpng')
    f8name=[filecode '_averageping_TS' printstr '.png'];
    print(f8,f8name,'-r300','-dpng')
    f9name=[filecode '_averageping_SV' printstr '.png'];
    print(f9,f9name,'-r300','-dpng')
    f10name=[filecode '_singleping_TS' printstr '.png'];
    print(f10,f10name,'-r300','-dpng')
    f11name=[filecode '_singleping_SV' printstr '.png'];
    print(f11,f11name,'-r300','-dpng')
    %f12name=[filecode '_3D_TS' printstr '.png'];
    %print(f12,f12name,'-r300','-dpng')
end % doprint