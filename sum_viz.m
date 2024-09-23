function sum_viz(filecode,pingtoplot,thisbeam,doprint,printstr,svclims,tsclims)
% visualize MBES data that has already been read from the .kmwcd data files
%
%   filecode = <kmwcd or kmall base filename>_kmwcd
%


%   set ringing level
indRing=1100;

% set minimun reasonable bottom level
minBot=10; % this depends on altitude of sonar above bottom

% load data for summary visualization
datadir='..\MBES_mat_files\';
vizdat=load(fullfile(datadir,[filecode '_viz.mat ']));
    XX=vizdat.XX;
    YY=vizdat.YY;
    ZZ=vizdat.ZZ;
    SV=vizdat.SV; % backscattering volume
    TT=vizdat.TT; % target strength????
    xBottom=vizdat.xBottom; % automatic floor detection results
    [lenBot,nbotbeams]=size(xBottom);
    yBottom=vizdat.yBottom(1:lenBot,:);
    zBottom=vizdat.zBottom(1:lenBot,:);
    [nBeams,mSamples,pDgms]=size(XX);
    fprintf('viz mat file has data for %d beams by %d samples for %d pings \n',...
        nBeams,mSamples,pDgms)
    fprintf('bottom picks done for %d pings and %d beams\n',lenBot,nbotbeams)
load(fullfile(datadir,[filecode '.mat ']),'wcdat');

% extract critical MBES and run information from mat file
Ndgm=length(wcdat); % this is number of WC datagrams
fprintf('total number WC datagrams in raw data file = %d \n',Ndgm)
fprintf('available datagrams in this file = %d \n',pDgms)
startDgm=vizdat.startDgm;
endDgm=vizdat.endDgm;
fprintf('using pings %d to %d \n',startDgm,endDgm)
if startDgm>1
    % adjust ping to plot
    pingtoplot=pingtoplot-startDgm+1;
end

startdate=datetime(wcdat(startDgm).header.time_sec,'ConvertFrom','posixtime');
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

% determine issues with bottom picks for later use
figure(17)
subplot(131)
    pcolor(zBottom)
    shading flat
    colorbar
    title('map of bottom picks - z value')
subplot(132)
pcolor(yBottom)
    shading flat
    colorbar
    title('map of bottom picks - y value')
subplot(133)
    pcolor(xBottom)
    shading flat
    colorbar
    title('map of bottom picks - x value')
% need to fill in all the missing points where there was no bottom pick
% try fillmissing2 function to see if this works
% maps look like linear interpolation will work fine for both z & y
% x values do not need interpolation

% locate missing locations 
%      z==0 will work because values unset where no bottom pick was valid
%      and otherwise z=0 is at the sonar so not a valid value
BotMissing=zBottom==0;

% fill them
[zBotFilled,FillFlags]=fillmissing2(zBottom,"linear",'MissingLocations',BotMissing);
[yBotFilled,~]=fillmissing2(yBottom,"linear",'MissingLocations',BotMissing);

% check new values
figure(18)
subplot(141)
    pcolor(zBotFilled)
    shading flat
    colorbar
    title('map of bottom picks - z value - after missing filled')
subplot(142)
    pcolor(yBotFilled)
    shading flat
    colorbar
    title('map of bottom picks - y value - after missing filled')
subplot(143)
    pcolor(FillFlags)
    shading flat
    colorbar
    title('map of filled points')
subplot(144)
    pcolor(zBotFilled-zBottom)
    shading flat
    colorbar
    title('z: original-filled')


% set the bottom picks to the filled versions
zBottom=zBotFilled;
yBottom=yBotFilled;
   
%plot along track profile
meanBottom=mean(zBottom(:,thisbeam));
botbuf=10;
f6=figure(6);
pcolor(squeeze(XX(thisbeam,:,:)),squeeze(ZZ(thisbeam,:,:)),squeeze(TT(thisbeam,:,:)))
ylim([0 meanBottom+botbuf])
shading flat
set(gca,'ydir','reverse');
set(gca,'fontname','Times'); 
if sscanf(version('-release'),'%d')<2022
    caxis([-140 -40]); colorbar %#ok<CAXIS>
else
    %clim([-140 -40]); colorbar    
    clim(tsclims); colorbar    
end
title(['run ' fileinfo{1} ': starts at ' startstr ': TS'])
hold on
 plot(xBottom(:,thisbeam),zBottom(:,thisbeam),'r.')
 plot(xBottom(pingtoplot,thisbeam).*ones(2,1),[0 meanBottom+botbuf],'k')
hold off
xlabel('estimated distance along track (m)')
ylabel('depth below sonar (m)')

f7=figure(7);
pcolor(squeeze(XX(thisbeam,:,:)),squeeze(ZZ(thisbeam,:,:)),squeeze(SV(thisbeam,:,:)))
ylim([0 meanBottom+botbuf])
shading flat
set(gca,'ydir','reverse');
set(gca,'fontname','Times'); 
if sscanf(version('-release'),'%d')<2022
    caxis([-140 -40]); colorbar %#ok<CAXIS>
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
pcolor(mean(YY(:,:,fix(pDgms/2):ceil(3*pDgms/4)),3),...
    mean(ZZ(:,:,fix(pDgms/2):ceil(3*pDgms/4)),3),...
    mean(TT(:,:,fix(pDgms/2):ceil(3*pDgms/4)),3))
shading flat
set(gca,'ydir','reverse');
set(gca,'fontname','Times'); 
if sscanf(version('-release'),'%d')<2022
    caxis(tsclims); colorbar %#ok<CAXIS>
else
    clim(tsclims); colorbar    
end
title(['run ' fileinfo{1} ': starts at ' startstr ': average Ping return: TS'])
hold on
 plot(mean(yBottom,1),mean(zBottom,1),'m.')
 plot(yBottom(pingtoplot,:),zBottom(pingtoplot,:),'r.')
 plot(zeros(indRing,1),ZZ(fix(nBeams/2),1:indRing,1),'w')
 ylim([0 meanBottom+botbuf])
hold off
xlabel('distance across swath (m)')
ylabel('depth below sonar (m)')

% plot a summary of pings
f9=figure(9);
pcolor(mean(YY(:,:,fix(pDgms/2):ceil(3*pDgms/4)),3),...
    mean(ZZ(:,:,fix(pDgms/2):ceil(3*pDgms/4)),3),...
    mean(SV(:,:,fix(pDgms/2):ceil(3*pDgms/4)),3))
shading flat
set(gca,'ydir','reverse');
set(gca,'fontname','Times'); 
if sscanf(version('-release'),'%d')<2022
    caxis(svclims); colorbar %#ok<CAXIS>
else
    clim(svclims); colorbar    
end
title(['run ' fileinfo{1} ': starts at ' startstr ': average Ping return: SV'])
hold on
 plot(mean(yBottom,1),mean(zBottom,1),'m.')
 plot(yBottom(pingtoplot,:),zBottom(pingtoplot,:),'r.')
 ylim([0 meanBottom+botbuf])
hold off
xlabel('distance across swath (m)')
ylabel('depth below sonar (m)')


% plot a specific ping
f10=figure(10);
pcolor(YY(:,:,pingtoplot),ZZ(:,:,pingtoplot),TT(:,:,pingtoplot))
ylim([0 meanBottom+botbuf])
shading flat
set(gca,'ydir','reverse');
set(gca,'fontname','Times'); 
if sscanf(version('-release'),'%d')<2022
    caxis([-140 -40]); colorbar %#ok<CAXIS>
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
ylim([0 meanBottom+botbuf])
shading flat
set(gca,'ydir','reverse');
set(gca,'fontname','Times'); 
if sscanf(version('-release'),'%d')<2022
    caxis([-140 -40]); colorbar %#ok<CAXIS>
else
    clim([-140 -40]); colorbar    
end
title(['run ' fileinfo{1} ': starts at ' startstr ': Ping ' num2str(pingtoplot) ': SV'])
hold on
 plot(yBottom(pingtoplot,:),zBottom(pingtoplot,:),'r.')
 plot(YY(thisbeam,:,pingtoplot),ZZ(thisbeam,:,pingtoplot),'k')
hold off
xlabel('distance across swath (m)')
ylabel('depth below sonar (m)')

% add a beam line and look at average values
sidebeam=50;
hold on
plot(YY(sidebeam,:,pingtoplot),ZZ(sidebeam,:,pingtoplot),'k')
hold off

figure(12)
Nsamp=length(ZZ(sidebeam,:,pingtoplot));
bot=Nsamp-100;
fprintf('number Z samps = %d should stop at %d\n',Nsamp,bot)
subplot(121)
    plot(SV(sidebeam,1:bot,pingtoplot),-ZZ(sidebeam,1:bot,pingtoplot),'k')
    ylabel('depth (m)')
    xlabel('scattering volume (dB)')
subplot(122)
    avgZdep=zeros(bot,1);
    for i=1:bot
        avgZdep(i)=mean(SV(sidebeam,1:i,pingtoplot));
    end
    plot(avgZdep,-ZZ(sidebeam,1:bot,pingtoplot),'k')
    ylabel('depth (m)')
    xlabel('cummulative average (dB)')

% compute and display average water column backscatter as map
xloc=zeros(pDgms,nBeams);
yloc=zeros(pDgms,nBeams);
zloc=zeros(pDgms,nBeams);
avgTS=zeros(pDgms,nBeams);
avgTSabovering=zeros(pDgms,nBeams);
maxTS=zeros(pDgms,nBeams);
maxTSabovering=zeros(pDgms,nBeams);
botInd=zeros(pDgms,nBeams);
copyTT=TT;
for iping=1:pDgms
    for ibeam=1:nBeams
        % get bottom pick
        % starting by assuming xBottom & zBottom have pDgmx pings and
        %   nBeams beams ==> this works!
        xloc(iping,ibeam)=xBottom(iping,ibeam);
        yloc(iping,ibeam)=yBottom(iping,ibeam);
        zloc(iping,ibeam)=zBottom(iping,ibeam);
        %fprintf('bottom pick for ping %d beam %d is %4.2f,%4.2f,%4.2f \n',...
        %    iping,ibeam,xloc(iping,ibeam),yloc(iping,ibeam),zloc(iping,ibeam))
        % check for valid bottom pick
        if zloc(iping,ibeam)<minBot
            fprintf('invalid bottom pick detected \n')
            % get zloc for next or previous beams
            switch ibeam
                case 1
                    alt_zloc=zloc(iping,ibeam+1);
                    alt_yloc=yloc(iping,ibeam+1);
                    % ARGH - need to reset yloc as well and this isn't
                    % correct as need value corresponding to beam
                case nBeams
                    alt_zloc=zloc(iping,ibeam-1);
                    alt_yloc=yloc(iping,ibeam-1);
                otherwise
                    alt_zloc=max(zloc(iping,ibeam-1:ibeam+1));
                    alt_yloc=max(zloc(iping,ibeam-1:ibeam+1));
            end
            % check new zloc and if good, set value
            if alt_zloc<minBot
                fprintf('still no valid bottom pick\n')
                %fprintf('alt bottom pick for ping %d beam %d is %4.2f,%4.2f,%4.2f \n',...
                %    iping,ibeam,xloc(iping,ibeam),yloc(iping,ibeam),zloc(iping,ibeam))
            else
                zloc(iping,ibeam)=alt_zloc;
                yloc(iping,ibeam)=alt_yloc;
                %fprintf('resetting bottom pick for ping %d beam %d is %4.2f,%4.2f,%4.2f \n',...
                %    iping,ibeam,xloc(iping,ibeam),yloc(iping,ibeam),zloc(iping,ibeam))
            end
        end

        % find index where zBottom == ZZ for this ping and beam
        thisZZ=squeeze(ZZ(ibeam,:,iping));
        indBot=find((thisZZ-zloc(iping,ibeam))<0,1,'last');
        botInd(iping,ibeam)=indBot;
        % average water column backscatter above bottom
        avgTS(iping,ibeam)=mean(TT(ibeam,1:indBot,iping));
        maxTS(iping,ibeam)=max(TT(ibeam,1:indBot-10,iping));
        copyTT(ibeam,indBot-40:end,iping)=NaN;
        % also average water column backscatter above first ringing ping
        avgTSabovering(iping,ibeam)=mean(TT(ibeam,1:indRing,iping));
        maxTSabovering(iping,ibeam)=max(TT(ibeam,1:indRing,iping));

    end
end


figure(13)
pcolor(xloc,yloc,avgTS)
shading flat
set(gca,'fontname','Times'); 
if sscanf(version('-release'),'%d')<2022
    caxis([-120 -60]); colorbar %#ok<CAXIS>
else
    clim([-120 -60]); colorbar    
    %clim(tsclims); colorbar    
end
title(['run ' fileinfo{1} ': start ' startstr ': average TS above bottom'])
xlabel('estimated distance along track (m)')
ylabel('cross distance (m)')


figure(14)
pcolor(xloc,yloc,avgTSabovering)
shading flat
set(gca,'fontname','Times'); 
if sscanf(version('-release'),'%d')<2022
    caxis([-120 -60]); colorbar %#ok<CAXIS>
else
    clim([-120 -60]); colorbar    
    %clim(tsclims); colorbar    
end
title(['run ' fileinfo{1} ': start ' startstr ': average TS above arc'])
xlabel('estimated distance along track (m)')
ylabel('cross distance (m)')


figure(15)
pcolor(xloc,yloc,maxTS)
shading flat
set(gca,'fontname','Times'); 
if sscanf(version('-release'),'%d')<2022
    caxis([-140 -40]); colorbar %#ok<CAXIS>
else
    %clim([-140 -40]); colorbar    
    clim(tsclims); colorbar    
end
title(['run ' fileinfo{1} ': start ' startstr ': max TS above bottom'])
xlabel('estimated distance along track (m)')
ylabel('cross distance (m)')


figure(16)
pcolor(xloc,yloc,maxTSabovering)
shading flat
set(gca,'fontname','Times'); 
if sscanf(version('-release'),'%d')<2022
    caxis([-140 -40]); colorbar %#ok<CAXIS>
else
    %clim([-140 -40]); colorbar    
    clim(tsclims); colorbar    
end
title(['run ' fileinfo{1} ': start ' startstr ': max TS above arc'])
xlabel('estimated distance along track (m)')
ylabel('cross distance (m)')

% figures 17 and 18 up top
figure(19)
surf(xloc,yloc,-zloc,maxTSabovering)
shading flat
colorbar
grid on
view(120,55)
xlabel('x dist')
ylabel('y dist')
zlabel('z dist')
title(['run ' fileinfo{1} ': start ' startstr ': max TS above arc'])

figure(20)
surf(xloc,yloc,-zloc,avgTSabovering)
shading flat
colorbar
grid on
view(120,55)
xlabel('x dist')
ylabel('y dist')
zlabel('z dist')
title(['run ' fileinfo{1} ': start ' startstr ': average TS above arc'])

figure(21)
surf(xloc,yloc,-zloc,maxTSabovering)
shading flat
hold on
dBlevel=-50;
hival=copyTT>dBlevel;
scatter3(XX(hival),YY(hival),-ZZ(hival),ZZ(hival),TT(hival))
grid on
hold off
view(120,55)
xlabel('x dist')
ylabel('y dist')
zlabel('z dist')
title(['run ' fileinfo{1} ': start ' startstr ': TS > ' num2str(dBlevel)])

%figure(xx)
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