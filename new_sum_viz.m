function new_sum_viz(beamdatafile,keymetafile)
% new_sum_viz(beamdatafile,keymetafile)
%   loads beam data and visualizes it with minimal processing

load(beamdatafile,'beamAmpData')
% beamAmpData(beamnumber, samplenumber, datagram)
[Nbeams,Nsamps,Ndgm]=size(beamAmpData);

% get metadata
metadata=load(keymetafile);
MSoundSpeed=metadata.SoundSpeed;
MSampFreq=metadata.SampFreq;
MBeamAngle=metadata.beamAngle;

% get data from one ping - for now the midpoint of the run
defaultping=fix(Ndgm/2);
oneping=beamAmpData(:,:,defaultping);
    SoundSpeed=MSoundSpeed(defaultping);
    SampFreq=MSampFreq(defaultping);
    beamAngle=MBeamAngle(defaultping,:);
%whos

% look at value distribution
figure(5)
histogram(oneping(:))

% setup coordinates for a single ping
% range in meters
% this assumes that the starting sample number is always zero
range = (1:Nsamps)'*SoundSpeed/2/SampFreq;
whos
% this is depth
z = range*cos(beamAngle*pi/180);
% this is across-track distance
y = range*sin(beamAngle*pi/180);

whos

% plot single ping
figure(6)
pcolor(y,-z,oneping')
shading flat
colorbar

% NEED to add
% - controlled scaling of plots!
% - labels & titles
% - image saving
% - stacked pings
% - 3D & along track views
% CONSIDER
%   what value does conversion to TS and SV add?
%   does it just make the numbers more "normal"?
%   are they any more or less calibrated?
%   what is relative change across view?