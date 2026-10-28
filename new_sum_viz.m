function new_sum_viz(beamdatafile,keymetafile)
% new_sum_viz(beamdatafile,keymetafile)
%   loads beam data and visualizes it with minimal processing

load(beamdatafile,'beamAmpData')
% beamAmpData(beamnumber, samplenumber, datagram)
[Nbeams,Nsamps,Ndgm]=size(beamAmpData);

% not sure if need this
%metadata=load(keymetafile);

defaultping=fix(Ndgm/2);
oneping=beamAmpData(:,:,defaultping);