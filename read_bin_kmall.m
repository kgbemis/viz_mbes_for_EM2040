function out_struct=read_bin_kmall(fname,wcdat)
% using commented code from CFF_read_EMdgmMWC.m to read amplitude and phase
% data

fid=fopen(fname,'r');
phaseFlag=wcdat.rxInfo.phaseFlag;
Nrx=wcdat.rxInfo.numBeams;
% preallocte out_struct
out_struct=struct('sampleAmplitude05dB_p', cell(1, Nrx), 'rxBeamPhase', cell(1, Nrx));
for ibeam=1:Nrx
    % get start of beam data section
    pif=wcdat.beamData_p.sampleDataPositionInFile(ibeam);
    fseek(fid,pif,'bof');
    % get number of samples
    Ns = wcdat.beamData_p.numSampleData(ibeam);

    % next section copied from CFF_read_EMdgmMWC.m
    % ------------------ OPTION 1: ACTUALLY READ DATA ---------------------
    %
    % % Pointer to start of array with Water Column data. Lenght of array =
    % % numSampleData. Sample amplitudes in 0.5 dB resolution. Size of
    % % array is numSampleData * int8_t. Amplitude array is followed by
    % % phase information if phaseFlag >0. Use (numSampleData * int8_t) to
    % % jump to next beam, or to start of phase info for this beam, if
    % % phase flag > 0.
    out_struct(ibeam).sampleAmplitude05dB_p = fread(fid,Ns,'int8');
    %
    switch phaseFlag
    %     % #MWC - Beam sample phase info, specific for each beam and water
    %     % column sample. numBeams * numSampleData = (Nrx * Ns) entries.
         case 1
    %         % Only added to datagram if phaseFlag = 1. Total size of
    %         % phase block is numSampleData * int8_t.
    %
    %         % Rx beam phase in 180/128 degree resolution.
             out_struct(ibeam).rxBeamPhase = fread(fid,Ns,'int8');
    %
         case 2
    %         % Only added to datagram if phaseFlag = 2. Total size of
    %         % phase block is numSampleData * int16_t.
    %
    %         % Rx beam phase in 0.01 degree resolution.
             out_struct(ibeam).rxBeamPhase = fread(fid,Ns,'int16');
    %
    end
    %
    % ------------------ END OF OPTION 1 ----------------------------------
end