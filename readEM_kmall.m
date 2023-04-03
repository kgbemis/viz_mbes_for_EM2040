% play with reading the kmall files

% used Liz's parsing codes as a start but have pulled in files from the
% CoFFee code to read the new format and better grasp structure of files

%% open an EM data file 

filelocation = './';
%filename='0009_20200918_094230.kmall';
filename='0010_20200918_095915.kmall';

fname = fullfile(filelocation,filename);

% open selected file
fid = fopen(fname,'r');

%% test reading the data file
% this part of code copied from CoFFee
%   Authors: Alex Schimel (NGU, alexandre.schimel@ngu.no) and Yoann
%   Ladroit (NIWA, yoann.ladroit@niwa.co.nz)
%   2017-2021; Last revision: 20-08-2021

    fprintf('Initial header info in file:\n')
% Datagram length in bytes. The length field at the start (4 bytes) and end
% of the datagram (4 bytes) are included in the length count.
out_struct.numBytesDgm = fread(fid,1,'uint32');
    fprintf('size of next datagram: %d\n',out_struct.numBytesDgm)

% Multi beam datagram type definition, e.g. #AAA
out_struct.dgmType = fscanf(fid,'%c',4);
    fprintf('datagram type: %s\n',out_struct.dgmType)

% Datagram version.
out_struct.dgmVersion = fread(fid,1,'uint8');
    fprintf('datagram version: %d\n',out_struct.dgmVersion)

% System ID. Parameter used for separating datagrams from different
% echosounders if more than one system is connected to SIS/K-Controller.
out_struct.systemID = fread(fid,1,'uint8');
    fprintf('system ID: %d\n',out_struct.systemID)

% Echo sounder identity, e.g. 124, 304, 712, 2040, 2045 (EM 2040C)
out_struct.echoSounderID = fread(fid,1,'uint16');
    fprintf('echo sounder identity: %d\n',out_struct.echoSounderID)

% UTC time in seconds. Epoch 1970-01-01. time_nanosec part to be added for
% more exact time.
out_struct.time_sec = fread(fid,1,'uint32');
    fprintf('time seconds: %d\n',out_struct.dgmVersion)

% Nano seconds remainder. time_nanosec part to be added to time_sec for
% more exact time.
out_struct.time_nanosec = fread(fid,1,'uint32');
    fprintf('time nanoseconds: %d\n',out_struct.time_nanosec)

    fprintf('\n')
% close file as next steps will reopen
fclose(fid);

%% Call CFF function to get datagram structure in kmall file
KMALLfileinfo = CFF_kmall_file_info(filename);
%   * |KMALLfileinfo|: structure about datagrams in KMALLfilename, with fields:
%     * |file_name|: input file name
%     * |fileSize|: file size in bytes
%     * |dgm_num|: number of datagram in file
%     * |dgm_type_code|: datagram type as string, e.g. '#IIP' (Kongsberg .kmall format)
%     * |dgm_type_text|: datagram type description (Kongsberg .kmall format) 
%     * |dgm_type_version|: version for this type of datagram, as int (Kongsberg .kmall format)
%     * |dgm_counter|: counter for this type and version of datagram in the
%           file. There should not be multiple versions of a same type in
%           a same file, but we never know...
%     * |dgm_start_pif|: position of beginning of datagram in file 
%     * |dgm_size|: datagram size in bytes
%     * |dgm_sys_ID|: System ID. Parameter used for separating datagrams from different echosounders.
%     * |dgm_EM_ID|: Echo sounder identity, e.g. 124, 304, 712, 2040, 2045 (EM 2040C)
%     * |sync_counter|: number of bytes found between this datagram and the
%           previous one (any number different than zero indicates a sync error)
%     * |date_time|: datagram date in datetime format
%     * |parsed|: flag for whether the datagram has been parsed. Initiated
%               at 0 at this stage. To be later turned to 1 for parsing.
% structure also has
%     * |list_dgm_type|: this has the unique dgm_type_codes
%     * |list_dgm_counter|: not sure if counts of instances?

%% explore what is in the kmall file
figure(1)
plot(KMALLfileinfo.dgm_num,'.')
fprintf('min dgm_num= %d; max dgm_num=%d\n',...
    min(KMALLfileinfo.dgm_num),max(KMALLfileinfo.dgm_num))
fprintf('min sync_counter= %d; max sync_counter=%d\n',...
    min(KMALLfileinfo.sync_counter),max(KMALLfileinfo.sync_counter))
fprintf('min parsed= %d; max parsed=%d\n',...
    min(KMALLfileinfo.parsed),max(KMALLfileinfo.parsed))

figure(2)
b1=barh(KMALLfileinfo.list_dgm_counter,...
    'FaceColor',[0.7 0.9 0.7],'EdgeColor',[0.7 0.9 0.7]);
yticklabels(KMALLfileinfo.list_dgm_type)
ylabel('Datagram Types')
xlabel('Counts')
datalabels=string(b1(1).YData);
infolabels={'IIP=Installation parameters and sensor setup',...
    'IOP= Runtime parameters as chosen by operator',...
    'SVP= sensor data from sound velocity profile or CTD',...
    'SKM= sensor KM binar sensor format',...
    'MWC= Multibeam water column datagram', ...
    'SVT= sensor data for sound velocity at transducer',...
    'CHE= Compatibility data for heave',...
    'SPO= sensor data for position',...
    'CPO= Compatibility data for position',...
    'MRZ= Multibeam raw range and depth datagram',...
    'SCL= sensor data from clock'};
text(200*ones(size(yticklabels)),b1(1).XEndPoints,datalabels,'VerticalAlignment','middle')
text(2000*ones(size(yticklabels)),b1(1).XEndPoints,infolabels,'VerticalAlignment','middle')

%% now we want to read some datagrams and figure out what to do with the data!
% this should read all the data
KMALLdata = CFF_read_kmall(filename);

% now have a set of stuctures in KMALLdata that contain the datagrams for
% each type plus an info strucure

wcdat=KMALLdata.EMdgmMWC;

return
% every thing beyond here is Liz's - mostly using as guide to which
% datagrams are needed and other steps

%% start reading data file

%tic
pingidx = 0;
rangeSB = [];
BufSize = 1000;
tsBuf1 = zeros(BufSize,256,2500);
stopp = 0;
bufCt1 = 0;

while ~stopp

    sz = fread(fid,1,'uint32');     % the size of the next datagram. note: the size does not include this field
        fprintf('size of next datagram: %d\n',sz)
    startid = fread(fid,1,'uint8'); % the start id of a datagram should be 2
        fprintf('start id: %d\n',startid)
    
    if startid == 2
       
        datatype = fread(fid,1,'char'); % this is the type of datagram that we are reading
         fprintf('data type id is %d',datatype)
       
        % char(83) = S, seabead image data
        if datatype == 83
            [DateVal,TimeVal,PingCount,MeanAbsorp,PulseLength,NormIncRange,TVGstartRange,...
                TVGstopRange,NormincBSN,OblqincBSO,TxBeamWidth,TVGcrossoverAng,NumBeamsHere,...
                beamIndexNum,sortDir,NumSampsPerBeam,BeamCenterSampNum,SimradSnippets,BeamDetAmp] = readEM_S(fid);
        pause
        % char(82) = R, runtime parameters data    
        elseif datatype == 82
            [ModelNum,DateVal,TimeVal,PingCount,SystemSerNum,Mode,MinDepth,MaxDepth,AbsCoef,TxPulseLength,...
                TxBeamWidth, TxPower,RxBeamWidth,RxBandWidth,RxFixedGain, TVGXover,MaxPortCoverage, MaxStbdCoverage,...
                MaxPortSwath,BeamSpacing,MaxStbdSwath] = readEM_R(fid);
            
        elseif datatype == 78
            [DateVal,TimeVal,PingCount,SSPD,SampFreq,TiltAng,FocRange,SignalLength,SectorXmitDelay,CentFreq,MeanAbsorpt,SigWavID,XmitSecNum,BW,BeamAngle,XmitSecNumByBeam,TwoWayTT,DetInfo,BS] = readEM_N(fid);
            %pause
        
        % char(107) = k, watercolumn data    
        elseif datatype == 107
            [DateVal,TimeVal,PingCount,NumDatagrams,DatagramNum,NumXmitSect,...
                NumRecBeams,NumBeamsHereTemp, SoundSpeed, SampFreq, TxTimeHeave, TVGFuncApplied, TVGOffset,...
                centFreqTemp,tiltAngle, xmitSectNumTemp,beamAngleTemp,startRangeSampNumTemp,numSampsTemp,xmitSectNumRecTemp,...
                beamNumTemp,beamAmpTemp,DRtemp] = readEM_k(fid);
          
            % convert datagram time to matlab datenumber
            YY = floor(DateVal/10000);
            MM = floor( (DateVal-YY*10000)/100);
            DD = floor( (DateVal-YY*10000-MM*100));
            dn = datenum(YY,MM,DD)+TimeVal/1000/3600/24;
            
            % there may be more than one water column datagram due to the
            % size limit on a single datagram, and we want to stitch all of
            % the data for a single ping together.  The offset calculated
            % below allows the data to be associated with its proper beam
            % number
            NumBeamsHere(DatagramNum) = NumBeamsHereTemp;
            if DatagramNum > 1
                offsetidx = sum(NumBeamsHere(1:(DatagramNum-1)));
            else
                offsetidx = 0;
            end
            
            indx = 1:length(beamAngleTemp);
            beamAngle(offsetidx + indx) = beamAngleTemp;
            startRangeSampNum(offsetidx + indx) = startRangeSampNumTemp;
            numSamps(offsetidx + indx) = numSampsTemp;
            xmitSectNum(offsetidx + indx) = xmitSectNumRecTemp;
            beamNum(offsetidx + indx) = beamNumTemp;
            beamAmp(offsetidx + indx,1:length(beamAmpTemp(1,:))) = beamAmpTemp;
            DR(offsetidx + indx) = DRtemp;
            
            % if we have read all of the data for a single ping (plot this
            % data)
            if DatagramNum == NumDatagrams
                %disp( [num2str(floor(TimeVal/1000/3600)) ':' num2str(floor(rem(TimeVal/1000/3600,1)*60),'%2d') ':' num2str(rem(rem(TimeVal/1000/3600,1)*60,1)*60,'%2.2f')])

                pingidx = pingidx + 1

                pingTime(pingidx) = TimeVal/1000;    % time in seconds since midnight
                
                beamAmp = [beamAmp zeros(NumRecBeams,maxWCSampIdx-length(beamAmp(1,:)))-999];

                %disp({'SoundSpeed=' SoundSpeed});
                %disp({'SampFreq=' SampFreq});
                
                % range in meters
                range =  (1:(length(beamAmp(1,:))))*SoundSpeed/10/2/(SampFreq/100);
                
                % pad with zeros out to the maximum range
                range = [range zeros(1,maxWCSampIdx-length(range))];   
                
                % this is depth
                z = cos(beamAngle'/100*pi/180)*range;  
                
                % this is across-trackd distance
                y = sin(beamAngle'/100*pi/180)*range;
                
                Awc = beamAmp/2;
                
                X = TVGFuncApplied;
                C = TVGOffset;
                clear TS
                TS = Awc + 10*log10((RxBeamWidth/10)*pi/1800*(TxBeamWidth/10)*pi/1800) - X*log10(ones(length(Awc(:,1)),1)*range) + 40*log10(ones(length(Awc(:,1)),1)*range) - C;
                tsBuf1(pingidx,1:length(TS(:,1)),1:length(TS(1,:))) = TS;
                
                steeringangle = beamAngle/100;
                recieveAngle = 1./cos(steeringangle*pi/180);
                
                RxRad = (RxBeamWidth/10)*pi/1800;
                Length = 2*range*sin(RxRad/2);
                
                TxRad = (TxBeamWidth/10)*pi/1800*recieveAngle;
                Width = 2*ones(length(TxRad),length(range)).*range.*(sin(TxRad./2)).';
                BeamArea = Length.*Width;
                
                Tau = 3500/1e5;
                for ii = 1:size(TS,2)
                    
                    Vol_log = 10*log10(BeamArea(:,ii)*Tau*SoundSpeed/2);
                    Sv(:,ii) = TS(:,ii) - Vol_log;
                end
                    
                for aa = 1:NumBeamsHere
                    if DR(aa) ~= 0
                        yBottom(pingidx,aa) = y(aa,DR(aa));
                        zBottom(pingidx,aa) = z(aa,DR(aa));
                    else
                        yBottom(pingidx,aa) = 0;
                        zBottom(pingidx,aa) = 0;
                    end
                end
                
                
                [~,nadir_idx] = min(abs(yBottom(pingidx,:)-35));
                
                if zBottom(nadir_idx) ~= 0
                    nadir_idx = nadir_idx;
                elseif zBottom(nadir_idx + 1) ~= 0
                    nadir_idx = nadir_idx + 1;
                else
                    nadir_idx = nadir_idx - 1;
                end
                
%                 subplot(211)
%                 pcolor(y,z,TS); shading flat; axis equal; set(gca,'ydir','reverse');
%                 set(gca,'fontname','Times'); caxis([-100 0]); colorbar
%                 title(['Ping ' num2str(pingidx) ': TS'])
%                 ylim([0 300])
%                 hold on
%                 plot(yBottom(pingidx,:),zBottom(pingidx,:),'r.')
%                 hold off 
%                 
%                 subplot(212)
%                 pcolor(y,z,Sv); shading flat; axis equal; set(gca,'ydir','reverse');
%                 set(gca,'fontname','Times'); caxis([-100 0]); colorbar
%                 title(['Ping ' num2str(pingidx) ' :Sv'])
%                 hold on
%                 plot(yBottom(pingidx,:),zBottom(pingidx,:),'r.')
%                 ylim([0 300])
%                 hold off
%                 drawnow
% %                 plot(y(bin_idx),z(bin_idx),'rs')
%                 pause(0.05)
                
                cut_Sv  = Sv(:,1:nadir_idx);
                cut_y   = y(:,1:nadir_idx);
                cut_z   = z(:,1:nadir_idx);
                
                bin_edges = -150:10:150;
                
%                 b_x = beam_x(pingidx,:);
%                 b_y = beam_y(pingidx,:);
%                 b_z = beam_z(pingidx,:);
                
                for i = 1:(length(bin_edges)-1)
                   bin_idx = cut_y < bin_edges(i+1) & cut_y > bin_edges(i);
          %          test_mean(i)    = mean(cut_Sv(bin_idx));
          %          test_max(i)     = max(cut_Sv(bin_idx));
                    test1 = movmean(cut_Sv(bin_idx),10);
                    test_mmax(i) = max(test1);
                    y_mid(i) =  mean(cut_y(bin_idx));
                end 
                
                y_mid_bin(pingidx,:) = y_mid;
                z_nadir(pingidx) = zBottom(nadir_idx);
                Sv_mmax(pingidx,:) = test_mmax;
                            
            end
        
        % all other datagrams, skip and move to next    
        else
            % char(73) = I, extra detects datagram
            % char(196) = A, attitude datagram
            
            fread(fid,sz-2,'uint8');
        end
    else
        % the start id was not equal to 2, so either this is not a file
        % containing simrad datagrams or we have made an error reading
        % it...
        disp('start id not equal to 2!')
        pause(.1)
        if feof(fid)
            stopp = 1;
        end
    end
end

fclose(fid)

%HeatMap_0030.nadir_z = z_nadir;
%HeatMap_0030.across_y = y_mid_bin;
%HeatMap_0030.Sv = Sv_mmax;

%save('HeatMap_0030.mat','HeatMap_0030')
   
