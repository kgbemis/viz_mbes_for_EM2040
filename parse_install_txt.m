function y=parse_install_txt(iipdat)
% parses the install text in a IIP datagram

install_txt=iipdat.install_txt;

newtxt=splitlines(install_txt);

% display(newtxt)

% find the Tx serial number listing and pull the beamwidth from here
%   will use this to verify sector data value
indSerial=find(contains(newtxt,'SERIALno:'));
if ~isempty(indSerial)
    txstr=newtxt{indSerial+1};
    if strcmp(txstr(1:2),'TX')
        A=sscanf(txstr,'TX:%d;%fdeg,');
    end
    tx_beam_width=A(2);
    fprintf('tx beam width %f\n',tx_beam_width)
else
    fprintf('error: no serial no listing\n')
end

% find the system name and run type in the system entry 
indSerial=find(contains(newtxt,'SYSTEM:'));
if ~isempty(indSerial)
    sysstr=newtxt{indSerial};
    systxt=split(sysstr,',');
    indTRAI_RX1=find(contains(systxt,'TRAI_RX1:'));
    rxstr=systxt{indTRAI_RX1};
    lasttxt=split(rxstr,';');
    indSysName=find(contains(lasttxt,'W='));
    namestr=lasttxt{indSysName};
    if strcmp(namestr(1:2),'W=')
        sys_name=namestr(3:end);
    else
        printf('system name field has unexpected beginning')
    end
    fprintf('system name: %s\n',sys_name)
    runstr=systxt{1};
    if strcmp(runstr(1:7),'SYSTEM:')
        run_type=runstr(8:end);
    else
        printf('run type field has unexpected beginning')
    end
    fprintf('run type: %s\n',run_type)    
else
    fprintf('error: no ssystem listing\n')
end
y.tx_check=tx_beam_width;
y.sys_name=sys_name;
y.run_type=run_type;

