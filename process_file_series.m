% script to process series of files

% set manually the file directory
datadir='../../UNH_data/UNH_EM2040_40_Apr_13_2023/';

% pull list of kmwcd files
filelist=dir(fullfile(datadir,'*.kmwcd'));
Nfiles=length(filelist);

% set file number to start or stop at
istart=50;
iend=Nfiles;

% this section pulls central frequency information
for i=istart:iend
    % get ith kmwcd filename
    kmwcdfile=filelist(i).name;
    [~,filecode,~]=fileparts(kmwcdfile);

    % check if matching kmall file and set accordingly
    if isfile(fullfile(datadir,filecode,'.kmall'))
        kmallfile=fullfile(datadir,filecode,'.kmall');
    else
        kmallfile=[];
    end

    % pull metadata and print to screen
    get_kmwcd_freq(datadir,kmwcdfile,kmallfile)

    % pause before next one so can enter data in spreadsheet
    % as of 6/29/26, get_kmwcd_freq only prints to screen
    pause
end

% need to add sections to extract data and pull characteristics from
% extracted data
% may want a switch structure to choose between


