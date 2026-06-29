% script to process series of files

% set manually the file directory
datadir='../../UNH_data/UNH_EM2040_40_Apr_13_2023/';

% pull list of kmwcd files
filelist=dir(fullfile(datadir,'*.kmwcd'));
Nfiles=length(filelist);

for i=1:Nfiles
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


