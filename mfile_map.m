% m-file map!

%% ORIGINAL setup
%   - one script does almost everything including basic visualization
%   - changing filename or data directory hardcoded and updates github
%   - files included:
%       readEM_kmall.m, sum_viz.m, read_bin_kmall.m, parse_metadata.m

% readEM_kmall.m - original m-file for processing data
%   calls read_bin_kmall.m
%   saves outmatfile - which has the water column data structure
%         outvizfile - which has grids for visualization
%         outallstruct - which saves the basic input metadata

% read_bin_kmall.m - reads the binary data in the .kmwcd file, returns
% binary data
%   calls n/a
%   saves n/a

% sum_viz.m
% parse_metadata.m 

%% NEW setup 
%   - changed to be more function based and flexible
%   - filenames and directories specified in function call not code
%   - files included:
%       read_kmall.m, check_sampfinfo.m, check_sonar_param.m,
%       parse_install_txt.m
%
%   - COFFEE library items needed:
%       CFF_read_kmall,

% read_kmall.m - new m-file to process data - reads meta data
%   calls parse_install_txt.m
%   saves 

% parse_install_txt.m              

% check_sampinfo.m       

% check_sonar_param.m

%% undetermined files:                    
detect_bottom.m        
process_file_series.m  
fill_grid.m                    
get_kmwcd_freq.m       