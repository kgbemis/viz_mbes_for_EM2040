The code in this package is intended to help read and visualize EM2040 water column data.  

It calls several functions from the package CoFFee (https://github.com/alexschimel/CoFFee) 
to read and parse the .kmall files.  It assumes that CoFFee has been downloaded and that the 
path to the CoFFee folder and its subfolders has been added to the Matlab path. I usually 
put it in the Matlab toolbox folder and add the path using the Matlab interface such that
it is always available.  But I think it will work if CoFFee is somewhere else as long as 
the path is properly added to Matlab.

Things you will need to change for it too work:
1 - update the path to your data files: filelocation='your/path/dir/';  
    Note that this can be either an absolute or a relative path.
2 - 

Implemented functionality
- gathers info from kmall file
- displays number and type of datagrams
- can extract start and end times
- can verify if any split datagrams
- 

----
latest readme file update - kgb 4/4/23
