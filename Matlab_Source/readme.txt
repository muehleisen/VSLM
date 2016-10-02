Download and Install Information for vslm the Virtual Sound Level Meter.

WINDOWS EXECUTABLE USERS:
Windows users without MATLAB and the signal processing toolbox installed should download the vslm_0_4_1_pkg.exe.  
This file includes the windows executables, help files, and MATLAB runtime libraries.

Move the vslm_0_4_1_pkg.exe file to the directory where you want the vslm executable to be installed (create the directory if necessary) and run the file.  This will extract the contents into the directory and install the MATLAB runtime library if necessary

To start vslm run the vslm.exe executable.  This is actually a wrapper program that runs the true executable, matgui, and shows asplash screen while the MATLAB runtime library loads (this can take a while).

The documentation includes vslm.chm, a windows helpfile format help file
vslm.htm and the /files subdirectory which contains an html format help file
vslm.pdf which is a PDF version of the same documentation
vslm.hnd is the HelpNDoc file which is the source file for all the documentation formats

The archive also includes vslm_0_4_1.fig and vslm_0_4_1.m which are the MATLAB source files

MATLAB USERS

Matlab users should download the vslm_0_4_1_matlab_source.zip and 
vslm_0_4_1_docs.zip and unzip both to the same directory.  
Users should then run vslm_0_4_1.m to run vslm.

vslm_0_4_1.m requires the signal processing toolbox.

vslm_0_4_1.m was developed using MATLAB R2010.  

I do not know what earlier versions with which it can run.