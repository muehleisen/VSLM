# VSLM
VSLM - The Virtual Sound Level Meter

This project is the MATLAB development of a virtual sound level meter. The program will read in a calibrated .wav file and allow the user to analyze it as one would analyze a sound field with a sound level meter. The software implements Fast, Slow, Impulse and LEQ; A, C and Flat Weighting; Ln and Noise Dose analysis; Octave and 1/3 Octave band analysis; high resolution FFT analysis and spectrograms. Band analysis can be made using fast FFT methods or slower, but ANSI standard methods.
Source code requires MATLAB + Signal Processing Toolbox.
Developed under MATLAB R2010a (but might work with some earlier versions)
Executable version does not require MATLAB and runs under Windows XP, Vista, and Windows 7.

This software does not do real-time sound level measurement. It does not read from a sound card.
It has been developed for the analysis of calibrated digital recordings.

## Latest Developments
The software archive has been moved from sourceforge to github in an attempt to increase collaborative additions to the software.
If you want to contribute to development of VSLM, please learn about forking and cloning of a git repository, fork or clone the VSLM repository,
and become part of the community.

Fork:
https://github.com/muehleisen/VSLM#fork-destination-box
Clone
https://github.com/muehleisen/VSLM.git

If you just want the most up-to-date version of the software to download as a zip file you can go to
https://github.com/muehleisen/VSLM/archive/master.zip

## Version

VSLM is currently version 0.4.2 but the documentation has not been updated yet.



## WINDOWS EXECUTABLE USERS:
Windows users without MATLAB and the signal processing toolbox installed should download the vslm_0_4_1_pkg.exe.  
This file includes the windows executables, help files, and MATLAB runtime libraries from the releases part of the github site.
https://github.com/muehleisen/VSLM/releases
This is not the most up-to-date version of the software and it is not known if this will install and execute properly on Windows 10.

Move the vslm_0_4_1_pkg.exe file to the directory where you want the vslm executable to be installed (create the directory if necessary) and run the file.  This will extract the contents into the directory and install the MATLAB runtime library if necessary

To start vslm run the vslm.exe executable.  This is actually a wrapper program that runs the true executable, matgui, and shows asplash screen while the MATLAB runtime library loads (this can take a while).

The documentation includes vslm.chm, a windows helpfile format help file
vslm.htm and the /files subdirectory which contains an html format help file
vslm.pdf which is a PDF version of the same documentation
vslm.hnd is the HelpNDoc file which is the source file for all the documentation formats

The archive also includes vslm_0_4_1.fig and vslm_0_4_1.m which are the MATLAB source files

## MATLAB USERS

Matlab users should download the vslm_0_4_1_matlab_source.zip and 
vslm_0_4_1_docs.zip and unzip both to the same directory.  
Users should then run vslm_0_4_1.m to run vslm.

vslm_0_4_1.m requires the signal processing toolbox.

vslm_0_4_1.m was developed using MATLAB R2010.  

I do not know what earlier versions with which it can run.
