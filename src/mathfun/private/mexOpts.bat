@echo off
rem MSVC50OPTS.BAT
rem
rem    Compile and link options used for building MEX-files
rem    using the Microsoft Visual C++ compiler version 5.0 
rem
rem    $Revision: 1.5 $  $Date: 1997/12/02 16:00:40 $
rem
rem ********************************************************************
rem General parameters
rem ********************************************************************

set MATLAB=%MATLAB%
set MSVC_ROOT=c:\DevStudio
set PATH=%MSVC_ROOT%\VC\BIN;%MSVC_ROOT%\sharedIDE\bin;%PATH%
set INCLUDE=%MSVC_ROOT%\VC\INCLUDE;%MSVC_ROOT%\VC\MFC\INCLUDE;%MSVC_ROOT%\VC\ATL\INCLUDE;%INCLUDE%
set LIB=%MSVC_ROOT%\VC\LIB;%MSVC_ROOT%\VC\MFC\LIB;%LIB%

rem *******************************************************************
rem Add include of the libraries used
rem *******************************************************************
set INCLUDE=c:\usr\danuser\soft\include;c:\usr\danuser\matlab\mexCalls;%INCLUDE%

rem ********************************************************************
rem Compiler parameters
rem ********************************************************************
set COMPILER=cl
set COMPFLAGS=-c -Zp8 -G5 -W3 -WX -DMATLAB_MEX_FILE 
set OPTIMFLAGS=-O2
set DEBUGFLAGS=-Zi

rem ********************************************************************
rem Library creation command
rem ********************************************************************
set PRELINK_CMDS=lib /def:"%MATLAB%\extern\include\matlab.def" /machine:ix86 /OUT:%LIB_NAME%1.lib 
set PRELINK_DLLS=lib /def:"%MATLAB%\extern\include\%DLL_NAME%.def" /machine:ix86 /OUT:%DLL_NAME%.lib	

rem ********************************************************************
rem Linker parameters
rem ********************************************************************
set LINKER=link
set LINKFLAGS=/dll /export:mexFunction %LIB_NAME%1.lib /implib:%LIB_NAME%l.lib /libpath:"c:\usr\danuser\soft\lib"
set LINKOPTIMFLAGS= 
set LINKDEBUGFLAGS=/debug
set LINK_FILE=
set LINK_LIB=
set NAME_OUTPUT=/out:%MEX_NAME%.dll

rem ********************************************************************
rem Resource compiler parameters
rem ********************************************************************
set RC_COMPILER=rc /fo mexversion.res
set RC_LINKER= 



