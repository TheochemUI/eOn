# Microsoft Developer Studio Generated NMAKE File, Format Version 4.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

!IF "$(CFG)" == ""
CFG=allib - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to allib - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "allib - Win32 Release" && "$(CFG)" != "allib - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "allib.mak" CFG="allib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "allib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "allib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 
################################################################################
# Begin Project
# PROP Target_Last_Scanned "allib - Win32 Debug"
F90=fl32.exe

!IF  "$(CFG)" == "allib - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
OUTDIR=.\Release
INTDIR=.\Release

ALL : "$(OUTDIR)\al.lib"

CLEAN : 
	-@erase "..\..\lib\al.lib"
	-@erase ".\Release\forces.obj"
	-@erase ".\Release\sumembforce.obj"
	-@erase ".\Release\embedenergy.obj"
	-@erase ".\Release\fofrhoDblexp.obj"
	-@erase ".\Release\potinit.obj"
	-@erase ".\Release\dfrhoDblexp.obj"
	-@erase ".\Release\gagafeDblexp.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE F90 /Ox /I "Release/" /c /nologo
# ADD F90 /Ox /4Yb /4L80 /I "Release/" /c /nologo /4Yaltparam /Ob2
F90_PROJ=/Ox /4Yb /4L80 /I "Release/" /c /nologo /4Yaltparam /Ob2 /Fo"Release/"\
 
F90_OBJS=.\Release/
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/allib.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"..\..\lib\al.lib"
LIB32_FLAGS=/nologo /out:"..\..\lib\al.lib" 
LIB32_OBJS= \
	"$(INTDIR)/forces.obj" \
	"$(INTDIR)/sumembforce.obj" \
	"$(INTDIR)/embedenergy.obj" \
	"$(INTDIR)/fofrhoDblexp.obj" \
	"$(INTDIR)/potinit.obj" \
	"$(INTDIR)/dfrhoDblexp.obj" \
	"$(INTDIR)/gagafeDblexp.obj"

"$(OUTDIR)\al.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "allib - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
OUTDIR=.\Debug
INTDIR=.\Debug

ALL : "$(OUTDIR)\allib.lib"

CLEAN : 
	-@erase ".\Debug\allib.lib"
	-@erase ".\Debug\gagafeDblexp.obj"
	-@erase ".\Debug\sumembforce.obj"
	-@erase ".\Debug\fofrhoDblexp.obj"
	-@erase ".\Debug\forces.obj"
	-@erase ".\Debug\embedenergy.obj"
	-@erase ".\Debug\potinit.obj"
	-@erase ".\Debug\dfrhoDblexp.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE F90 /Z7 /I "Debug/" /c /nologo
# ADD F90 /Z7 /I "Debug/" /c /nologo /Ob2
F90_PROJ=/Z7 /I "Debug/" /c /nologo /Ob2 /Fo"Debug/" 
F90_OBJS=.\Debug/
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/allib.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/allib.lib" 
LIB32_OBJS= \
	"$(INTDIR)/gagafeDblexp.obj" \
	"$(INTDIR)/sumembforce.obj" \
	"$(INTDIR)/fofrhoDblexp.obj" \
	"$(INTDIR)/forces.obj" \
	"$(INTDIR)/embedenergy.obj" \
	"$(INTDIR)/potinit.obj" \
	"$(INTDIR)/dfrhoDblexp.obj"

"$(OUTDIR)\allib.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

################################################################################
# Begin Target

# Name "allib - Win32 Release"
# Name "allib - Win32 Debug"

!IF  "$(CFG)" == "allib - Win32 Release"

!ELSEIF  "$(CFG)" == "allib - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=.\sumembforce.f
DEP_F90_SUMEM=\
	".\commonblks\parameters.cmn"\
	".\commonblks\comintlis.cmn"\
	

"$(INTDIR)\sumembforce.obj" : $(SOURCE) $(DEP_F90_SUMEM) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\forces.f
DEP_F90_FORCE=\
	".\commonblks\parameters.cmn"\
	".\commonblks\combaths.cmn"\
	".\commonblks\comgeom.cmn"\
	".\commonblks\comconf.cmn"\
	".\commonblks\comenergy.cmn"\
	".\commonblks\comenperat.cmn"\
	".\commonblks\compotent.cmn"\
	".\commonblks\comluns.cmn"\
	".\commonblks\comintlis.cmn"\
	

"$(INTDIR)\forces.obj" : $(SOURCE) $(DEP_F90_FORCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\fofrhoDblexp.f
DEP_F90_FOFRH=\
	".\commonblks\parameters.cmn"\
	".\commonblks\comconf.cmn"\
	".\commonblks\compotent.cmn"\
	

"$(INTDIR)\fofrhoDblexp.obj" : $(SOURCE) $(DEP_F90_FOFRH) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\embedenergy.f
DEP_F90_EMBED=\
	".\commonblks\parameters.cmn"\
	".\commonblks\comconf.cmn"\
	".\commonblks\comluns.cmn"\
	".\commonblks\comintlis.cmn"\
	

"$(INTDIR)\embedenergy.obj" : $(SOURCE) $(DEP_F90_EMBED) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\dfrhoDblexp.f
DEP_F90_DFRHO=\
	".\commonblks\parameters.cmn"\
	".\commonblks\comconf.cmn"\
	".\commonblks\compotent.cmn"\
	

"$(INTDIR)\dfrhoDblexp.obj" : $(SOURCE) $(DEP_F90_DFRHO) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\gagafeDblexp.f
DEP_F90_GAGAF=\
	".\commonblks\parameters.cmn"\
	".\commonblks\comgeom.cmn"\
	".\commonblks\combaths.cmn"\
	".\commonblks\compotent.cmn"\
	".\commonblks\comconf.cmn"\
	".\commonblks\comenperat.cmn"\
	".\commonblks\comluns.cmn"\
	".\commonblks\comintlis.cmn"\
	

"$(INTDIR)\gagafeDblexp.obj" : $(SOURCE) $(DEP_F90_GAGAF) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\potinit.f
DEP_F90_POTIN=\
	".\commonblks\parameters.cmn"\
	".\commonblks\combaths.cmn"\
	".\commonblks\comgeom.cmn"\
	".\commonblks\comconf.cmn"\
	".\commonblks\comenergy.cmn"\
	".\commonblks\comenperat.cmn"\
	".\commonblks\compotent.cmn"\
	".\commonblks\comluns.cmn"\
	".\commonblks\comintlis.cmn"\
	

"$(INTDIR)\potinit.obj" : $(SOURCE) $(DEP_F90_POTIN) "$(INTDIR)"


# End Source File
# End Target
# End Project
################################################################################
