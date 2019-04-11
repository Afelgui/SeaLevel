# Microsoft Developer Studio Project File - Name="arc_ocn" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=arc_ocn - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "arc_ocn.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "arc_ocn.mak" CFG="arc_ocn - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "arc_ocn - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "arc_ocn - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "arc_ocn - Win32 Release"

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
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /real_size:64 /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x419 /d "NDEBUG"
# ADD RSC /l 0x419 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "arc_ocn - Win32 Debug"

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
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /debug:full /nologo /real_size:64 /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x419 /d "_DEBUG"
# ADD RSC /l 0x419 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /incremental:no /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "arc_ocn - Win32 Release"
# Name "arc_ocn - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\bkc.F
DEP_F90_BKC_F=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\density.F
DEP_F90_DENSI=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\difxz.F
DEP_F90_DIFXZ=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\difyz.F
DEP_F90_DIFYZ=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\fact1.F
DEP_F90_FACT1=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\factts.F
# End Source File
# Begin Source File

SOURCE=.\factuv.F
DEP_F90_FACTU=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\filter.for
DEP_F90_FILTE=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\filtr_lvl.for
DEP_F90_FILTR=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\filtr_ts.for
DEP_F90_FILTR_=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\integ.F
DEP_F90_INTEG=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\integ1.F
DEP_F90_INTEG1=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\n_ad.for
DEP_F90_N_AD_=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\n_ad1.F
DEP_F90_N_AD1=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\ocn.f
DEP_F90_OCN_F=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\org.F
DEP_F90_ORG_F=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\qstx.f
DEP_F90_QSTX_=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\qstx1.F
DEP_F90_QSTX1=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\quickest_uv_dina.F
DEP_F90_QUICK=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\quickest_uv_dina1.F
DEP_F90_QUICKE=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\rddisk.F
DEP_F90_RDDIS=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\right1.for
DEP_F90_RIGHT=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\rivers.F
DEP_F90_RIVER=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\stbet0.F
DEP_F90_STBET=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\stbet0_backup.F
DEP_F90_STBET0=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\tras1.for
DEP_F90_TRAS1=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\tras23.for
DEP_F90_TRAS2=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\tras_difxz.F
DEP_F90_TRAS_=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\tras_difyz.F
DEP_F90_TRAS_D=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\vyaz_uv_dina.F
DEP_F90_VYAZ_=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\wind_bg.F
DEP_F90_WIND_=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\wtdisk.F
DEP_F90_WTDIS=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\yflux.f
DEP_F90_YFLUX=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\yflux1.F
DEP_F90_YFLUX1=\
	".\model.par"\
	
# End Source File
# Begin Source File

SOURCE=.\z_adv.F
DEP_F90_Z_ADV=\
	".\model.par"\
	
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# Begin Source File

SOURCE=.\model.par
# End Source File
# End Target
# End Project
