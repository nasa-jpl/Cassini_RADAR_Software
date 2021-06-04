#======================================================================
# @(#) $Id: Makefile,v 11.6 2017/04/07 17:03:08 richw Exp $
# Makefile for RAS support library
#======================================================================

SHELL = /bin/csh
INCDIR = ../../include
TARGET_LIB = ../../lib/libsupport.a
SPICEDIR = ../../cspice/src/cspice
INCFLAGS = -I$(SPICEDIR)
DEBUGFLAGS = -ggdb2 

CPPFLAGS = -c $(SPECIAL_USER_FLAGS) $(DEBUGFLAGS) $(INCFLAGS) -Wall
SYSLIBS = -lm
CCC = gcc
CC = g++

# Object file list
OBJS = Io.o Units.o Sab.o Config.o Plot.o RadiometerData.o TextMgr.o \
  Frame.o Time.o TargetGeom.o Utils.o Beam.o Ivd.o Ckernel.o  \
  CassiniSim.o GPHData.o AmbiguityData.o Ieb.o BurstData.o \
  L1B.o Flyby.o IebList.o Ambiguity.o AmbiguityFile.o \
  PointTargetSim.o PerfFile.o SARFunctions.o Profiling.o Baq.o \
  Distributions.o L1A.o MFile.o rmss_functions.o IOinput.o InstructDensity.o \
  Interpolate.o	IebProfile.o DebugInfo.o BIDR.o Projections.o LatLonGrid.o \
  PDSLabel.o SARProcParams.o SincTable.o L1I.o RasList.o RasProc.o SimpleArray.o \
  DopplerCentroid.o CAGTable.o RADFunctions.o RangeDopplerFile.o BIDRFile.o \
  PDSVolume.o Spline.o TopoMap.o

# Include file list
INCS = Io.h Units.h Sab.h Config.h Plot.h RadiometerData.h TextMgr.h \
  Frame.h Time.h TargetGeom.h Utils.h Beam.h Array.h Array_impl.h \
  Constants.h Error.h Units_impl1.h Units_impl2.h Ivd.h Ckernel.h   \
  CassiniSim.h GPHData.h AmbiguityData.h TemplateUtils.h \
  Ieb.h BurstData.h L1B.h config_keywords.h Flyby.h IebList.h \
  Ambiguity.h AmbiguityFile.h PointTargetSim.h PerfFile.h \
  SARFunctions.h Baq.h Profiling.h Distributions.h L1A.h MFile.h \
  rmss_functions.h IOinput.h EGSEdefinitions.h InstructDensity.h \
  Interpolate.h IebProfile.h SARProcParams.h DebugInfo.h BIDR.h LatLonGrid.h \
  Projections.h PDSLabel.h  SincTable.h L1I.h RasList.h RasProc.h SimpleArray.h \
  DopplerCentroid.h CAGTable.h RADFunctions.h Radiometer.h RangeDopplerFile.h \
  BIDRFile.h PDSVolume.h Spline.h TopoMap.h

.PRECIOUS: %.cpp
.SUFFIXES: .cpp

# Implicit rules

.cpp.o:
	$(CC) $(CPPFLAGS) $<

# Explicit rules

default: $(TARGET_LIB)

$(TARGET_LIB): $(OBJS) $(INCS)
	- /bin/rm -f $(TARGET_LIB)
	ar -r $(TARGET_LIB) $(OBJS)
	ranlib $(TARGET_LIB)
	- (set incfiles=($(INCS)); cd $(INCDIR); /bin/rm -f $$incfiles)
	cp $(INCS) $(INCDIR)

#  Obsolete rules for extracting from RCS revision control
#%.cpp: $(RAS_CENTRAL_TREE)/src/objs/RCS/%.cpp,v
#	co $(RAS_CENTRAL_TREE)/src/objs/RCS/$@,v
#
#%.h: $(RAS_CENTRAL_TREE)/src/objs/RCS/%.h,v
#	co $(RAS_CENTRAL_TREE)/src/objs/RCS/$@,v

clean:
	- /bin/rm -f $(TARGET_LIB)
	- /bin/rm -f $(OBJS)
	- /bin/rm -f core
	- (set incfiles=($(INCS)); cd $(INCDIR); /bin/rm -f $$incfiles)

# List out dependencies
depend:
	@@ inclist `echo $(OBJS) | sed -e 's/[.]o/.cpp/g'`

# Dependencies from make depend
Io.o: Io.cpp Io.h Units.h Error.h Units_impl1.h
Units.o: Units.cpp Units.h Io.h Constants.h Utils.h Error.h Units_impl1.h Frame.h Array.h Time.h Config.h Ieb.h Array_impl.h CAGTable.h
Sab.o: Sab.cpp Sab.h Array.h Error.h Io.h Units.h Plot.h Constants.h Ieb.h CAGTable.h Utils.h Array_impl.h Units_impl1.h Frame.h Time.h Config.h
Config.o: Config.cpp Units.h Error.h Config.h Constants.h Utils.h DebugInfo.h Io.h Units_impl1.h Frame.h Array.h Time.h Ieb.h Array_impl.h CAGTable.h
Plot.o: Plot.cpp Plot.h Io.h Array.h Error.h Units.h Utils.h Array_impl.h Units_impl1.h Frame.h Time.h Config.h Ieb.h CAGTable.h
RadiometerData.o: RadiometerData.cpp RadiometerData.h L1B.h Array.h Units.h Error.h Sab.h BurstData.h GPHData.h Frame.h Io.h Utils.h Array_impl.h Units_impl1.h Ieb.h CAGTable.h Time.h TopoMap.h Config.h
TextMgr.o: TextMgr.cpp TextMgr.h Error.h
Frame.o: Frame.cpp Frame.h Time.h Constants.h Profiling.h Interpolate.h TargetGeom.h Flyby.h SARProcParams.h DebugInfo.h Config.h Error.h Units.h Array.h BurstData.h Projections.h Ambiguity.h RasList.h Io.h Units_impl1.h Utils.h Array_impl.h GPHData.h CAGTable.h TopoMap.h Beam.h SimpleArray.h Ieb.h Sab.h
Time.o: Time.cpp Error.h Time.h Units.h TemplateUtils.h Io.h Units_impl1.h Frame.h Array.h Config.h Utils.h Array_impl.h Ieb.h CAGTable.h
TargetGeom.o: TargetGeom.cpp TargetGeom.h Frame.h Constants.h TemplateUtils.h Error.h Units.h Time.h Config.h Array.h Io.h Units_impl1.h Utils.h Array_impl.h Ieb.h CAGTable.h
Utils.o: Utils.cpp Frame.h Units.h Time.h Constants.h TemplateUtils.h SARFunctions.h TargetGeom.h Plot.h SimpleArray.h Utils.h Config.h Error.h Array.h Io.h Units_impl1.h Ieb.h Array_impl.h CAGTable.h
Beam.o: Beam.cpp Error.h Units.h Beam.h Array.h Constants.h DebugInfo.h Utils.h Io.h Units_impl1.h Config.h Array_impl.h Frame.h Time.h Ieb.h CAGTable.h
Ivd.o: Ivd.cpp Ivd.h Array.h Config.h Io.h Units.h Error.h Frame.h Time.h Constants.h TemplateUtils.h Utils.h Array_impl.h Units_impl1.h Ieb.h CAGTable.h
Ckernel.o: Ckernel.cpp Ivd.h Array.h Config.h Io.h Units.h Error.h Frame.h Time.h Constants.h Ckernel.h Utils.h Array_impl.h Units_impl1.h Ieb.h CAGTable.h
CassiniSim.o: CassiniSim.cpp CassiniSim.h Constants.h Array.h DebugInfo.h Frame.h Beam.h Units.h Config.h Time.h TargetGeom.h SARFunctions.h Error.h Io.h Utils.h Array_impl.h Units_impl1.h Ieb.h CAGTable.h
GPHData.o: GPHData.cpp GPHData.h Units.h Io.h Error.h Array.h Time.h Units_impl1.h Utils.h Array_impl.h Frame.h Config.h Ieb.h CAGTable.h
AmbiguityData.o: AmbiguityData.cpp AmbiguityData.h TemplateUtils.h Array.h Units.h Error.h Frame.h Time.h Config.h Io.h Utils.h Array_impl.h Units_impl1.h Ieb.h CAGTable.h
Ieb.o: Ieb.cpp Error.h Units.h Ieb.h Time.h Array.h Constants.h Utils.h Io.h Units_impl1.h Config.h CAGTable.h Array_impl.h Frame.h
BurstData.o: BurstData.cpp BurstData.h Error.h Constants.h config_keywords.h TargetGeom.h Frame.h Flyby.h Units.h Array.h GPHData.h Time.h CAGTable.h TopoMap.h Config.h Io.h Units_impl1.h Utils.h Array_impl.h Ieb.h
L1B.o: L1B.cpp L1B.h Error.h config_keywords.h L1I.h PDSLabel.h Constants.h IebProfile.h Units.h Sab.h Array.h BurstData.h GPHData.h Frame.h SimpleArray.h CassiniSim.h MFile.h SARProcParams.h SincTable.h DebugInfo.h RangeDopplerFile.h Baq.h RasList.h TopoMap.h Config.h Io.h Units_impl1.h Ieb.h CAGTable.h Utils.h Array_impl.h Time.h Beam.h TargetGeom.h SARFunctions.h Projections.h Flyby.h Ambiguity.h
Flyby.o: Flyby.cpp Units.h Error.h Config.h Frame.h Time.h TargetGeom.h TemplateUtils.h Ckernel.h Flyby.h Plot.h Utils.h Io.h Units_impl1.h Array.h Ivd.h Constants.h Ieb.h Array_impl.h CAGTable.h
IebList.o: IebList.cpp Error.h Utils.h Units.h Ieb.h Time.h Array.h Constants.h Config.h IebList.h Ambiguity.h Plot.h IebProfile.h Flyby.h Frame.h Io.h Units_impl1.h CAGTable.h Array_impl.h TargetGeom.h Beam.h SimpleArray.h
Ambiguity.o: Ambiguity.cpp Frame.h TargetGeom.h Units.h Time.h Constants.h TemplateUtils.h Ambiguity.h Config.h Array.h Io.h Utils.h Error.h Units_impl1.h Beam.h SimpleArray.h Ieb.h Array_impl.h CAGTable.h
AmbiguityFile.o: AmbiguityFile.cpp Error.h Units.h Time.h Array.h Constants.h Config.h AmbiguityFile.h Io.h Utils.h Units_impl1.h Array_impl.h Frame.h TargetGeom.h Beam.h Ambiguity.h Ieb.h SimpleArray.h CAGTable.h
PointTargetSim.o: PointTargetSim.cpp PointTargetSim.h Constants.h config_keywords.h Units.h Array.h SARFunctions.h Distributions.h DebugInfo.h CassiniSim.h L1B.h Profiling.h Baq.h Ieb.h IebList.h Flyby.h MFile.h CAGTable.h Frame.h Error.h Io.h Units_impl1.h Utils.h Array_impl.h Config.h Beam.h Time.h TargetGeom.h Sab.h BurstData.h GPHData.h Ambiguity.h IebProfile.h TopoMap.h SimpleArray.h
PerfFile.o: PerfFile.cpp Error.h Units.h Ambiguity.h Time.h Array.h Constants.h Config.h PerfFile.h Utils.h Io.h Units_impl1.h Beam.h Frame.h SimpleArray.h Ieb.h Array_impl.h CAGTable.h
SARFunctions.o: SARFunctions.cpp SARFunctions.h Constants.h TargetGeom.h Plot.h Baq.h Ieb.h SimpleArray.h Units.h Error.h Array.h Frame.h Time.h Config.h Beam.h Sab.h CAGTable.h Io.h Units_impl1.h Utils.h Array_impl.h
Profiling.o: Profiling.cpp 
Baq.o: Baq.cpp Array.h Error.h Units.h Time.h Config.h Frame.h TargetGeom.h Beam.h Baq.h Io.h Ieb.h Sab.h Utils.h Array_impl.h Units_impl1.h CAGTable.h
Distributions.o: Distributions.cpp Distributions.h Constants.h Units.h Frame.h Error.h Io.h Units_impl1.h Config.h Time.h Array.h Utils.h Array_impl.h Ieb.h CAGTable.h
L1A.o: L1A.cpp L1A.h Error.h Units.h Sab.h Array.h BurstData.h Io.h Units_impl1.h Ieb.h CAGTable.h Utils.h Array_impl.h GPHData.h Frame.h Time.h TopoMap.h Config.h
MFile.o: MFile.cpp MFile.h Error.h Frame.h Utils.h Config.h Units.h Time.h Array.h Ieb.h Io.h Units_impl1.h Array_impl.h CAGTable.h
rmss_functions.o: rmss_functions.cpp rmss_functions.h
IOinput.o: IOinput.cpp IOinput.h
InstructDensity.o: InstructDensity.cpp Array.h Error.h Units.h Time.h InstructDensity.h TargetGeom.h Frame.h Io.h Utils.h Array_impl.h Units_impl1.h Config.h Ieb.h CAGTable.h
Interpolate.o: Interpolate.cpp Interpolate.h
IebProfile.o: IebProfile.cpp Units.h Error.h Config.h IebProfile.h Plot.h Utils.h TargetGeom.h Frame.h Io.h Units_impl1.h Array.h Constants.h Time.h Ieb.h Array_impl.h CAGTable.h
DebugInfo.o: DebugInfo.cpp DebugInfo.h Error.h Frame.h Utils.h Config.h Units.h Time.h Array.h Ieb.h Io.h Units_impl1.h Array_impl.h CAGTable.h
BIDR.o: BIDR.cpp BIDR.h config_keywords.h Constants.h SimpleArray.h Utils.h PDSLabel.h L1I.h SARProcParams.h Time.h Io.h LatLonGrid.h Projections.h Units.h Frame.h Error.h Array.h Config.h Ieb.h L1B.h BurstData.h GPHData.h CassiniSim.h MFile.h SincTable.h DebugInfo.h RangeDopplerFile.h Baq.h RasList.h TopoMap.h Flyby.h Ambiguity.h Units_impl1.h Array_impl.h CAGTable.h Sab.h Beam.h TargetGeom.h SARFunctions.h
Projections.o: Projections.cpp SARFunctions.h Projections.h Constants.h DebugInfo.h Units.h Error.h Array.h Ieb.h Frame.h Config.h Io.h Units_impl1.h Utils.h Array_impl.h Time.h CAGTable.h
LatLonGrid.o: LatLonGrid.cpp LatLonGrid.h Error.h Constants.h DebugInfo.h Units.h Frame.h Config.h Io.h Units_impl1.h Time.h Array.h Utils.h Array_impl.h Ieb.h CAGTable.h
PDSLabel.o: PDSLabel.cpp DebugInfo.h Error.h PDSLabel.h Units.h Frame.h Utils.h Config.h Io.h Units_impl1.h Time.h Array.h Ieb.h Array_impl.h CAGTable.h
SARProcParams.o: SARProcParams.cpp SARProcParams.h DebugInfo.h Utils.h Constants.h BurstData.h L1I.h Distributions.h Config.h Units.h Projections.h Flyby.h Ambiguity.h RasList.h Error.h Frame.h Array.h Time.h Ieb.h GPHData.h CAGTable.h TopoMap.h SimpleArray.h L1B.h CassiniSim.h MFile.h SincTable.h RangeDopplerFile.h Baq.h Io.h Units_impl1.h TargetGeom.h Beam.h Sab.h Array_impl.h SARFunctions.h
SincTable.o: SincTable.cpp Config.h Constants.h DebugInfo.h Error.h SimpleArray.h SincTable.h Utils.h Units.h Frame.h Array.h Time.h Ieb.h Io.h Units_impl1.h Array_impl.h CAGTable.h
L1I.o: L1I.cpp L1I.h Error.h SARFunctions.h RangeDopplerFile.h config_keywords.h Constants.h DebugInfo.h BIDR.h BIDRFile.h Units.h SimpleArray.h L1B.h BurstData.h GPHData.h Frame.h CassiniSim.h MFile.h SARProcParams.h SincTable.h Baq.h RasList.h TopoMap.h Array.h Ieb.h Time.h Io.h Config.h PDSLabel.h LatLonGrid.h Projections.h Plot.h Units_impl1.h Sab.h CAGTable.h Beam.h TargetGeom.h Flyby.h Ambiguity.h Utils.h Array_impl.h
RasList.o: RasList.cpp Error.h Units.h Ieb.h Time.h Array.h Constants.h Config.h RasList.h Plot.h SARFunctions.h Baq.h Frame.h Utils.h Io.h Units_impl1.h CAGTable.h Array_impl.h Sab.h TargetGeom.h Beam.h
RasProc.o: RasProc.cpp Array.h Error.h Units.h Time.h RasProc.h SARFunctions.h Plot.h Io.h Utils.h Array_impl.h Units_impl1.h RasList.h Ieb.h Frame.h Config.h Constants.h Sab.h CAGTable.h
SimpleArray.o: SimpleArray.cpp SimpleArray.h
DopplerCentroid.o: DopplerCentroid.cpp Array.h Error.h Units.h DopplerCentroid.h Plot.h Constants.h TargetGeom.h Utils.h Ivd.h Ckernel.h Io.h Array_impl.h Units_impl1.h Frame.h Time.h Config.h Ieb.h CAGTable.h
CAGTable.o: CAGTable.cpp CAGTable.h Error.h
RADFunctions.o: RADFunctions.cpp RADFunctions.h Units.h Radiometer.h TemplateUtils.h Error.h Array.h Io.h Units_impl1.h Frame.h Time.h Config.h Utils.h Array_impl.h Ieb.h CAGTable.h
RangeDopplerFile.o: RangeDopplerFile.cpp Error.h Units.h Time.h Array.h Constants.h Config.h Plot.h Utils.h RangeDopplerFile.h Io.h Units_impl1.h Array_impl.h Frame.h Ieb.h CAGTable.h
BIDRFile.o: BIDRFile.cpp BIDRFile.h DebugInfo.h Error.h PDSLabel.h SimpleArray.h BIDR.h Constants.h Io.h Plot.h SARFunctions.h Config.h L1I.h SARProcParams.h Time.h LatLonGrid.h Projections.h Units.h Frame.h Array.h Ieb.h L1B.h BurstData.h GPHData.h CassiniSim.h MFile.h SincTable.h RangeDopplerFile.h Baq.h RasList.h TopoMap.h Flyby.h Ambiguity.h Units_impl1.h Utils.h Array_impl.h CAGTable.h Sab.h Beam.h TargetGeom.h
PDSVolume.o: PDSVolume.cpp DebugInfo.h Error.h PDSLabel.h PDSVolume.h Units.h Frame.h Utils.h Config.h Io.h Units_impl1.h Time.h Array.h Ieb.h Array_impl.h CAGTable.h
Spline.o: Spline.cpp Spline.h Units.h Error.h Array.h Io.h Units_impl1.h Utils.h Array_impl.h Frame.h Time.h Config.h Ieb.h CAGTable.h
TopoMap.o: TopoMap.cpp TopoMap.h Utils.h SimpleArray.h Constants.h BIDRFile.h Config.h Error.h Units.h Frame.h Array.h Time.h Ieb.h BIDR.h Io.h PDSLabel.h Plot.h SARFunctions.h Units_impl1.h Array_impl.h CAGTable.h L1I.h SARProcParams.h LatLonGrid.h Projections.h L1B.h BurstData.h GPHData.h CassiniSim.h MFile.h SincTable.h DebugInfo.h RangeDopplerFile.h Baq.h RasList.h Flyby.h Ambiguity.h Sab.h Beam.h TargetGeom.h
