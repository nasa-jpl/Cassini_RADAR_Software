#======================================================================
# @(#) $Id: Makefile,v 11.7 2017/04/10 16:49:00 cveerama Exp $
# Makefile for RAS programs
#======================================================================

SHELL = /bin/csh
INCDIR = ../../include
SPICEDIR = ../../cspice/src/cspice
INCFLAGS = -I$(INCDIR) -I$(SPICEDIR) -I.
DEBUGFLAGS = -ggdb2

LIBS = ../../lib/libsupport.a ../../cspice/lib/cspice.a
LDLIBS = $(LIBS) -lm -ljpeg -ltiff
INSTALLDIR = ../../bin

CPPFLAGS = $(SPECIAL_USER_FLAGS) $(DEBUGFLAGS) $(INCFLAGS) -Wall
FFLAGS = -ffixed-line-length-none
SYSLIBS = -lm
CCC = g++
CC = gcc
FF = g77


# Main downlink software
MAIN_DOWNLINK = \
  decode_ccsd \
  preproc \
  sar_proc \
  subsample_bidr \
  remove_nan \
  flip_bidr \
  estrd_to_cfg \
  average_down_bidr \
  sar_topo \
  sar_stereo 


# Downlink support software
DOWNLINK_SUPPORT = \
  generate_geometry_file \
  clip_header \
  compute_kp_bidr \
  compute_snr_bidr \
  cut_binary_record \
  sar_index_to_dopran \
  find_pixel \
  getint_from_bidr \
  getfloat_from_bidr \
  rowcol_from_backplane \
  mean_std_from_bidr \
  backplane_to_plot \
  backplane_to_plot_vs_distance \
  backplane_to_plot_vs_time \
  compute_saralt_bidr \
  estrd_anal \
  topomap2jpeg \
  compute_bidr_bound_box \
  make_bounds_str \
  topomap \
  get_tiepoint_info \
  get_tiepoint_info_from_lbdr \
  make_extra_backplane 


# Downlink testing software
DOWNLINK_TEST = \
  burst_data_extract \
  monitor \
  read_lbdr \
  read_sbdr \
  sar_stereo




# Miscellaneous
MISC = \
  getClosestApproach \
  bidr_latlon \
  bidr_rowcol \
  bidr2jpeg \
  bidr2tiff \
  subsample_beammask \
  update_config \
  zero_run_float_zip


ALL = \
         $(MAIN_DOWNLINK) $(DOWNLINK_SUPPORT) $(DOWNLINK_TEST) \
         $(MISC) 

.PRECIOUS: %.cpp
.PRECIOUS: %.F
.PRECIOUS: %.pl
.SUFFIXES: .cpp
.SUFFIXES: .F

# Implicit rules

.cpp.o:
	$(CC) $(CPPFLAGS) $<

.F.o:
	$(FF) $(FFLAGS) $<

# High level make targets

default: $(ALL)
down: $(MAIN_DOWNLINK) $(DOWNLINK_SUPPORT) $(DOWNLINK_TEST)
other: $(MISC)


# Explicit rules


%: %.cpp $(LIBS) 
	$(CCC) $(CPPFLAGS) -o $@ $@.cpp $(LDLIBS)

%: %.F  
	$(FF) $(FFLAGS) -o $@ $@.F 

install: 
	mv -f $(ALL) $(INSTALLDIR)

clean:
	- /bin/rm -f $(ALL)
	- /bin/rm -f core


# handle executables that differ only by command name
bidr_rowcol: bidr_latlon
	\cp bidr_latlon bidr_rowcol

get_tiepoint_info_from_lbdr: get_tiepoint_info
	\cp get_tiepoint_info get_tiepoint_info_from_lbdr
