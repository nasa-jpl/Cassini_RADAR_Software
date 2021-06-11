#
# $Id: Makefile,v 11.1 2011/09/16 00:03:36 richw Exp $
#
# Makefile for CALIB directory
#

CHMOD	= /bin/chmod
MKDIR	= /bin/mkdir
RM	= /bin/rm -f
RMDIR	= /bin/rmdir

ANTENNA_PATTERN_VERSION	= 01
ANTENNA_PATTERN_STATIC_FILES = \
	BEAM1_V$(ANTENNA_PATTERN_VERSION).LBL \
	BEAM2_V$(ANTENNA_PATTERN_VERSION).LBL \
	BEAM3_V$(ANTENNA_PATTERN_VERSION).LBL \
	BEAM4_V$(ANTENNA_PATTERN_VERSION).LBL \
	BEAM5_V$(ANTENNA_PATTERN_VERSION).LBL \
	PATINFO.TXT
ANTENNA_PATTERN_DIR	= BEAMPAT
CONFIG_DIR	= CONFIG
CONFIG_DYNAMIC_FILES	= CONFIG_FILE.LBL
CONFIG_STATIC_FILES	= CFGINFO.TXT
DYNAMIC_FILES		=
STATIC_FILES		= CALINFO.TXT
TARGET_DIR	= CALIB

all:	tree

tree:	beampat config dynamic static

beampat:
	@if (test ! -d $(ANTENNA_PATTERN_DIR)); \
		then echo "Creating $(ANTENNA_PATTERN_DIR)"; \
		$(MKDIR) $(ANTENNA_PATTERN_DIR); \
	fi;
	@ cd $(ANTENNA_PATTERN_DIR); \
	for file in $(ANTENNA_PATTERN_STATIC_FILES); \
	do if (test ! -e $$file); \
		then echo "Installing $$file"; \
		co $(RAS_CENTRAL_TREE)/PDS/$(TARGET_DIR)/$(ANTENNA_PATTERN_DIR)/RCS/$$file,v; \
		fi; \
	done

config:
	@if (test ! -d $(CONFIG_DIR)); \
		then echo "Creating $(CONFIG_DIR)"; \
		$(MKDIR) $(CONFIG_DIR); \
	fi;
	@ cd $(CONFIG_DIR); \
	for file in $(CONFIG_DYNAMIC_FILES); \
	do if (test ! -e $$file); \
		then echo "Installing $$file"; \
		co $(RAS_CENTRAL_TREE)/PDS/$(TARGET_DIR)/$(CONFIG_DIR)/RCS/$$file,v; \
		$(CHMOD) 644 $$file; \
		fi; \
	done
	@ cd $(CONFIG_DIR); \
	for file in $(CONFIG_STATIC_FILES); \
	do if (test ! -e $$file); \
		then echo "Installing $$file"; \
		co $(RAS_CENTRAL_TREE)/PDS/$(TARGET_DIR)/$(CONFIG_DIR)/RCS/$$file,v; \
		fi; \
	done

dynamic:
	@ for file in $(DYNAMIC_FILES); \
	do if (test ! -e $$file); \
		then echo "Installing $$file"; \
		co $(RAS_CENTRAL_TREE)/PDS/$(TARGET_DIR)/RCS/$$file,v; \
		$(CHMOD) 644 $$file; \
		fi; \
	done

static:
	@ for file in $(STATIC_FILES); \
	do if (test ! -e $$file); \
		then echo "Installing $$file"; \
		co $(RAS_CENTRAL_TREE)/PDS/$(TARGET_DIR)/RCS/$$file,v; \
		fi; \
	done

clean:
	@ if (test -d $(ANTENNA_PATTERN_DIR)); \
		then echo Cleaning $(ANTENNA_PATTERN_DIR); \
		(cd $(ANTENNA_PATTERN_DIR); \
		for file in $(ANTENNA_PATTERN_STATIC_FILES); \
		do if (test -e $$file); \
			then echo "Deleting $$file"; \
			$(RM) $$file; \
			fi; \
		done); \
		echo "Deleting $(ANTENNA_PATTERN_DIR)"; \
		$(RMDIR) $(ANTENNA_PATTERN_DIR); \
	fi
	@ if (test -d $(CONFIG_DIR)); \
		then echo Cleaning $(CONFIG_DIR); \
		(cd $(CONFIG_DIR); \
		for file in $(CONFIG_DYNAMIC_FILES); \
		do if (test -e $$file); \
			then echo "Deleting $$file"; \
			$(RM) $$file; \
			fi; \
		done); \
		for file in $(CONFIG_STATIC_FILES); \
		do if (test -e $$file); \
			then echo "Deleting $$file"; \
			$(RM) $$file; \
			fi; \
		done); \
		echo "Deleting $(CONFIG_DIR)"; \
		$(RMDIR) $(CONFIG_DIR); \
	fi
	@ for file in $(DYNAMIC_FILES) $(STATIC_FILES); \
	do if (test -e $$file); \
		then echo "Deleting $$file"; \
		$(RM) $$file; \
		fi; \
	done
