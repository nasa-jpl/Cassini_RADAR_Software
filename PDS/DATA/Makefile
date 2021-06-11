#
# $Id: Makefile,v 11.3 2015/09/24 01:42:33 cveerama Exp $
#
# Makefile for DATA directory
#

CHMOD	= /bin/chmod
MKDIR	= /bin/mkdir
RM	= /bin/rm -f
RMDIR	= /bin/rmdir

DYNAMIC_FILES	= \
	ABDR_SUMMARY.FMT \
        STDR.FMT

STATIC_FILES	= \
	ABDR.FMT \
	LBDR.FMT \
	SBDR.FMT

TARGET_DIR	= DATA

all:	tree

tree:
	@ for file in $(DYNAMIC_FILES) $(STATIC_FILES); \
	do if (test ! -e $$file); \
		then echo "Installing $$file"; \
		co $(RAS_CENTRAL_TREE)/PDS/$(TARGET_DIR)/RCS/$$file,v; \
		fi; \
	done;
	@ for file in $(DYNAMIC_FILES); \
	do \
		$(CHMOD) 644 $$file; \
	done;

clean:
	@ for file in $(DYNAMIC_FILES) $(STATIC_FILES); \
	do if (test -e $$file); \
		then echo "Deleting $$file"; \
		$(RM) $$file; \
		fi; \
	done
