#
# $Id: Makefile,v 11.1 2011/09/16 00:03:37 richw Exp $
#
# Makefile for SOFTWARE directory
#

CHMOD	= /bin/chmod
CP	= /bin/cp -i
MKDIR	= /bin/mkdir
RM	= /bin/rm -f
RMDIR	= /bin/rmdir

DYNAMIC_FILES	=
STATIC_FILES	= \
	SOFTINFO.TXT
SUBDIRS		=
TARGET_DIR	= SOFTWARE

all:	tree

tree:	static dynamic

dynamic:
	@ for file in $(DYNAMIC_FILES); \
	do if (test ! -e $$file); \
		then echo "Installing $$file"; \
		co $(RAS_CENTRAL_TREE)/PDS/$(TARGET_DIR)/RCS/$$file,v; \
		$(CHMOD) 644 $$file; \
		fi; \
	done
	@ for dir in $(SUBDIRS); \
	do if (test ! -d $$dir); \
		then echo "Creating $$dir"; \
		$(MKDIR) $$dir; \
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
	@ for file in $(DYNAMIC_FILES) $(STATIC_FILES); \
	do if (test -e $$file); \
		then echo "Deleting $$file"; \
		$(RM) $$file; \
		fi; \
	done
	@ for dir in $(SUBDIRS); \
	do if (test -e $$dir); \
		then echo "Deleting $$dir"; \
		$(RMDIR) $$dir; \
		fi; \
	done
