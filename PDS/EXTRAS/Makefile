#
# $Id: Makefile,v 11.1 2011/09/16 00:03:37 richw Exp $
#
# Makefile for EXTRAS directory
#

CHMOD	= /bin/chmod
MKDIR	= /bin/mkdir
RM	= /bin/rm -f
RMDIR	= /bin/rmdir

DYNAMIC_FILES	= \
	EXTRINFO.TXT
STATIC_FILES	=
SUBDIRS		=
TARGET_DIR	= EXTRAS

all:	tree

tree:	dirs files

dirs:
	@ for dir in $(SUBDIRS); \
	do if (test ! -e $$dir); \
			then echo "Creating $$dir"; \
		$(MKDIR) $$dir; \
		fi; \
	done

files:
	@ for file in $(DYNAMIC_FILES) $(STATIC_FILES); \
	do if (test ! -e $$file); \
			then echo "Installing $$file"; \
		co $(RAS_CENTRAL_TREE)/PDS/$(TARGET_DIR)/RCS/$$file,v; \
		fi; \
	done
	@ for file in $(DYNAMIC_FILES); \
	do \
		$(CHMOD) 644 $$file; \
	done

clean:
	@ for dir in $(SUBDIRS); \
	do if (test ! -e $$dir); \
			then $(RMDIR) $$dir; \
		fi; \
	done
	@ for file in $(DYNAMIC_FILES) $(STATIC_FILES); \
	do if (test -e $$file); \
			then echo "Deleting $$file"; \
			$(RM) $$file; \
		fi; \
	done
