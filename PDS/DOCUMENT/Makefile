#
# $Id: Makefile,v 11.1 2011/09/16 00:03:37 richw Exp $
#
# Makefile for DOCUMENT directory
#

CHMOD	= /bin/chmod
CP	= /bin/cp -i
MKDIR	= /bin/mkdir
RM	= /bin/rm -f
RMDIR	= /bin/rmdir

DYNAMIC_FILES	=
STATIC_FILES	= \
	BIDRSIS.HTML \
	BIDRSIS.LBL \
	BODPSIS.HTML \
	BODPSIS.LBL \
	VOLSIS.HTML \
	VOLSIS.LBL \
	DOCINFO.TXT
GIF_DIR		= GIF/CURRENT
PDF_DIR		= PDF/CURRENT
PNG_DIR		= PNG/CURRENT
JPG_DIR		= JPG/CURRENT
SUBDIRS		=
TARGET_DIR	= DOCUMENT

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
	@ echo "Installing PDFs";
	@ $(CP) $(RAS_CENTRAL_TREE)/PDS/$(TARGET_DIR)/$(PDF_DIR)/*.PDF .
	@ $(CHMOD) 444 *.PDF
	@ echo "Installing PNGs";
	@ $(CP) $(RAS_CENTRAL_TREE)/PDS/$(TARGET_DIR)/$(PNG_DIR)/*.PNG .
	@ $(CHMOD) 444 *.PNG
	@ echo "Installing JPGs";
	@ $(CP) $(RAS_CENTRAL_TREE)/PDS/$(TARGET_DIR)/$(JPG_DIR)/*.JPG .
	@ $(CHMOD) 444 *.JPG

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
	@ for file in *.PDF; \
	do if (test -e $$file); \
		then echo "Deleting $$file"; \
		$(RM) $$file; \
		fi; \
	done
	@ for file in *.PNG; \
	do if (test -e $$file); \
		then echo "Deleting $$file"; \
		$(RM) $$file; \
		fi; \
	done
	@ for file in *.JPG; \
	do if (test -e $$file); \
		then echo "Deleting $$file"; \
		$(RM) $$file; \
		fi; \
	done
