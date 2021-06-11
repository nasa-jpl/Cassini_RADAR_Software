#
# $Id: Makefile,v 11.6 2015/09/28 19:18:07 cveerama Exp $
#
# Makefile for CATALOG directory
#

CHMOD	= /bin/chmod
MKDIR	= /bin/mkdir
RM	= /bin/rm -f
RMDIR	= /bin/rmdir

DYNAMIC_FILES	= \
	ABDRDS.CAT \
	ABDRSUMDS.CAT \
	BIDRDS.CAT \
	LBDRDS.CAT \
	SBDRDS.CAT \
        STDRDS.CAT

STATIC_FILES	= \
	CATINFO.TXT \
	DSMAP.CAT \
	INST.CAT \
	INSTHOST.CAT \
	MISSION.CAT \
	PERSON.CAT \
	PROJREF.CAT \
	REF.CAT \
	TARGETALBIORIX.CAT \
	TARGETBESTLA.CAT \
	TARGETERRIAPUS.CAT \
	TARGETKIVIUQ.CAT \
	TARGETPAALIAQ.CAT \
	TARGETSIARNAQ.CAT \
	TARGETSKATHI.CAT \
	TARGETTARVOS.CAT \
	TARGETYMIR.CAT \
	TARGETHYROKKIN.CAT \
	SOFTWARE.CAT

TARGET_DIR	= CATALOG

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
