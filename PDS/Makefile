#======================================================================
# @(#) $Id: Makefile,v 11.1 2011/09/16 00:03:36 richw Exp $
# Top level Makefile for PDS tree
#======================================================================

CHMOD	= /bin/chmod
MKDIR	= /bin/mkdir
RM	= /bin/rm -f
RMDIR	= /bin/rmdir

DYNAMIC_FILES	= \
	AAREADME.TXT \
	ERRATA.TXT \
	VOLDESC.CAT
STATIC_FILES	=

SUBDIRS	= \
	CALIB \
	CATALOG \
	DATA \
	DOCUMENT \
	EXTRAS \
	INDEX \
	SOFTWARE

all:	tree

#----------------------------------------------------------------------
# tree: make the tree
#----------------------------------------------------------------------

tree:
	@ for dir in $(SUBDIRS); \
	do if (test ! -d $$dir); \
			then echo "Creating $$dir"; \
			$(MKDIR) $$dir; \
		fi; \
		(cd $$dir; \
		co $(RAS_CENTRAL_TREE)/PDS/$$dir/RCS/Makefile,v; \
		make tree; \
		$(RM) Makefile); \
	done
	@ echo "Installing $(DYNAMIC_FILES) $(STATIC_FILES)"
	@ for file in $(DYNAMIC_FILES) $(STATIC_FILES); \
	do \
		co $(RAS_CENTRAL_TREE)/PDS/RCS/$$file,v; \
	done
	@ for file in $(DYNAMIC_FILES); \
	do \
		$(CHMOD) 644 $$file; \
	done

#----------------------------------------------------------------------
# clean: remove all files
#----------------------------------------------------------------------
clean:
	@ for dir in $(SUBDIRS); \
	do if (test -d $$dir); \
		  then (cd $$dir; \
		  echo "Making clean in `pwd`"; \
		  co $(RAS_CENTRAL_TREE)/PDS/$$dir/RCS/Makefile,v; \
		  make clean; \
		  $(RM) Makefile; \
		  cd ..; \
		  $(RMDIR) $$dir); \
		fi; \
	done
	@ for file in $(DYNAMIC_FILES) $(STATIC_FILES); \
	do if (test -e $$file); \
		then echo "Deleting $$file"; \
		$(RM) $$file; \
		fi; \
	done

# rule for subtree makefiles
%/Makefile: $(RAS_CENTRAL_TREE)/PDS/%/RCS/Makefile,v
	(cd $*; co $(RAS_CENTRAL_TREE)/PDS/$*/RCS/Makefile,v)
