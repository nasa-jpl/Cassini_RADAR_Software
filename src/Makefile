#======================================================================
# @(#) $Id: Makefile,v 11.5 2011/09/16 00:03:29 richw Exp $
# Makefile for all RAS source code
#======================================================================

#----------------------------------------------------------------------
# default: make all object and executable files
#----------------------------------------------------------------------

default: objs/Makefile programs/Makefile
	@ for dir in objs programs; \
		do (cd $$dir; \
			echo "Making default in `pwd`"; \
			make default); \
	done

#----------------------------------------------------------------------
# clean: remove all object, executable and  unlocked make files
#----------------------------------------------------------------------

clean:
	@ for dir in objs programs; \
		do (cd $$dir; \
			echo "Making clean in `pwd`"; \
			make clean ;); \
	done


#----------------------------------------------------------------------
# install: make all object and executable files and install
#----------------------------------------------------------------------

install: objs/Makefile programs/Makefile
	@ for dir in objs programs; \
		do (cd $$dir; \
			echo "Making install in `pwd`"; \
			make install); \
	done

tree:


