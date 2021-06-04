#======================================================================
# @(#) $Id: Makefile,v 11.5 2011/09/16 00:03:29 richw Exp $
# Top level Makefile for all RAS software
#======================================================================

#----------------------------------------------------------------------
# default: make all object and executable files
#----------------------------------------------------------------------

default: src/Makefile
	@ for dir in src; \
		do (cd $$dir; \
			echo "Making default in `pwd`"; \
			make default); \
	done

#----------------------------------------------------------------------
# tree: make a level of the tree
#----------------------------------------------------------------------

tree:
	@ for dir in doc include lib scripts src src/objs \
	src/programs; \
	do if (test ! -d $$dir); \
			then echo "Creating $$dir"; \
			mkdir $$dir; \
		fi; \
	done

#----------------------------------------------------------------------
# clean: remove all object and executable files
#----------------------------------------------------------------------

clean: src/Makefile
	@ /bin/rm -f include/{*.h,*.cpp}
	@ /bin/rm -f lib/*
	@ for dir in src; \
		do (cd $$dir; \
			echo "Making clean in `pwd`"; \
			make clean;); \
	done

#----------------------------------------------------------------------
# install: make all object and executable files and install
#----------------------------------------------------------------------

install: src/Makefile
	@ for dir in src; \
		do (cd $$dir; \
			echo "Making install in `pwd`"; \
			make install); \
	done

