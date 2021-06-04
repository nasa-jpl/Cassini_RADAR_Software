//==============================================================//
// Copyright (C) 1997-2002, California Institute of Technology.    //
// U.S. Government sponsorship acknowledged.                    //
//==============================================================//

static const char rcs_id_array_c[] =
    "@(#) $Id: SimpleArray.cpp,v 11.5 2011/09/16 00:03:30 richw Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "SimpleArray.h"

//----------
// This module contains functions to allocate and free multidimensional
// arrays stored as a tree of pointers with the actual data stored in
// vectors at the fringe of the tree. These arrays are accessed by
// dereferencing the pointer. For example: a[3][2][5][1] for a
// 4-dimensional array.  No bounds checking is provided, and the
// pointer returned by make_array needs to be cast into the proper
// form. (eg., float**** for a 4 dimensional array of floats.)
// Following is a usage description for each function.
//
// make_array(int type_size, int ndims, int dimsize1, int dimsize2, ...)
//    type_size - the number of bytes (returned by sizeof()) for the desired
//        data type.
//    ndims - the number of dimensions desired.
//    dimsize1 - the size of the first dimension.
//    dimsizen - the size of the nth dimension.
//        A variable length argument list is used to allow an arbitrary number
//        of dimensions.
//
// free_array(void* ptr, int ndims, int dimsize1, int dimsize2, ...)
//    ptr - the pointer returned by a preceeding call to make_array.
//    ndims - the number of dimensions desired.
//    dimsize1 - the size of the first dimension.
//    dimsizen - the size of the nth dimension.
//        The arguments supplied to free_array must be consistent with the
//        preceeding call to make_array, otherwise memory leaks and/or damage
//        to the system memory heap will occur.
//----------

//------------//
// make_array //
//------------//

void*
make_array(
    int  type_size,
    int  ndims,
    ...)
{
    va_list ap;
    va_start(ap, ndims);

    // The first argument is the size (in bytes) of the data type desired.
    if (type_size <= 0)
    {
        va_end(ap);
        return(NULL);
    }

    // The second argument is the number of dimensions.
    if (ndims <= 0)
    {
        va_end(ap);
        return(NULL);
    }

    // Make an array of dimension sizes, and read from the argument list.
    int* dimsize = (int*)malloc((size_t)(ndims * sizeof(int)));
    if (dimsize == NULL)
    {
        va_end(ap);
        return(NULL);
    }
    for (int i=0; i < ndims; i++)
    {
        dimsize[i] = va_arg(ap, int);
    }

    // Initial call in the recursive allocation tree
    void* ptr = dim_setup(0, ndims, dimsize, type_size);

    free(dimsize);
    va_end(ap);
    return(ptr);
}

int    
write_array(FILE* ofp, void* ptr, int type_size, int ndims, ...){
      va_list ap;
    va_start(ap, ndims);

    // The first argument is the size (in bytes) of the data type desired.
    if (type_size <= 0)
    {
        va_end(ap);
        return(0);
    }

    // The second argument is the number of dimensions.
    if (ndims <= 0)
    {
        va_end(ap);
        return(0);
    }

    // Make an array of dimension sizes, and read from the argument list.
    int* dimsize = (int*)malloc((size_t)(ndims * sizeof(int)));
    if (dimsize == NULL)
    {
        va_end(ap);
        return(0);
    }
    for (int i=0; i < ndims; i++)
    {
        dimsize[i] = va_arg(ap, int);
    }

   int  retval=dim_write(ofp,ptr,0,ndims,dimsize,type_size);   
    free(dimsize);
    va_end(ap);
    return(retval);
}


int    dim_write(FILE* ofp, void* ptr, int level, int ndims, int dimsize[], int type_size){
  // bottom level
  if(level==ndims-1){
    return(fwrite(ptr,1,type_size*dimsize[level],ofp)==(unsigned)(type_size*dimsize[level]));    
  }
  else{
    void** subptr= (void**)ptr;
    for(int c=0;c<dimsize[level];c++){
      if(!dim_write(ofp,subptr[c],level+1,ndims,dimsize,type_size)){
	return(0);
      }
    }
  }
  return(1);
}

int    
read_array(FILE* ifp, void* ptr, int type_size, int ndims, ...){
      va_list ap;
    va_start(ap, ndims);

    // The first argument is the size (in bytes) of the data type desired.
    if (type_size <= 0)
    {
        va_end(ap);
        return(0);
    }

    // The second argument is the number of dimensions.
    if (ndims <= 0)
    {
        va_end(ap);
        return(0);
    }

    // Make an array of dimension sizes, and read from the argument list.
    int* dimsize = (int*)malloc((size_t)(ndims * sizeof(int)));
    if (dimsize == NULL)
    {
        va_end(ap);
        return(0);
    }
    for (int i=0; i < ndims; i++)
    {
        dimsize[i] = va_arg(ap, int);
    }

   int  retval=dim_read(ifp,ptr,0,ndims,dimsize,type_size);   
    free(dimsize);
    va_end(ap);
    return(retval);
}


int    dim_read(FILE* ifp, void* ptr, int level, int ndims, int dimsize[], int type_size){
  // bottom level
  if(level==ndims-1){
    return(fread(ptr,1,type_size*dimsize[level],ifp)==(unsigned)(type_size*dimsize[level]));    
  }
  else{
    void** subptr= (void**)ptr;
    for(int c=0;c<dimsize[level];c++){
      if(!dim_read(ifp,subptr[c],level+1,ndims,dimsize,type_size)){
	return(0);
      }
    }
  }
  return(1);
}
//------------//
// size_array //
//------------//

// computes the total amount of memory allocated by a make_array command
// including pointers and data
unsigned int
size_array(
    int  type_size,
    int  ndims,
    ...)
{
    va_list ap;
    va_start(ap, ndims);

    // The first argument is the size (in bytes) of the data type desired.
    if (type_size <= 0)
    {
        va_end(ap);
        return(0);
    }

    // The second argument is the number of dimensions.
    if (ndims <= 0)
    {
        va_end(ap);
        return(0);
    }

    // Make an array of dimension sizes, and read from the argument list.
    int* dimsize = (int*)malloc((size_t)(ndims * sizeof(int)));
    if (dimsize == NULL)
    {
        va_end(ap);
        fprintf(stderr,"Internal error in size_array()");
        exit(1);
    }

    unsigned int retval=0;
    for (int i=0; i < ndims; i++)
    {
        dimsize[i] = va_arg(ap, int);
        unsigned int size_this_level=1;
        for (int j=i; j>=0; j--){
	  size_this_level*=dimsize[j];
	}
	if(i==ndims-1) size_this_level*=type_size;
	else size_this_level*=sizeof(void*);
	retval+=size_this_level;
    }

    free(dimsize);
    va_end(ap);
    return(retval);
}

//------------//
// free_array //
//------------//

void
free_array(
    void*  ptr,
    int    ndims,
    ...)
{
    va_list ap;
    va_start(ap, ndims);

    // The first argument is the pointer to the array to be deallocated.
    if (ptr == NULL)
    {
        va_end(ap);
        return;
    }

    // The second argument is the number of dimensions.
    if (ndims <= 0)
    {
        va_end(ap);
        return;
    }

    // Make an array of dimension sizes, and read from the argument list.
    int* dimsize = (int*)malloc((size_t)(ndims * sizeof(int)));
    if (dimsize == NULL)
    {
        va_end(ap);
        return;
    }
    for (int i = 0; i < ndims; i++)
    {
        dimsize[i] = va_arg(ap, int);
    }

    // Initial call in the recursive deallocation tree
    dim_free(ptr, 0, ndims, dimsize);

    free(dimsize);
    va_end(ap);
    return;
}

//-----------//
// dim_setup //
//-----------//
// the recursive function that performs the allocation.
//
// level - the current dimension level (0 for the 1st dimension).
// ndims - the total number of dimensions.
// dimsize - a vector of dimension sizes (length = ndims).
// type_size - number of bytes in each data element.

void*
dim_setup(
    int  level,
    int  ndims,
    int  dimsize[],
    int  type_size)
{
    void** sub_ptr;

    if (level == ndims - 1)
    {    // reached lowest level, so allocate the actual data vector
        sub_ptr = (void**)malloc((size_t)(dimsize[level] * type_size));
        return((void*)sub_ptr);
    }

    // Intermediate levels require a vector of pointers each of which need
    // to be set by a call to dim_setup for the lower levels.
    sub_ptr = (void**)malloc((size_t)(dimsize[level] * sizeof(void*)));
    if (sub_ptr == NULL)
        return(NULL);

    for (int i = 0; i < dimsize[level]; i++)
    {
        // set each pointer to point to all the lower levels.
        sub_ptr[i] = dim_setup(level + 1, ndims, dimsize, type_size);
        if (sub_ptr[i] == NULL)
        {
            for (int j = i - 1; j >= 0; j--)
            {
                // free the pointers and data already allocated
                dim_free(sub_ptr[j], level + 1, ndims, dimsize);
            }
            return(NULL);
        }
    }

    return((void*)sub_ptr);
}

//----------//
// dim_free //
//----------//
// the recursive function that performs deallocation.
//
// ptr - pointer to the current subtree to be free'd.
// level - the current dimension level (0 for the 1st dimension).
// ndims - the total number of dimensions.
// dimsize - a vector of dimension sizes (length = ndims).

void
dim_free(
    void*  ptr,
    int    level,
    int    ndims,
    int    dimsize[])
{
    if (ptr == NULL)
        return;

    if (level == ndims - 1)
    {
        // reached lowest level, so deallocate the actual data vector
        free(ptr);
        return;
    }

    //
    // Intermediate levels require freeing of lower levels before
    // freeing the vector of pointers at this level.
    //

    for (int i = 0; i < dimsize[level]; i++)
    {
        dim_free(((void**)ptr)[i], level + 1, ndims, dimsize);
    }

    // Now free the vector of pointers at this level.
    free(ptr);

    return;
}
