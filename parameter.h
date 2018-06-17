//////////////////////////////////////////////////////////////////////////////////////////
//////   header file, input parameters, partition domain
//////   Oringinal structed by Xu Han at 2007-8  
//////   National University of Defence Technology, ChangSha, China
//////   Thanks to Zuo HongBin's for  useful discussion
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef _PARAMETERS_H
#define _PARAMETERS_H

#ifdef PARALLEL
#include "mpi.h"
#endif

#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <sys/timeb.h>
#include "readfile.h"
#include "common.h"
//#include "oerror.h"

#ifdef PARALLEL
#define IGNORE MPI_PROC_NULL
#else
#define IGNORE -1
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif 

class Parameter
{
public:
    readfile rf;
    char input_filename[MAXLENGTH];    
    Parameter(int argc, char **argv);
    
	~Parameter();

    Scalar get_time();  // static calculating time.
};

#endif
      
