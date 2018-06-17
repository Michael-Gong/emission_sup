//////////////////////////////////////////////////////////////////////////////
///////// physical constants in CGS units
///////// mathematical constants
///////// precision controls
///////// designed by XuHan, in Nation University of Defence Techonology, 
///////// 2007-8-4
/////////////////////////////////////////////////////////////////////////////

#ifndef _MYCOMMON_H
#define _MYCOMMON_H
#include <math.h>

#define  Vc0       2.9979246e+10  // cm/s     light velocity
#define  Qe0       4.8065295e-10  // staticC  electron charge
#define  Me0       9.1093819e-28  // g        electron mass

#define  EPSILON   1.0e-8 
#define  MAXLENGTH 128            // maximum characters of file name

#ifndef PI
#define PI    3.1415926
#endif 
#ifndef M_PI
#define M_PI  3.1415926
#endif


#ifdef PARALLEL
#include <mpi.h>
#endif 

#ifdef LOWPRECISION
#define Scalar float
#define H5T_SCALAR H5T_NATIVE_FLOAT
#ifdef PARALLEL
#define MPI_SCALAR MPI_FLOAT
#endif
#else
#define Scalar double
#define H5T_SCALAR H5T_NATIVE_DOUBLE
#ifdef PARALLEL
#define MPI_SCALAR MPI_DOUBLE
#endif
#endif

int part_number;
int set;
//const double time0=0.0;
double a0;//=150.0;
double r0;//=5;
double tau;//=3.0;
double gamma0; //=1000.0
double px0; //=1000.0
double py0; //=1000.0
double pz0; //=1000.0
double total_t;//=100.0*2*M_PI;
double dt;//=2*M_PI/100.0;
double dump_dt; // 
double a0_out; // amplitude for self-generated field
double ratio_out; // ration for self-generated field
//int total_steps=10010;

#endif
