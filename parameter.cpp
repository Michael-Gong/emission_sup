///////////////////////////////////////////////////////////////////////////
#include "parameter.h"
using namespace std;

Parameter::Parameter(int argc, char **argv)
{

//// deal with the command line argument
    if(argc==1) { // no argument
	strcpy(input_filename,"./input/input.deck");
    }
    printf("%d%s%s%s\n",argc,argv[0],argv[1],argv[2]);

///// read parameters from input file

    rf.openinput(input_filename);

    a0    = atof( rf.setget("@common","a0") );
    r0    = atof( rf.setget("@common","r0") );
    tau   = atof( rf.setget("@common","tau") );
    total_t = atof( rf.setget("@common","total_t") );
    dt    = atof( rf.setget("@common","dt") );

}


////////////////////////////////////////////////////////////////////////////////////
Scalar Parameter::get_time()
{
#ifdef PARALLEL
    return MPI_Wtime();
#else
    timeb tm;
    ftime(&tm);
    return tm.time+0.001*tm.millitm;
#endif
}


////////////////////////////////////////////////////////////////////////////////////
