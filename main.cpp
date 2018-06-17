#include <iostream>
#include <mpi.h>
#include <fstream>
#include <math.h>
#include <ctime>
#include <string>
#include <sstream>
#include <stdlib.h>
#include "common.h"
#include "ovector.h"
#include "givefield.h"
#include "readfile.h"
//#include "parameter.h"

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"



using namespace std;

template <typename T>
  string NumberToString ( T Number )
  {
     ostringstream ss;
     ss << Number;
     return ss.str();
  }

int main(int argc, char **argv)
{ 
    double startwtime = 0.0;
    double endwtime;
    double temp_factor, suppression_factor, a_random_num;  // for emission suppresion
 //    stringstream ss;
 //    string str;
 // 
 //    Vector3 dx1,dx2,dx3,dx4,dp1,dp2,dp3,dp4,constant0;
 //    int part_number=10; // the number of the electron 
 //    int set=1; //set = 1 means the radiation reaction is switched on.
 //  
 //    Scalar time0,a0,r0,tau, t_start, t_stop;
 //    time0  = 0.0;
 //    a0	   = 100.0;
 //    r0	   = 5.0;
 //    tau    = 3.0; 
 //    t_start= get_time();
 //
 //    Scalar total_t, dt;
 //    int total_steps;
 //    total_t= 100.0*2*M_PI;     
 //    dt     = 2*M_PI/100.0;
 //    total_steps = 10010;
 //    
 //    Particle p;
    //int part_number=10; // the number of the electron 
    //int set=2; //set = 0 means without RR; 1 means the classical radiation reaction; 2 means the QED recoil.
//###################################################
    readfile rf;
    char input_filename[100];
    if(argc==1) { // no argument
        strcpy(input_filename,"./Data/input.deck");
    }
    printf("%d%s%s%s\n",argc,argv[0],argv[1],argv[2]);

///// read parameters from input file

    rf.openinput(input_filename);

    set   = atoi( rf.setget("@common","set") );
    part_number = atoi( rf.setget("@common","part_number"));
    a0    = atof( rf.setget("@common","a0") )/100.0;
    r0    = atof( rf.setget("@common","r0") )/10.0;
    tau   = atof( rf.setget("@common","tau") )/10.0;
    px0= atof( rf.setget("@common","px0") );
    py0= atof( rf.setget("@common","py0") );
    pz0= atof( rf.setget("@common","pz0") );
    total_t = atof( rf.setget("@common","total_t") );
    dt    = atof( rf.setget("@common","dt") );
    dump_dt = atof( rf.setget("@common","dump_dt") );

    a0_out  = atof( rf.setget("@common","a0_out") )/10000.0;
    //a0_out  = pow(10.0, a0_out-3.5);
    ratio_out = atof( rf.setget("@common","ratio_out") )/1000.0;
//std::cout<<"r0="<<r0<<endl;
//#################################################

    ifstream ifile("./TABLES/hsokolov.table");
    if(! ifile){cerr<<"error."<<endl;return -1;}
     int h_sample;
     double h_min, h_max;
    ifile>>h_sample;
     double T1_1[500];
     double T1_2[500];
    ifile>>h_min; ifile>>h_max;
    for(int i=0;i<=h_sample-1;i++)
	{ 
	    ifile>>T1_1[i];
	    ifile>>T1_2[i];
	}	
    ifile.close();
    ifile.clear();
//###################################################
    ifile.open("./TABLES/ksi_sokolov.table");
    if(! ifile){cerr<<"error."<<endl; return -1;}
     int n_sample_eta;
     int n_sample_chi;
    ifile>>n_sample_eta; 
    ifile>>n_sample_chi; 
     double T2[100][100];
     double etalog_min, etalog_max;
    ifile>>etalog_min; 
    ifile>>etalog_max;
    for(int i=0;i<=n_sample_eta-1;i++)
	  { 
	    for(int j=0;j<=n_sample_chi-1;j++)
	      {
	  	ifile>>T2[i][j];
	      }
	  }	
    ifile.close();
    ifile.clear();
//###################################################
    ifile.open("./TABLES/chimin.table");
    if(! ifile){cerr<<"error."<<endl;return -1;}
     double T3[100];//n_sample_eta
    for(int i=0;i<=n_sample_eta-1;i++)
          { 
	    ifile>>T3[i];
	  }	
   ifile.close();
   ifile.clear(); 
//##################################################
//##################################################
    double log_eta[100]; //n_sample_eta  
    double log_chi[100][100]; //n_sample_eta, n_sample_chi
   double etalog_dx,chi_min,chi_dx;
   etalog_dx = (etalog_max-etalog_min)/((double)(n_sample_eta-1));
   for (int ieta=0; ieta<n_sample_eta; ieta++)
	{
		log_eta[ieta]=etalog_min+(double)(ieta-1)*etalog_dx;
		chi_min=log10(T3[ieta]);
		chi_dx=(log_eta[ieta]-log10(2.0)-chi_min)/((double)(n_sample_chi-1));
 		for(int ichi=0; ichi<n_sample_chi; ichi++)
		{  log_chi[ieta][ichi]=chi_min+double(ichi-1)*chi_dx;}          
 	}   



//##################################################
//##################################################
//    cout<<"T1_1 and T1_2!"<<endl;
//    for (int i=0; i<=499; i++)
//	{cout<<T1_1[i]<<" ; "<<T1_2[i]<<endl;}
//    cout<<"T2"<<endl;
//    for (int i=0; i<=99; i++)
//	{
//	 for (int j=0; j<=99; j++)
//	     {
//		cout<<T2[i][j]<<" ";
//	     }    
//	 cout<<endl;
//	}
//    cout<<"log_eta!"<<endl;
//    for (int i=0; i<=99; i++)
//	{cout<<log_eta[i]<<endl;}
//    cout<<"log_chi"<<endl;
//    for (int i=0; i<=99; i++)
//	{
//	 for (int j=0; j<=99; j++)
//	     {
//		cout<<log_chi[i][j]<<" ";
//	     }    
//	 cout<<endl;
//	}

//##################################################
//##################################################
//##################################################
//##################################################


    //Scalar total_t, dt;
    //int total_steps;
    //total_t= 100.0*2*M_PI;     
    //dt     = 2*M_PI/100.0;
    //total_steps = 10010;
    Scalar initial_x[part_number], initial_y[part_number], initial_z[part_number];
    //for(int in=0; in<part_number; in++)
    //	{
    //	  initial_x[in]=(95.0+(double)(in/part_number)*5.0)*2*M_PI;
    //	  initial_y[in]=(-2.50+(double)(in%part_number)/part_number.0*5.0)*2*M_PI;
    //	}
    int iin =0;
    double r0_t, y0, z0;
    while(iin < part_number)
        {    
            initial_y[iin]=double(iin)/double(part_number)*1.0*2*M_PI;
    	    initial_z[iin]=0.0;
            initial_x[iin]=100.0*2*M_PI; //(100.0+rand()/double(RAND_MAX)*4.0)*2*M_PI;
            iin++;
	        //y0 = (rand()/double(RAND_MAX)*5.0-2.5)*2*M_PI;
            //z0 = (rand()/double(RAND_MAX)*5.0-2.5)*2*M_PI;
           //r0_t = pow(pow(y0,2.0),0.5);
            //r0_t = pow(pow(y0,2.0)+pow(z0,2.0),0.5);
           //if ( (r0_t > 2.5*2*M_PI))
           //  { continue; }
           //else
           //  {
           //    initial_x[iin]=20.0*2*M_PI; //(100.0+rand()/double(RAND_MAX)*4.0)*2*M_PI;
 	       //    initial_y[iin]=y0+3.0*2*M_PI;
 	       //    initial_z[iin]=z0;
           //    iin++;
           //  }
        }


    if (set==0){ cout<<"This is no RR"<<endl;}
    if (set==1){ cout<<"This is rr RR"<<endl;}
    if (set==2){ cout<<"This is qe RR"<<endl;}
    Scalar t_start, t_stop;
    t_start= get_time();

    // Initialize the MPI environment
	MPI_Init(&argc, &argv);
	// Find out rank, size
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	
	if(world_rank==0)
	startwtime = MPI_Wtime();
    
	int in;
	int interval = part_number/world_size;
	if(part_number%world_size) interval++;
    for(in=interval*world_rank ; in<world_rank*(part_number/world_size)+interval; in++)
    {
	if(in>=part_number)
	       break;
    Scalar time0; //,a0,r0,tau, t_start, t_stop;
    time0  = 0.0;
    //a0	   = 100.0;
    //r0	   = 5.0;
    //tau    = 3.0; 
    stringstream ss;
    string str;
 
    Vector3 dx1,dx2,dx3,dx4,dp1,dp2,dp3,dp4,constant0;
  
    Particle p;
    Vector3 efield, bfield; 
       
       //   Scalar t_px[total_steps];
       //   Scalar t_py[total_steps];    
       //   Scalar t_x[total_steps];
       //   Scalar t_y[total_steps];
       //   Scalar t_t[total_steps];
      
          p.grid.x1=initial_x[in]; p.grid.x2=initial_y[in]; p.grid.x3=initial_z[in]; p.uu.x1=px0; p.uu.x2=py0; p.uu.x3=pz0;
          p.gamma=factor(p.uu); p.ID=in; 
#ifdef QED_BLOCK
	  p.eta=0.0;
 	  p.optical_depth=initial_optical_depth();
	  p.radn=0;
	  p.radt=0.0;
	  p.rad_px=0.0;
	  p.rad_py=0.0;
	  p.rad_pz=0.0;
	  Scalar photon_energy=0.0;
#endif
      //    t_t[0]=time0; t_x[0]=p.grid.x1; t_y[0]=p.grid.x2; t_px[0]=p.uu.x1; t_py[0]=p.uu.x2; 
         char filename[MAXLENGTH];
         sprintf(filename,"./Data/px_%d.txt",world_rank);
         FILE *fp_px=fopen(filename,"a");
         sprintf(filename,"./Data/py_%d.txt",world_rank);
         FILE *fp_py=fopen(filename,"a");
         sprintf(filename,"./Data/pz_%d.txt",world_rank);
         FILE *fp_pz=fopen(filename,"a");
         sprintf(filename,"./Data/x_%d.txt",world_rank);
         FILE *fp_x=fopen(filename,"a");
         sprintf(filename,"./Data/y_%d.txt",world_rank);
         FILE *fp_y=fopen(filename,"a");
         sprintf(filename,"./Data/z_%d.txt",world_rank);
         FILE *fp_z=fopen(filename,"a");
         sprintf(filename,"./Data/t_%d.txt",world_rank);
         FILE *fp_t=fopen(filename,"a");
         sprintf(filename,"./Data/e_part_%d.txt",world_rank);
         FILE *fp_e_part=fopen(filename,"a");
         sprintf(filename,"./Data/b_part_%d.txt",world_rank);
         FILE *fp_b_part=fopen(filename,"a");
         //sprintf(filename,"./Data/px_%d.txt",world_rank);
         //FILE *fp_a_part=fopen("./Data/a_part.txt","a");
         sprintf(filename,"./Data/eta_%d.txt",world_rank);
         FILE *fp_eta=fopen(filename,"a");
         sprintf(filename,"./Data/opt_%d.txt",world_rank);
         FILE *fp_opt=fopen(filename,"a");
         sprintf(filename,"./Data/radt_%d.txt",world_rank);
         FILE *fp_radt=fopen(filename,"a");
         sprintf(filename,"./Data/rad_px_%d.txt",world_rank);
         FILE *fp_rad_px=fopen(filename,"a");
         sprintf(filename,"./Data/rad_py_%d.txt",world_rank);
         FILE *fp_rad_py=fopen(filename,"a");
         sprintf(filename,"./Data/rad_pz_%d.txt",world_rank);
         FILE *fp_rad_pz=fopen(filename,"a");
         sprintf(filename,"./Data/radn_%d.txt",world_rank);
         FILE *fp_radn=fopen(filename,"a");


         sprintf(filename,"./Data/photon/no_s_px_%d.txt",world_rank); // photon emission suppression
         FILE *fp_no_s_px=fopen(filename,"a");
         sprintf(filename,"./Data/photon/no_s_py_%d.txt",world_rank);
         FILE *fp_no_s_py=fopen(filename,"a");
         sprintf(filename,"./Data/photon/no_s_pz_%d.txt",world_rank);
         FILE *fp_no_s_pz=fopen(filename,"a");
         sprintf(filename,"./Data/photon/try_s_px_%d.txt",world_rank); // photon emission suppression
         FILE *fp_try_s_px=fopen(filename,"a");
         sprintf(filename,"./Data/photon/try_s_py_%d.txt",world_rank);
         FILE *fp_try_s_py=fopen(filename,"a");
         sprintf(filename,"./Data/photon/try_s_pz_%d.txt",world_rank);
         FILE *fp_try_s_pz=fopen(filename,"a");
         sprintf(filename,"./Data/photon/s_px_%d.txt",world_rank); // photon emission suppression
         FILE *fp_s_px=fopen(filename,"a");
         sprintf(filename,"./Data/photon/s_py_%d.txt",world_rank);
         FILE *fp_s_py=fopen(filename,"a");
         sprintf(filename,"./Data/photon/s_pz_%d.txt",world_rank);
         FILE *fp_s_pz=fopen(filename,"a");


         //fprintf(fp_px,"id\n");// is %d\n",p.ID);
         //fprintf(fp_py,"id\n");// is %d\n",p.ID);
         //fprintf(fp_x,"id\n");// is %d\n",p.ID);
         //fprintf(fp_y,"id\n");// is %d\n",p.ID);
         //fprintf(fp_t,"id\n");// is %d\n",p.ID);
         //fprintf(fp_eta,"id\n");// is %d\n",p.ID);
         //fprintf(fp_opt,"id\n");// is %d\n",p.ID);
            fprintf(fp_px,"%lf\n",p.uu.x1);
            fprintf(fp_py,"%lf\n",p.uu.x2);
            fprintf(fp_pz,"%lf\n",p.uu.x3);
            fprintf(fp_x,"%lf\n",p.grid.x1);
            fprintf(fp_y,"%lf\n",p.grid.x2);
            fprintf(fp_z,"%lf\n",p.grid.x3);
            fprintf(fp_t,"%lf\n",time0);
            fprintf(fp_e_part,"%lf\n",p.e_part.x2);
            fprintf(fp_b_part,"%lf\n",p.b_part.x3);
            //fprintf(fp_a_part,"%lf\n",p.a_part.x2);
            fprintf(fp_eta,"%lf\n",p.eta);
            fprintf(fp_opt,"%lf\n",p.optical_depth);
            fprintf(fp_radn,"%d\n",p.radn);
            fprintf(fp_radt,"%lf\n",p.radt);
            fprintf(fp_rad_px,"%lf\n",p.rad_px);
            fprintf(fp_rad_py,"%lf\n",p.rad_py);
            fprintf(fp_rad_pz,"%lf\n",p.rad_pz);
//    std::cout<<"r0="<<r0<<endl;
          int dump_n=0;
          for(int it=0; time0<=total_t; it++)  // domain iteration 
          {// cout<<time0<<"\n";
           // cout<<dt<<"\n";
            efield = E_field(p.grid,time0,a0,tau,r0);
            bfield = B_field(p.grid,time0,a0,tau,r0)+out_B_field(p.grid,a0_out,ratio_out);
            dx1 = velocity(p.uu);
            dp1 = constant0-(efield+dx1.cross(bfield));
            if (set==1){ dp1=dp1-dx1*1.2e-8*pow(factor(p.uu),2)*(dp1*dp1-pow(efield*dx1,2));}
      
            efield = E_field(p.grid+dx1*0.5*dt,time0+0.5*dt,a0,tau,r0);
            bfield = B_field(p.grid+dx1*0.5*dt,time0+0.5*dt,a0,tau,r0)+out_B_field(p.grid+dx1*0.5*dt,a0_out,ratio_out);
            dx2 = velocity(p.uu+dp1*0.5*dt);
            dp2 = constant0-(efield+dx2.cross(bfield));
            if (set==1){ dp2=dp2-dx2*1.2e-8*pow(factor(p.uu+dp1*0.5*dt),2)*(dp2*dp2-pow(efield*dx2,2));}
      
      
            efield = E_field(p.grid+dx2*0.5*dt,time0+0.5*dt,a0,tau,r0);
            bfield = B_field(p.grid+dx2*0.5*dt,time0+0.5*dt,a0,tau,r0)+out_B_field(p.grid+dx2*0.5*dt,a0_out,ratio_out);
            dx3 = velocity(p.uu+dp2*0.5*dt);
            dp3 = constant0-(efield+dx3.cross(bfield));
            if (set==1){ dp3=dp3-dx3*1.2e-8*pow(factor(p.uu+dp2*0.5*dt),2)*(dp3*dp3-pow(efield*dx3,2));}
      
            efield = E_field(p.grid+dx3*dt,time0+1.0*dt,a0,tau,r0);
            bfield = B_field(p.grid+dx3*dt,time0+1.0*dt,a0,tau,r0)+out_B_field(p.grid+dx3*dt,a0_out,ratio_out);
            dx4 = velocity(p.uu+dp3*dt);
            dp4 = constant0-(efield+dx4.cross(bfield));
            if (set==1){ dp4=dp4-dx4*1.2e-8*pow(factor(p.uu+dp3*dt),2)*(dp4*dp4-pow(efield*dx4,2));}
      
            p.grid   = p.grid+(dx1+dx2*2+dx3*2+dx4)*dt/6.0;
            p.uu     = p.uu  +(dp1+dp2*2+dp3*2+dp4)*dt/6.0;
            p.gamma  = factor(p.uu);
            time0    = time0+dt;
            p.e_part = E_field(p.grid,time0,a0,tau,r0);
            p.b_part = B_field(p.grid,time0,a0,tau,r0)+out_B_field(p.grid,a0_out,ratio_out);
            //p.a_part = A_field(p.grid,time0,a0,tau,r0);
#ifdef QED_BLOCK            
	    if (set==2)
 	    {   //cout<<"calculation qed!"<<endl; 
		p.eta = calculate_eta(p);
	  	//cout<<"calculate_eta "<<p.eta<<endl;
	        p.optical_depth = p.optical_depth - delta_optical_depth(p.eta, p.gamma,h_sample,T1_1,T1_2);
	  	//cout<<"optical_depth is"<<p.optical_depth<<endl;
	 	if (p.optical_depth <= 0.0)
		  {  
                     //  cout<<"radiation !!!^o^"<<endl;
		     photon_energy=calculate_photon_energy(p.eta,p.gamma,n_sample_eta,n_sample_chi,log_eta,log_chi,T2);
                     temp_factor = photon_energy/p.eta/(p.gamma-photon_energy); // emission suppression
                     suppression_factor=pow(temp_factor,0.666667);
                     a_random_num = (rand()/(double)RAND_MAX);
                     if ( suppression_factor >= a_random_num )
                       {     
                             if ( suppression_factor >= 1.0 )
              		        {
              		           fprintf(stderr, ANSI_COLOR_BLUE  "Totally no suppression!            suppression_factor is %.3f;\t random number is %.3f;\t emitted photon is %.3f;\n" ANSI_COLOR_RESET , suppression_factor, a_random_num, photon_energy);
                                   fprintf(fp_no_s_px,"%lf\n",photon_energy*p.uu.x1/pow(p.gamma*p.gamma-1.0,0.5));
                                   fprintf(fp_no_s_py,"%lf\n",photon_energy*p.uu.x2/pow(p.gamma*p.gamma-1.0,0.5));
                                   fprintf(fp_no_s_pz,"%lf\n",photon_energy*p.uu.x3/pow(p.gamma*p.gamma-1.0,0.5));
                                }
                             else
              		        {
              		           fprintf(stderr, ANSI_COLOR_YELLOW"emiision try to suppress, but not! suppression_factor is %.3f;\t random number is %.3f;\t emitted photon is %.3f;\n" ANSI_COLOR_RESET , suppression_factor, a_random_num, photon_energy);
                                   fprintf(fp_try_s_px,"%lf\n",photon_energy*p.uu.x1/pow(p.gamma*p.gamma-1.0,0.5));
                                   fprintf(fp_try_s_py,"%lf\n",photon_energy*p.uu.x2/pow(p.gamma*p.gamma-1.0,0.5));
                                   fprintf(fp_try_s_pz,"%lf\n",photon_energy*p.uu.x3/pow(p.gamma*p.gamma-1.0,0.5));
                                }
	 		     p.radn=p.radn+1;
		      	     p.radt=p.radt+photon_energy;
			     p.rad_px=p.rad_px+photon_energy*p.uu.x1/pow(p.gamma*p.gamma-1.0,0.5);
			     p.rad_py=p.rad_py+photon_energy*p.uu.x2/pow(p.gamma*p.gamma-1.0,0.5);
		             p.rad_pz=p.rad_pz+photon_energy*p.uu.x3/pow(p.gamma*p.gamma-1.0,0.5);
                             p.uu  = p.uu*(p.gamma-photon_energy)/p.gamma;
                 	     p.optical_depth=initial_optical_depth();	
                       }
                     else 
                       {     
              		     fprintf(stderr, ANSI_COLOR_RED         "emiision is suppressed!            suppression_factor is %.3f;\t random number is %.3f;\t suppressed photon is %.3f;\n" ANSI_COLOR_RESET , suppression_factor, a_random_num, photon_energy);
                             fprintf(fp_s_px,"%lf\n",photon_energy*p.uu.x1/pow(p.gamma*p.gamma-1.0,0.5));
                             fprintf(fp_s_py,"%lf\n",photon_energy*p.uu.x2/pow(p.gamma*p.gamma-1.0,0.5));
                             fprintf(fp_s_pz,"%lf\n",photon_energy*p.uu.x3/pow(p.gamma*p.gamma-1.0,0.5));
                 	     p.optical_depth=initial_optical_depth();	
                       }
	          } 	      
	    } 
#endif
            //t_px[it+1]=p.uu.x1; t_py[it+1]=p.uu.x2; t_x[it+1]=p.grid.x1; t_y[it+1]=p.grid.x2; t_t[it+1]=time0;
            t_stop = get_time();
            //fprintf(stderr, "iteration %d steps, finished %.4f %s ! And realtime is %.3f s\n", it+1, 100.0*time0/total_t,"%", t_stop-t_start);
      //      cout<<t_stop<<endl;
	        if(time0>=dump_dt*(1+dump_n))
			{
            fprintf(fp_px,"%lf\n",p.uu.x1);
            fprintf(fp_py,"%lf\n",p.uu.x2);
			fprintf(fp_pz,"%lf\n",p.uu.x3);
            fprintf(fp_x,"%lf\n",p.grid.x1);
            fprintf(fp_y,"%lf\n",p.grid.x2);
            fprintf(fp_z,"%lf\n",p.grid.x3);
            fprintf(fp_t,"%lf\n",time0);
            fprintf(fp_e_part,"%lf\n",p.e_part.x2);
            fprintf(fp_b_part,"%lf\n",p.b_part.x3);
            //fprintf(fp_a_part,"%lf\n",p.a_part.x2);
            fprintf(fp_eta,"%lf\n",p.eta);
            fprintf(fp_opt,"%lf\n",p.optical_depth);
            fprintf(fp_radn,"%d\n",p.radn);
            fprintf(fp_radt,"%lf\n",p.radt);
            fprintf(fp_rad_px,"%lf\n",p.rad_px);
            fprintf(fp_rad_py,"%lf\n",p.rad_py);
            fprintf(fp_rad_pz,"%lf\n",p.rad_pz);
			dump_n++;
			}
          } // domain iteration
//             fprintf(fp_px,"%lf\n",p.uu.x1);
//             fprintf(fp_py,"%lf\n",p.uu.x2);
//             fprintf(fp_pz,"%lf\n",p.uu.x3);
//             fprintf(fp_x,"%lf\n",p.grid.x1);
//             fprintf(fp_y,"%lf\n",p.grid.x2);
//             fprintf(fp_z,"%lf\n",p.grid.x3);
//             fprintf(fp_t,"%lf\n",time0);
//             fprintf(fp_e_part,"%lf\n",p.e_part.x2);
//             fprintf(fp_b_part,"%lf\n",p.b_part.x3);
//             //fprintf(fp_a_part,"%lf\n",p.a_part.x2);
//             fprintf(fp_eta,"%lf\n",p.eta)dump_dt;
//             fprintf(fp_opt,"%lf\n",p.optical_depth);
//             fprintf(fp_radn,"%d\n",p.radn);
//             fprintf(fp_radt,"%lf\n",p.radt);
//             fprintf(fp_rad_px,"%lf\n",p.rad_px);
//             fprintf(fp_rad_py,"%lf\n",p.rad_py);
//             fprintf(fp_rad_pz,"%lf\n",p.rad_pz);
          fclose(fp_px);
          fclose(fp_py);
          fclose(fp_pz);
          fclose(fp_x);
          fclose(fp_y);
          fclose(fp_z);
          fclose(fp_t);
          fclose(fp_e_part);
          fclose(fp_b_part);
       //   fclose(fp_a_part);
          fclose(fp_eta);
          fclose(fp_opt);
          fclose(fp_radn);
          fclose(fp_radt);
          fclose(fp_rad_px);
          fclose(fp_rad_py);
          fclose(fp_rad_pz);

 
          fclose(fp_no_s_px);
          fclose(fp_no_s_py);
          fclose(fp_no_s_pz);
          fclose(fp_try_s_px);
          fclose(fp_try_s_py);
          fclose(fp_try_s_pz);
          fclose(fp_s_px);
          fclose(fp_s_py);
          fclose(fp_s_pz);


          t_stop = get_time();
          fprintf(stdout, "This is %dth electron of core %dth; finished %.4f %s ! And realtime is %.3f s\n", in, world_rank, 100.0*((double)in)/part_number,"%", t_stop-t_start);
    }// each particle
    if(world_rank==0)
	{
	    endwtime = MPI_Wtime();
	    printf("Time elapsed: %fs\n", (endwtime-startwtime));
	}			   
	MPI_Finalize();
    return 0;
}
   
   
   
   

