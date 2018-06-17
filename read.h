#include<iostream>
#include<fstream>
#include<string>
#include<iomanip>

#ifdef QED_BLOCK
class Table1
{
  public:
    int 	h_sample;
    Scalar     	min;
    Scalar	max;
    Scalar      x1[500];
    Scalar	x2[500];
   
    Table1(Table1 * _t1){
	h_sample  = _t1->h_sample;
	min	  = _t1->min;
	max	  = _t1->max;
	x1	  = _t1->x1;
	x2	  = _t1->x2;
};  

  Table1(){
  };
};
 
class Table2
{
  public:
    int		n_sample_eta;
    int		n_sample_chi;
    Scalar	etalog_min;
    Scalar	etalog_max;
    Scalar	x[100][100];

  Table2(){
  };  
};

