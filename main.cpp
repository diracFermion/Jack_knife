#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <math.h>
#include <string.h>
#include <time.h>
#include <limits>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <stdarg.h>
#include <unistd.h>
#include <errno.h>
#include "variables.h"
#include "jack_knife.h"


int DATA_COUNT,BIN_SIZE,JK_BIN_COUNT,RUNS;
char parameter_file[1024],result_file[1024];
FILE *Finput,*Fresult;
int NX,NY,STEPS,PERIOD,STRIP_SIZE,DATA_DISCARD,LOGGING,LINEARBIN,LOG2BIN;
double EPSILON,KAPPA,k_bT,DISP;

int main( int argc, char **argv )
{
  char param_file[256],dumpFile[256]; 


  //char *param_file = argv[2];
  sprintf(param_file,"parameter.dat");
  //printf("Parameter file : %s\n",param_file);

     

  /*	Logging toggle to print to terminal	*/
  switch (argc){
   case 2:
       LOGGING=atoi(argv[1]);
       break;
   default:
       print_and_exit("Usage Pass command line arguments 0 or 1 for LOGGING and change parameter.dat  \n");
   }

  /*	Reading the parameter file	*/
  Finput=fopen(param_file,"r");
  if (Finput==NULL)
      print_and_exit(" ERROR Could NOT open parameter.dat file.\n");
  
  fgets(parameter_file,1024,Finput);
  sscanf(parameter_file,"NX %d",&NX);
  fgets(parameter_file,1024,Finput);
  sscanf(parameter_file,"NY %d",&NY);
  fgets(parameter_file,1024,Finput);
  sscanf(parameter_file,"STEPS %d",&STEPS);
  fgets(parameter_file,1024,Finput);
  sscanf(parameter_file,"PERIOD %d",&PERIOD);
  fgets(parameter_file,1024,Finput);
  sscanf(parameter_file,"STRIPSIZE %d",&STRIP_SIZE);
  fgets(parameter_file,1024,Finput);
  sscanf(parameter_file,"EPSILON %lf",&EPSILON);
  fgets(parameter_file,1024,Finput);
  sscanf(parameter_file,"KAPPA %lf",&KAPPA);
  fgets(parameter_file,1024,Finput);
  sscanf(parameter_file,"k_BT %lf",&k_bT);
  fgets(parameter_file,1024,Finput);
  sscanf(parameter_file,"RUNS %d",&RUNS);
  fclose(Finput); 
  
  
  /*	Total MD Steps and Jack Knife Bin Count	*/
  DATA_COUNT = STEPS/PERIOD;
  //JK_BIN_COUNT = atoi(argv[1]);/*	 Number of runs (each run is a block)	*/
  JK_BIN_COUNT = RUNS;
 
  /*    Printing read paramater.dat data        */
  if (LOGGING == 1)
  {
	  printf("\nReading parameter.dat\n");
	  printf("parameter file NX = %d\n",NX);
	  printf("parameter file NY = %d\n",NY);
	  printf("strip size\t%d \n",STRIP_SIZE);
	  printf("STEPS\t%d \n",STEPS);
	  printf("PERIOD\t%d \n",PERIOD);
	  printf("EPSILON\t%f \n",EPSILON);
	  printf("KAPPA\t%f \n",KAPPA);
	  printf("DATA_COUNT\t%d \n",DATA_COUNT);
	  printf("JK_BIN_COUNT\t%d \n",JK_BIN_COUNT);
	  printf("RUNS\t%d\n",RUNS);
  }

 /*	Opening the output file for Simulation Average and Error        */
/* sprintf(result_file,"../Results/Average/avg_linear.dat");
 Fresult = fopen(result_file, "a");
 if (Fresult == NULL)
	{
	    print_and_exit("Could Not Open Result Output File: Results/Average/avg_linear.dat");
	}
 fprintf(Fresult, "%d\t%d\t%.1f\t%.1f\t%d\t%.4f\t%.1f",NX,NY,EPSILON,KAPPA,STRIP_SIZE,DISP,k_bT);
 fclose(Fresult);
*/
 /*	Error Analysis Function Calls	*/
// char dumpFile[]="observable.log";
 //sprintf(dumpFile,"../Sim_dump_ribbon/TE_L%d_W%d_k%.1f.log",NX,NY,KAPPA);
 //printf("Time Evolution of observables, dumpFile=%s\n",dumpFile);

 log2_single_observable_time_evolution(DATA_COUNT,JK_BIN_COUNT);
 
 return 0;

}

