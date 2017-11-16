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
int NX,NY,STEPS,LOGGING,LOG2BIN;
double EPSILON,KAPPA,k_bT;

int main( int argc, char **argv )
{
  /*	Logging toggle to print to terminal	*/
  switch (argc){
   case 8:
       sscanf(argv[1],"%d",&NX);
       sscanf(argv[2],"%d",&NY);
       sscanf(argv[3],"%lf",&KAPPA);
       sscanf(argv[4],"%d",&STEPS);
       sscanf(argv[5],"%d",&RUNS);
       sscanf(argv[6],"%lf",&k_bT);
       sscanf(argv[7],"%d",&LOGGING);
       break;
   default:
       print_and_exit("Usage Pass command line arguments: NX NY KAPPA STEPS #RUNS kBT 0 or 1 for LOGGING\n");
   }
 
  /*	Epsilon	*/
  EPSILON = 720.0 * KAPPA;
 
  /*	Total MD Steps and Jack Knife Bin Count	*/
  DATA_COUNT = STEPS/PERIOD;
  JK_BIN_COUNT = RUNS;
 
  /*    Printing read paramater.dat data        */
  if (LOGGING == 1)
  {
	  printf("NX = %d\n",NX);
	  printf("NY = %d\n",NY);
	  printf("STEPS\t%d \n",STEPS);
	  printf("PERIOD\t%d \n",PERIOD);
	  printf("EPSILON\t%f \n",EPSILON);
	  printf("KAPPA\t%f \n",KAPPA);
	  printf("DATA_COUNT\t%d \n",DATA_COUNT);
	  printf("JK_BIN_COUNT\t%d \n",JK_BIN_COUNT);
	  printf("RUNS\t%d\n",RUNS);
  }

 log2_single_observable_time_evolution(DATA_COUNT,JK_BIN_COUNT);
 
 return 0;

}

