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
#include "variables.h"
#include "jack_knife.h"


int raw_data_marker,log_bin_size;
unsigned long long raw_step[MAXMES];
char observable_file[1024],time_evolution_dir[1024],time_evolution_file[1024];
double raw_data[MAXRUNS][MAXMES],jk_blocks[MAXLOG][MAXRUNS+1];
double sum_log2_bin;
double log2_bin[MAXRUNS][MAXLOG],step_log2[MAXLOG];
double jk_avg[5][MAXLOG];
double error[5][MAXLOG],error_term1[5][MAXLOG],error_term2[5][MAXLOG];
char read_line[1024];


 /*	Print & Exit function	*/
 void print_and_exit(char *format, ...)
 {
    va_list list;
    va_start(list,format);
    vfprintf(stderr,format,list);
    va_end(list);
    exit(1);
 }

 /*	LOG2 BINNING TIME EVOLUTION	*/

int log2_single_observable_time_evolution(int data_size,int run_size,char dumpFile[])
{
        FILE *Finput,*Foutput;
	printf("HOOMD dumpFile inside function=%s\n",dumpFile);
	printf("Data_Size=%d,run_size=%d\n",data_size,run_size);
        int RUN,log_bin_cnt;
	double raw_observable[5];
	unsigned long long read_step;
        log_bin_cnt=log(double(data_size))/log(2);

   /*      Creating 2D array with raw observable data */
   for(int iobser=0;iobser<5;iobser++)
   {     
	for(RUN=0;RUN<JK_BIN_COUNT;RUN++)
        {
                sprintf(observable_file,"../Sim_dump_frame/Frame_%d/E_%.1f_K_%.2f_S_%d_R_%d/%s",NX,EPSILON,KAPPA,STRIP_SIZE,RUN+1,dumpFile);
                //printf("Observable Input File:%s\n",observable_file);
		Finput=fopen(observable_file,"r");
                if (Finput==NULL)
                        print_and_exit("**************Could NOT open HOOMD Observable.log file:%s*****************/n",observable_file);
                fgets(read_line, 1024,Finput); //Skip first line
                for(int j=0;j<data_size;j++)
                {
			fscanf(Finput,"%llu\t%lf\t%lf\t%lf\t%lf\t%lf\n",&read_step,raw_observable,raw_observable+1,raw_observable+2,raw_observable+3,raw_observable+4);
			raw_data[RUN][j]=raw_observable[iobser];
			raw_step[j]=read_step;
		}
	}

	/*	Geometric Average of Time Steps		*/
        raw_data_marker = data_size;	
	for(int ilog=0;ilog<log_bin_cnt;ilog++)
 	{
		for(int j=DATA_COUNT/int(pow(2,ilog+1));j<DATA_COUNT/int(pow(2,ilog));j++)
		{
			if (j==DATA_COUNT/int(pow(2,ilog+1)))
                                {
                                        step_log2[ilog]=raw_step[j];
                                        //printf("step_log2[%d]=%f\n",ilog,step_log2[ilog]);
                                }
                        if (j==(DATA_COUNT/int(pow(2,ilog))-1))
                                {
                                        step_log2[ilog]=sqrt(step_log2[ilog] * double(raw_step[j]));
                                        //printf("step_log2[%d]=%f\n",ilog,step_log2[ilog]);
                                        //printf("raw_step[%d]=%llu\n",j,raw_step[j]);
                                }	
		}
	}
	//step_log2[0]=0;
	//step_log2[1]=0;
	
	/*	Creating Log2 binned data	*/
	for(RUN=0;RUN<JK_BIN_COUNT;RUN++)
	{
		for(int ilog=0;ilog<log_bin_cnt;ilog++)
		{
			/*	Geometric Average of Time Steps         */
			//step_log2[ilog] = sqrt(double(DATA_COUNT/int(pow(2,ilog+1)))*double(DATA_COUNT/int(pow(2,ilog))));
			/*	Initializing Log2 bins array	*/
			log2_bin[RUN][ilog]=0;
			log_bin_size=0;
			for(int j=DATA_COUNT/int(pow(2,ilog+1));j<DATA_COUNT/int(pow(2,ilog));j++)
			{
				log2_bin[RUN][ilog]+=raw_data[RUN][j];
				log_bin_size++;
			/*	if (j==DATA_COUNT/int(pow(2,ilog+1)))
				{	
					step_log2[ilog]=raw_step[j];
					printf("step_log2[%d]=%llu\n",ilog,step_log2[ilog]);
				}
				if (j-1==DATA_COUNT/int(pow(2,ilog)))
				{	
					step_log2[ilog]=sqrt(step_log2[ilog] * raw_step[j]);
					printf("step_log2[%d]=%llu\n",ilog,step_log2[ilog]);
					printf("raw_step[%d]=%llu\n",j,raw_step[j]);
				}
			*/			
				//printf("RUN %d\tLog2_Bin[%d]\t%.8f\n",RUN,i-1,raw_data[RUN][raw_data_marker - j]);
			}
			log2_bin[RUN][ilog]=log2_bin[RUN][ilog]/double(log_bin_size);
			if(LOGGING == 1)
                        {
                                //printf("Log2_Bin[%d]\tBin_Size = %f\n",i,(pow(2,i)-pow(2,i-1)));
				printf("At Observable %d RUN %d\tLog2_Bin[%d][%d] = %.8f\n",iobser,RUN+1,RUN,ilog,log2_bin[RUN][ilog]);
			}
		}
	} 			

        /*      Total of Jack Knife Blocks at each Log2 interval     */
        for(int i=0;i<log_bin_cnt;i++)
        {
                jk_blocks[i][JK_BIN_COUNT]=0;
                for(int j=0;j<JK_BIN_COUNT;j++)
                {
                      jk_blocks[i][JK_BIN_COUNT] += log2_bin[j][i];
                }
                if(LOGGING==1)
                {
                        printf("At Observable %d Total JK at Log2 time bin[%d] = %.8f\n",iobser,i,jk_blocks[i][JK_BIN_COUNT]);
                }
		/*	Average of JK Blocks	*/
		jk_avg[iobser][i]=jk_blocks[i][JK_BIN_COUNT]/JK_BIN_COUNT;
        }			

        /*      Jack Knife Blocking at each Log2 time bin    */
        for(int i=0;i<log_bin_cnt;i++)
        {
		for(int j=0;j<JK_BIN_COUNT;j++)
                {
                	jk_blocks[i][j]=(jk_blocks[i][JK_BIN_COUNT]-log2_bin[j][i])/double(JK_BIN_COUNT-1);
                        if (LOGGING == 1)
                        {       
                                //printf("jk_blocks[%d][%d]=(total_jk_log2[%d]-log2_bin[%d][%d])/double(%d-1)\n",i,j,i,i,j,JK_BIN_COUNT);
                                //printf("%.8f\t%.8f\t%.8f\n",jk_blocks[i][j],total_jk_log2[i],log2_bin[j][i]);
                                printf("Jack Knife Block %d at time bin %f = %.8f\n",j+1,step_log2[i],jk_blocks[i][j]);
                        }
                }
	}

        /*      Master Jack Knife Average at each Log2 time bin      */
/*        for(int i=0;i<log_bin_cnt;i++)
        {
                jk_blocks[i][JK_BIN_COUNT]=0;
                for(int j=0;j<JK_BIN_COUNT;j++)
                {
                        jk_blocks[i][JK_BIN_COUNT]=jk_blocks[i][JK_BIN_COUNT]+jk_blocks[i][j];
                }
*/                /*      Jack Knife Average      */
/*                jk_blocks[i][JK_BIN_COUNT]=jk_blocks[i][JK_BIN_COUNT]/JK_BIN_COUNT;
                if (LOGGING == 1)
                {
                        printf("Jack Knife Average at Log2 time bin %f = %.8f\n",step_log2[i],jk_blocks[i][JK_BIN_COUNT]);
                }
        }

*/        /*      Error Calculation at each time bin      */
        for(int i=0;i<log_bin_cnt;i++)
        {
                error[iobser][i]=0;
                error_term1[iobser][i]=0;
                error_term2[iobser][i]=0;
        }

        /*      Term 1 in (square) error        */
        for(int i=0;i<log_bin_cnt;i++)
        {
                for(int j=0;j<JK_BIN_COUNT;j++)
                {
                        error_term1[iobser][i]=error_term1[iobser][i] + jk_blocks[i][j] * jk_blocks[i][j];
                }
                error_term1[iobser][i]=1/double(JK_BIN_COUNT) * error_term1[iobser][i];
                if(LOGGING == 1)
                {
                        printf("error_term1[%d][%d] = %.8f\n",iobser,i,error_term1[iobser][i]);
                }
        }
        /*      Term 2 in (square) error        */
        for(int i=0;i<log_bin_cnt;i++)
        {
                for(int j=0;j<JK_BIN_COUNT;j++)
                {
                        error_term2[iobser][i]=error_term2[iobser][i] + jk_blocks[i][j];
                }
                error_term2[iobser][i]=(1/double(JK_BIN_COUNT) * error_term2[iobser][i]) *(1/double(JK_BIN_COUNT) * error_term2[iobser][i]);
                if(LOGGING == 1)
                {
                        printf("error_term2[%d][%d] = %.8f\n",iobser,i,error_term2[iobser][i]);
                }
                /*      Error at each time bin  */
                error[iobser][i] = sqrt((JK_BIN_COUNT-1)*(error_term1[iobser][i] - error_term2[iobser][i]));
                if(LOGGING == 1)
                {
                        printf("error[%d][%d] = %.8f\n",iobser,i,error[iobser][i]);
                }
        }

    }

        /*      Writing time evolution of observable to file    */
        
	sprintf(time_evolution_file,"../Sim_dump_frame/Frame_%d/time_evolution.log",NX);
        if(LOGGING == 1)
        {
                printf("Log2 Time Evolution Output File: %s\n",time_evolution_file);
        }

        Foutput = fopen(time_evolution_file, "a");
        if (Foutput == NULL)
        {
                print_and_exit("***************Could Not Open Log2 Time Evolution Output File: %s************************",time_evolution_file);
        }

        for(int i=0;i<log_bin_cnt;i++)
        {
                if(LOGGING == 1)
                {
                        printf("%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%d\t%d\t%.1f\t%.1f\t%d\t%.4f\n",i,step_log2[i],jk_avg[0][i],error[0][i],jk_avg[1][i],error[1][i],jk_avg[2][i],error[2][i],jk_avg[3][i],error[3][i],NX,NY,EPSILON,KAPPA,STRIP_SIZE,k_bT);
                }
                fprintf(Foutput,"%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%d\t%d\t%.1f\t%.1f\t%d\t%.4f\n",i,step_log2[i],jk_avg[0][i],error[0][i],jk_avg[1][i],error[1][i],jk_avg[2][i],error[2][i],jk_avg[3][i],error[3][i],NX,NY,EPSILON,KAPPA,STRIP_SIZE,k_bT);
        }
        fclose(Foutput);

return 0;
}
