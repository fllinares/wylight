#ifndef _time_keeping_c_
#define _time_keeping_c_

/* LIBRARY INCLUDES FOR MEASURING EXECUTION TIME */
#include <time.h> //Already included in original LCM source code above
#include <sys/time.h>
#include <sys/resource.h>
/* END OF LIBRARY INCLUDES FOR MEASURING EXECUTION TIME */

/* CODE DEPENDENCIES */
#include"var_declare.h"

/* GLOBAL VARIABLES (TIME SPENT) */
double tic, toc; //Auxiliary holders of time variables
double time_program_init, time_program_end, time_program_init_ch, time_program_end_ch;
double time_initialisation = 0;
double time_termination = 0;
double time_mining = 0;
double time_fileio = 0;
double time_minpval = 0;

// Measure running time
double measureTime(){
  struct rusage t;
  struct timeval tv,ts;
  getrusage(RUSAGE_SELF, &t);
  tv = t.ru_utime;
  ts = t.ru_stime;
  return tv.tv_sec + ts.tv_sec + ((double)tv.tv_usec + (double)ts.tv_usec) * 1e-6;
}

// Measure running time of children processes (calls to LCM)
double measureTimeChildren(){
  struct rusage t;
  struct timeval tv,ts;
  getrusage(RUSAGE_CHILDREN, &t);
  tv = t.ru_utime;
  ts = t.ru_stime;
  return tv.tv_sec + ts.tv_sec + ((double)tv.tv_usec + (double)ts.tv_usec) * 1e-6;
}

// Measure peak memory usage (not counting the memory used by the children)
size_t measurePeakMemory(){
  struct rusage t;
  getrusage(RUSAGE_SELF, &t);
  return (size_t)t.ru_maxrss;
}

// Measure peak memory usage of all children processes (i.e. report the peak memory usage
// of the child process which had the highest peak memory usage)
size_t measurePeakMemoryChildren(){
  struct rusage t;
  getrusage(RUSAGE_CHILDREN, &t);
  return (size_t)t.ru_maxrss;
}

// Display execution time and memory consumption
void profileCode(){
	size_t peak_memory, peak_memory_ch;

	printf("\nCODE PROFILING\n");

	printf("Total execution time: %f (s).\n",(time_program_end-time_program_init) + (time_program_end_ch-time_program_init_ch));
	printf("\t Initialisation time: %f (s).\n",time_initialisation);
	printf("\t Mining time: %f (s).\n",time_mining);
	printf("\t File I/O time: %f (s).\n",time_fileio);
	printf("\t Time computing cell-counts and P-values: %f (s).\n",time_minpval);
	printf("\t\t Number of cell-counts computed: %lld.\n",n_cellcounts_computed);
	printf("\t\t Total dataset frequency: %lld.\n", effective_total_dataset_frq);
	printf("\t\t Number of P-values computed: %lld.\n",n_pvalues_computed);
	printf("\t Time to terminate algorithm: %f (s).\n",time_termination);
	printf("\n");

	peak_memory = measurePeakMemory(); peak_memory_ch = measurePeakMemoryChildren();
	printf("Peak memory consumption: %lld (KB in Linux, B in Mac OS X).\n",peak_memory);
	printf("\t Number of transaction lists in memory: %lld\n",n_itemsets_in_memory);
	printf("\t Number of transactions in memory: %lld\n",n_transactions_in_memory);
	printf("\t\t Memory consumption of transaction keeping: %f (MB).\n",((N+1)*sizeof(int **)+n_itemsets_in_memory*sizeof(int *)+n_transactions_in_memory*sizeof(int))/((double)1048576));
	printf("\t Number of rows in Hypergeom PDF cache: %lld\n",n_hypergeom_rows_cache);
	printf("\t Number of rows in P-value cache: %lld\n",n_pvalue_rows_cache);
	printf("\t\t Memory consumption of caches: %f (MB).\n",(2*(N+1)*sizeof(double *)+(n_hypergeom_rows_cache+n_pvalue_rows_cache)*(n+1)*sizeof(double))/((double)1048576));
	printf("Peak memory consumption of children processes: %lld (KB in Linux, B in Mac OS X).\n",peak_memory_ch);
	printf("\n");
}

#endif
