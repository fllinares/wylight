#ifndef _fastwy_c_
#define _fastwy_c_

/* LIBRARY INCLUDES */
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

/* CODE DEPENDENCIES */
#include"time_keeping.c"
#include"permutation.c"
#include"fre_pattern.c"
#include"var_declare.h"

/* CONSTANT DEFINES */

/* GLOBAL VARIABLES */
// Number of observations, N; and midpoint of interval [0,N], floor(N/2)
int N, N_over_2;
// Number of observations in positive class
int n;
// Target FWER
double alpha;
// Current permutation
int j_perm;

// Minimum P-value for each permutation
double *min_pval;
// Region thresholds: Sigma_k = [sl1,sl2] U [N-sl2,N-sl1]
int *sl1, *sl2;
// Current P-value threshold (index and actual value)
int k_th;
double delta;
// Flag variable to keep track of the last change done to region Sigma_k
// If flag==1, the last change was sl1++ (shrink on extremes of the W)
// If flag==0, the last change was sl2-- (shrink on center of the W)
int flag;

// Array with all values of log(n!) in [0,N] pre-computed
double *loggamma;
// Logarithm of 1/binom(N,n). This terms appears for every evaluation of the hypergeometric
// PDF, so it makes sense to precompute
double log_inv_binom_N_n;
// Array with all values of minimum attainable P-value in [0,N] pre-computed
double *psi, *psi_sorted;

// Array of flags of length (N+1). flag_frq_examined[i] = 1 is p-values corresponding
// to items with frequency i have already been examined
char *flag_frq_examined;

//P-value cache variables
double **hypergeom_cache;
double **pval_cache;

// Profiling variables
int min_support_across_permutations;
long long n_hypergeom_rows_cache;
long long n_pvalue_rows_cache;
long long n_pvalues_computed;
long long n_cellcounts_computed;
long long effective_total_dataset_frq;

/* FUNCTION DECLARATIONS */
void loggamma_init();
void psi_init();
void decrease_threshold();
int doublecomp(const void*,const void*);

/* -------------------------------- INITIALISATION AND TERMINATION FUNCTIONS ----------------------------------------- */

/* Initialise the Westfall-Young permutation code
 * Input arguments are self-explanatory
 * */
void wy_init(double target_fwer){
	int j; //Loop variable

	// Store core constants
	N_over_2 = (N % 2) ? (N-1)/2 : N/2;//floor(N/2)
	alpha = target_fwer;

	// Initialise cache for log(x!) and psi(x)
	loggamma_init();
	psi_init();

	// Initialise arrays to keep "testable" regions, raising an error if it fails
	sl1 = (int *)malloc((N_over_2+1)*sizeof(double));
	if(!sl1){
		fprintf(stderr,"Error in function wy_init: couldn't allocate memory for array sl1\n");
		exit(1);
	}
	sl2 = (int *)malloc((N_over_2+1)*sizeof(double));
	if(!sl2){
		fprintf(stderr,"Error in function wy_init: couldn't allocate memory for array sl2\n");
		exit(1);
	}
	// Set the (trivial) largest regions
	flag = 1;
	k_th = N_over_2;
	sl1[k_th] = 0; sl2[k_th] = N_over_2;
	delta = psi[sl1[k_th]];
	psi_sorted[k_th] = delta;
	// Compute all other regions, from largest to smallest
	while(k_th > 0) decrease_threshold();

	// Allocate memory for minimum p-values, raising error if it fails
	min_pval = (double *)malloc(J*sizeof(double));
	if(!min_pval){
		fprintf(stderr,"Error in function wy_init: couldn't allocate memory for array min_pval\n");
		exit(1);
	}
	// Initialise all p-values to 1
	for(j=0; j<J; j++) min_pval[j] = 1;

	// Allocate memory for array of flags flag_frq_examined, raising an error if it fails
	flag_frq_examined = (char *)malloc((N+1)*sizeof(char));
	if(!flag_frq_examined){
		fprintf(stderr,"Error in function wy_init: couldn't allocate memory for array flag_frq_examined\n");
		exit(1);
	}

	// Allocate memory for hypergeom cache, raising an error if it fails
	hypergeom_cache = (double **)malloc((N+1)*sizeof(double *));
	if(!hypergeom_cache){
		fprintf(stderr,"Error in function wy_init: couldn't allocate memory for array hypergeom_cache\n");
		exit(1);
	}
	for(j=0; j<=N; j++) hypergeom_cache[j] = ((double *)0);
	n_hypergeom_rows_cache = 0; //Init profiling variable

	// Allocate memory for pval cache, raising an error if it fails
	pval_cache = (double **)malloc((N+1)*sizeof(double *));
	if(!pval_cache){
		fprintf(stderr,"Error in function wy_init: couldn't allocate memory for array hypergeom_cache\n");
		exit(1);
	}
	for(j=0; j<=N; j++) pval_cache[j] = ((double *)0);
	n_pvalue_rows_cache = 0; //Init profiling variable

	n_pvalues_computed = 0; n_cellcounts_computed = 0; effective_total_dataset_frq = 0; //Init profiling variables
	min_support_across_permutations = N; //Init profiling variables
}

/* Precompute values of log(x!) storing them in the array loggamma */
void loggamma_init(){
	int x;
	// Allocate memory for log-gamma cache, raising error if it fails
	loggamma = (double *)malloc((N+1)*sizeof(double));
	if(!loggamma){
		fprintf(stderr,"Error in function loggamma_init: couldn't allocate memory for array loggamma\n");
		exit(1);
	}
	// Initialise cache with appropriate values
	for(x=0;x<=N;x++) loggamma[x] = lgamma(x+1);//Gamma(x) = (x-1)!
	// Initialise log_inv_binom_N_n
	log_inv_binom_N_n = loggamma[n] + loggamma[N-n] - loggamma[N];
}

/* Precompute minimum attainable P-values $\psi(x)$ for all x in [0,N] and store them in array psi */
void psi_init(){
	double xi1;
	int x, x_init;
	// Allocate memory for psi, raising error if it fails
	psi = (double *)malloc((N+1)*sizeof(double));
	if(!psi){
		fprintf(stderr,"Error in function psi_init: couldn't allocate memory for array psi\n");
		exit(1);
	}
	// Allocate memory for psi_sorted, raising error if it fails
	psi_sorted = (double *)malloc((N_over_2+1)*sizeof(double));
	if(!psi_sorted){
		fprintf(stderr,"Error in function psi_init: couldn't allocate memory for array psi_sorted\n");
		exit(1);
	}

	/* Initialise caches with appropriate values */

	// First compute the left side of "the W", i.e. the range [0,n]
	psi[0] = 1;
	//In this range, the recursion $\psi(x)$=$\psi(x-1)$*[(n-(x-1))/(N-(x-1))] can be seen to be correct
	for(x=1; x<=n; x++) psi[x] = (((double)(n-(x-1)))/(N-(x-1)))*psi[x-1];

	// Now, start computing xi1 in the range [N-N_over_2,N] using another recursion, this time
	// starting in N
	// Note that we don't need to store all values, since this will be used only to initialise
	// psi[N_over_2]
	x_init = N-N_over_2;
	xi1 = 1;
	//In this range, the recursion $\xi_{1}(x)$=$\xi_{1}(x+1)$*[((x-1)-n)/(x+1)] can be seen to be correct
	for(x=(N-1); x>=x_init; x--) xi1 = (((double)((x+1)-n))/(x+1))*xi1;

	// Now, use the value of $\xi_{1}(N-N_over_2)$=xi1[0] to get $\psi(N_over_2)$=psi[N_over_2] using the
	// same recursion if N is odd, or simply copy the value of xi1[0] since $\xi_{1}(N-N_over_2)=\psi(N_over_2)$
	// if N is even
	if (N % 2) psi[N_over_2] = (((double)(x_init-n))/x_init)*xi1;
	else psi[N_over_2] = xi1;

	// Now, using psi[N_over_2] compute the right side of "the W", i.e. the range [n+1,N_over_2]
	// using the same recursion as for $\xi_{1}$
	for(x=(N_over_2-1); x > n; x--) psi[x] = (((double)((x+1)-n))/(x+1))*psi[x+1];

	// Finally, since $\psi(x)$ is symmetric around N_over_2, complete the right half by symmetry
	for(x=x_init; x<=N; x++) psi[x] = psi[N-x];

	// Correct minimum attainable P-value in some edge-cases
	if((N % 2)==0){
		if (n == (N/2)) for(x=1; x<N; x++) psi[x] *= 2;
		else psi[N/2] *= 2;
	}
}

/* Decrease the minimum p-value threshold one level
 * Main operations that need to be performed are:
 * 1) Figure out whether we have to shrink "the W" on the left side or the right side, that is, if the current region
 *    is Sigma_{k} = [sl1,sl2] U [N-sl2,N-sl1], we need to figure out if Sigma_{k+1} is of the form
 *    Sigma_{k+1} = [sl1+1,sl2] U [N-sl2,N-sl1-1] (shrink left side) or Sigma_{k+1} = [sl1,sl2-1] U [N-sl2+1,N-sl1-1]
 *    (shrink right side). This is done with help of a binary flag that remembers which of the two types of region
 *    change happened the last time the threshold was decreased.
 * 2) Update variables sl1, sl2 and delta accordingly
 * */
void decrease_threshold(){
	// Flag==1 means the last call to decrease_threshold() shrunk "the W" on the left side
	k_th--;
	if(flag){
		sl1[k_th] = sl1[k_th+1] + 1; sl2[k_th] = sl2[k_th+1];
		// Check what the new case will be
		if (psi[sl1[k_th]] >= psi[sl2[k_th]]) delta = psi[sl1[k_th]];
		else{ delta = psi[sl2[k_th]]; flag = 0; }
	}else{ // Flag==0 means the last call to decrease_threshold() shrunk "the W" on the right side
		sl1[k_th] = sl1[k_th+1]; sl2[k_th] = sl2[k_th+1]-1;
		// Check what the new case will be
		if (psi[sl1[k_th]] >= psi[sl2[k_th]]){ delta = psi[sl1[k_th]]; flag = 1; }
		else delta = psi[sl2[k_th]];
		//No need to update LCM minimum support in this case, since sl1 remains the same
	}
	psi_sorted[k_th] = delta;
}

/* Free all allocated memory and give some output for debugging purposes */
void wy_end(){
	int j, idx_max;
	double delta_corrected;
	// Sort p-values
	qsort(min_pval,J,sizeof(double),doublecomp);
	// Tentative index to corrected significance threshold
	idx_max = (int)floor(alpha*J)-1; delta_corrected = min_pval[idx_max];
	// Check and correct (if necessary) boundary cases
	if(delta_corrected==min_pval[idx_max+1]){
		while(min_pval[--idx_max]==delta_corrected);
		delta_corrected = min_pval[idx_max];
	}
	// Print results
	printf("\nRESULTS\n");
	printf("\t Corrected significance threshold: %e\n",delta_corrected);
	printf("\t FWER at corrected significance threshold: %e\n",floor(idx_max+1)/J);
	printf("\t Minimum support used across all permutations: %d\n",min_support_across_permutations);
	printf("\nMINIMUM P-VALS (%d PERMUTATIONS)\n",J);
	for(j=0;j<(J-1);j++) printf("%e,",min_pval[j]);
	printf("%e\n",min_pval[J-1]);

	// Free allocated memory
	free(loggamma);
	free(psi); free(psi_sorted);
	free(sl1); free(sl2);
	free(min_pval);
	free(flag_frq_examined);
	for(j=0; j<=N; j++) if(hypergeom_cache[j]) free(hypergeom_cache[j]);
	free(hypergeom_cache);
	for(j=0; j<=N; j++) if(pval_cache[j]) free(pval_cache[j]);
	free(pval_cache);
}

/* -------------------------------FUNCTIONS TO COMPUTE FISHER EXACT TEST P-VALUES----------------------------------- */

/* Return PDF of Hypergeometric R.V of parameters x,n,N at a
 * a = Cell count X=1, Y=1
 * x = Frequency of X=1
 * Keeps a cache with all precomputed values to save computing time at the expense of memory
 * */
double hypergeom_pdf(int a, int x){
	int i;
	if(!hypergeom_cache[x]) {
		// Allocate an entire row of the cache matrix
		hypergeom_cache[x] = (double *)malloc((n+1)*sizeof(double));
		if(!hypergeom_cache[x]){
			fprintf(stderr,"Error in function hypergeom_pdf: couldn't allocate memory for array hypergeom_cache[%d]\n",x);
			exit(1);
		}
		for(i=0; i<=n; i++) hypergeom_cache[x][i] = -1;
		n_hypergeom_rows_cache++;
	}
	if(hypergeom_cache[x][a] == -1) hypergeom_cache[x][a] = exp(log_inv_binom_N_n + loggamma[x] + loggamma[N-x] - (loggamma[a] + loggamma[n-a] + loggamma[x-a] + loggamma[(N-n)-(x-a)]));
	return hypergeom_cache[x][a];
}

/* Evaluate Fisher's exact test P-value in a 2x2 contingency table
 * a = Cell count X=1, Y=1
 * x = Frequency of X=1
 * Keeps a cache with all precomputed values to save computing time at the expense of memory
 * */
double fisher_pval(int a, int x){
	int i;
	int a_min, a_max, k;
	double p_left, p_right;
	if(!pval_cache[x]){
		// Allocate an entire row of the cache matrix
		pval_cache[x] = (double *)malloc((n+1)*sizeof(double));
		if(!pval_cache[x]){
			fprintf(stderr,"Error in function fisher_pval: couldn't allocate memory for array pval_cache[%d]\n",x);
			exit(1);
		}
		for(i=0; i<=n; i++) pval_cache[x][i] = -1;
		n_pvalue_rows_cache++;
	}
	if(pval_cache[x][a] == -1){
		a_min = ((n+x-N) > 0) ? (n+x-N) : 0;//max(0,n+x-N)
		a_max = (x < n) ? x : n;//min(x,n)
		for(k=a_min,p_left=0; k <= a; k++) p_left += hypergeom_pdf(k,x);
		for(k=a,p_right=0; k <= a_max; k++) p_right += hypergeom_pdf(k,x);
		pval_cache[x][a] = (p_left < p_right) ? p_left : p_right;//min(p_left,p_right)
		n_pvalues_computed++; //Update profiling variable
	}
	return pval_cache[x][a];
}

/* --------------------------------CORE FAST WESTFALL-YOUNG PERMUTATION FUNCTIONS------------------------------------ */

/* Compute minimum P-value across all itemsets with frequency=x */
double processItemsets(int x){
	int i,j; // Loop iterator variables
	int a; // Itemset cell-count
	double pval; // Itemset p-value
	double pval_min; // Minimum p-value across itemsets of frequency x
	int *tmp_ptr;

	// Check if we have already examined itemsets with that frequency
	// and, if yes, skip unnecessary recalculation
	if(flag_frq_examined[x]) return -1;
	flag_frq_examined[x] = 1; // Mark as examined for future iterations
	// Check if the patterns are not already in memory and, in that case,
	// call LCM and load results
	if(!frq_queried_lcm[x]) load_patterns(x);
	// Now we know patterns are in memory and can proceed with calculation of
	// their respective p-values
	pval_min = 1;
	for(i=0; i<frq_cnt[x]; i++){
		// Reset cell count
		a = 0;
		// Compute cell count
		//tmp_ptr = trans_per_frq[x][i];
		//for(j=0; j<x; j++) a += labels_perm[tmp_ptr[j]];
		for(j=0; j<x; j++) a += labels_perm[trans_per_frq[x][i][j]];
		n_cellcounts_computed++; //Update profiling variable
		// Obtain the corresponding P-value
		pval = fisher_pval(a,x);
		// Update minimum P-value if necessary
		if(pval < pval_min) pval_min = pval;
	}
	effective_total_dataset_frq += ((long long)frq_cnt[x])*x;
	return pval_min;
}

#if 0
void checkTransactions(){
	FILE *out_str_debug2;
	int x, i, j, aux, n_incorrect;
	printf("Checking transactions\n");
	n_incorrect = 0;
	if(!(out_str_debug2 = fopen("out_str_debug2.txt","w"))){
		fprintf(stderr, "Error in function checkTransactions when opening file out_str_debug2\n");
		exit(1);
	}
	for(x=0; x<=N; x++){
		if(!frq_queried_lcm[x]) { fprintf(out_str_debug2,"x=%d not queried yet: frq_cnt[%d]=%d\n",x,x,frq_cnt[x]); continue; }
		fprintf(out_str_debug2,"x=%d already queried: frq_cnt[%d]=%d\n",x,x,frq_cnt[x]);
		for(i=0; i<frq_cnt[x]; i++){
			for(j=0; j<x; j++){
				aux = trans_per_frq[x][i][j];
				if((aux < 0) | (aux>=N)) {fprintf(out_str_debug2,"Incorrect transaction: trans_per_frq[%d][%d][%d]=%d\n",x,i,j,trans_per_frq[x][i][j]); n_incorrect++;}
			}
		}
	}

	fprintf(out_str_debug2,"Transactions checked: n_incorrect=%d\n",n_incorrect);
	fclose(out_str_debug2);
}
#endif

void computeMinPval(){
	int i,j; // Loop iterator variables
	int a,x; // Current cell-count and frequency
	double pval; // Current p-value and minimum p-value

	double time_perm_init, time_perm_end; // Variables to measure time to compute minimum P-value of current permutation
	double time_minpval_init; // Accumulated time computing cell-counts and P-values before the computation of the current minimum P-val starts
	double time_fileio_init; // Accumulated time reading LCM output files before the computation of the current minimum P-val starts
	double time_mining_init; // Accumulated time in executions of LCM before the computation of the current minimum P-val starts

	// Initialise time counter
	time_minpval_init = time_minpval; time_fileio_init = time_fileio; time_mining_init = time_mining;
	time_perm_init = measureTime();

	// Generate a fresh permutation of the labels
	randperm(labels_perm, labels, N, f_perm);
	// Clear flags
	for(i=0; i<=N; i++) flag_frq_examined[i] = 0;
	// Start from the smallest "testable" region, i.e. smallest attainable p-value
	k_th = 0; delta = psi_sorted[k_th];
	// Find the minimum p-value via decremental search
	while(1){
		// Start inspecting patterns in the "left side of the W"
		for(x=sl1[k_th]; x<=sl2[k_th]; x++){
			pval = processItemsets(x);
			if(pval < 0) continue;
			if(pval < min_pval[j_perm]) min_pval[j_perm] = pval;
		}
		// Now continue with "the right side"
		for(x=(N-sl2[k_th]); x<=(N-sl1[k_th]); x++){
			pval = processItemsets(x);
			if(pval < 0) continue;
			if(pval < min_pval[j_perm]) min_pval[j_perm] = pval;
		}
		// Check if the current minimum p-value is already below the bound
		if((min_pval[j_perm] <= delta) || (k_th==N_over_2)) break;
		// Otherwise increase threshold and continue
		k_th++; delta = psi_sorted[k_th];
	}

	// Update minimum support across permutations for profiling purposes
	if(sl1[k_th] < min_support_across_permutations) min_support_across_permutations = sl1[k_th];

	// Compute elapsed time
	time_perm_end = measureTime(); // This time does NOT include time elapsed in the execution of children processes, such as "black-box" calls to LCM
	time_minpval += (time_perm_end-time_perm_init) - (time_fileio-time_fileio_init); //Obtain time elapsed in computation of P-value without taking into account time elapsed reading LCM output files

	// The function has finished at this point
	printf("Permutation=%d \t Minimum P-value=%e \t Threshold_index = %d \t Minimum Support=%d\n",j_perm,min_pval[j_perm],k_th,sl1[k_th]);
	printf("Total time elapsed in permutation %d: %f (s).\n",j_perm,(time_perm_end-time_perm_init)+(time_mining-time_mining_init));
	printf("Mining time in permutation %d: %f (s).\n",j_perm,time_mining-time_mining_init);
	printf("File I/O time in permutation %d: %f (s).\n",j_perm,time_fileio-time_fileio_init);
	printf("Time computing cell-counts and P-values in permutation %d: %f (s).\n",j_perm,time_minpval-time_minpval_init);
}


/* -----------------------------------------AUXILIARY FUNCTIONS------------------------------------------------------ */
// Comparison function used by quicksort implementation in C
int doublecomp(const void* elem1, const void* elem2){
    if(*(const double*)elem1 < *(const double*)elem2)
        return -1;
    return *(const double*)elem1 > *(const double*)elem2;
}

int main(int argc, char *argv[]){
	int i;
	char tmp_random_filename[128];
	int seed_idx;

	if (argc < 8){
		printf("LCM_WY: J alpha sigma_lamp class_labels_file transactions_file tmp_filename seed\n");
		exit(1);
	}

	// Get input arguments
	int n_perm = atoi(argv[1]);
	double target_fwer = atof(argv[2]);
	int sigma_lamp = atoi(argv[3]);
	char *class_labels_filename = argv[4];
	char *trans_filename = argv[5];
	strcpy(tmp_random_filename,argv[6]);
	seed_idx = atoi(argv[7]);

	// Get time tics for both main process and children processes
	time_program_init = measureTime(); time_program_init_ch = measureTimeChildren();

	// Initialise code
	printf("Initialising...\n");
	tic = measureTime();
	permutation_init(n_perm,class_labels_filename,seed_idx);
	wy_init(target_fwer);
	init_frequent_patterns(tmp_random_filename,trans_filename);
	toc = measureTime(); time_initialisation = toc-tic;
	printf("Done!\n");

	// Load initial patterns (sigma_lamp chosen by executing LAMP a priori and inputed as an argument)
	load_patterns_first(sigma_lamp);

	// Compute minimum P-values for each permutation
	for(j_perm=0; j_perm < J; j_perm++) computeMinPval();

	// Get time tocs for both main process and children processes
	time_program_end = measureTime(); time_program_end_ch = measureTimeChildren();

	// Termination
	tic = measureTime();
	permutation_end();
	wy_end();
	end_frequent_patterns();
	toc = measureTime(); time_termination = toc-tic;

	// Output profiling information
	profileCode();

	exit(0);
}

#endif
