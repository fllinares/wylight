#ifndef _wy_c_
#define _wy_c_

/* LIBRARY INCLUDES */
#include<math.h>
#include"../cephes_double/chdtr.c"

/* CODE DEPENDENCIES */
#include"time_keeping.c"
#include"var_declare.h"
#include"permutation.c"
#include"transaction_keeping.c"
#include"lcm_var.c"

/* CONSTANT DEFINES */


/* GLOBAL VARIABLES */
FILE* results_file, *minpvals_file;
// Number of observations, N; and midpoint of interval [0,N], floor(N/2)
int N, N_over_2;
// Number of observations in positive class
int n;
// Target FWER
double alpha;
// Current FWER
double FWER;

// Maximum Chi2 statistic for every permutation
double *max_Tval;
// Region thresholds: Sigma_k = [sl1,sl2] U [N-sl2,N-sl1]
int sl1, sl2;
// Current critical value for Chi2 test
// NOTE: In this implementation, we work with test statistics instead of p-values to avoid having to
// evaluate the CDF of the Chi2 distribution many many times.
double Tth;
// Flag variable to keep track of the last change done to region Sigma_k
// If flag==1, the last change was sl1++ (shrink on extremes of the W)
// If flag==0, the last change was sl2-- (shrink on center of the W)
int flag;

// Array with all values of maximum attainable Chi2 statistic in [0,N] pre-computed
// NOTE: In this implementation, we work with test statistics instead of p-values to avoid having to
// evaluate the CDF of the Chi2 distribution many many times.
double *psi;
// And some constants which are useful to precompute
double class_ratio, class_ratio_bin;

// Cell-count counter
int *a_cnt;

/* FUNCTION DECLARATIONS */
void psi_init();
int doublecomp(const void*,const void*);
// Two extern functions related to the Chi2 distribution
extern double chdtrc(double,double); //Chi2 Survival function, from Cephes library

// Profiling variables
long long n_pvalues_computed;
long long n_cellcounts_computed;
long long effective_total_dataset_frq;

/* -------------------------------- INITIALISATION AND TERMINATION FUNCTIONS ----------------------------------------- */

/* Initialise the Westfall-Young permutation code
 * Input arguments are self-explanatory
 * */
void wy_init(double target_fwer){
	int j; //Loop variable

	// Store core constants
	N_over_2 = (N % 2) ? (N-1)/2 : N/2;//floor(N/2)
	alpha = target_fwer;

	// Initialise cache for psi(x)
	psi_init();

	// And initialise some others
	sl1 = 1; sl2 = N_over_2; Tth = psi[sl1];
	flag = 1;
	FWER = 0;

	// Allocate memory for maximum Chi2 statistics, raising error if it fails
	max_Tval = (double *)malloc(J*sizeof(double));
	if(!max_Tval){
		fprintf(stderr,"Error in function wy_init: couldn't allocate memory for array max_Tval\n");
		exit(1);
	}
	// Initialise maximum Chi2 statistics to 0
	for(j=0; j<J; j++) max_Tval[j] = 0;

	// Allocate memory for cell counts, raising an error if it fails
	a_cnt = (int *)malloc(J*sizeof(int));
	if(!a_cnt){
		fprintf(stderr,"Error in function wy_init: couldn't allocate memory for array a_cnt\n");
		exit(1);
	}
	for(j=0; j<J; j++) a_cnt[j] = 0;

	n_pvalues_computed = 0; n_cellcounts_computed = 0; effective_total_dataset_frq = 0; //Init profiling variables
}


/* Precompute minimum attainable P-values $\psi(x)$ for all x in [0,N] and store them in array psi */
void psi_init(){
	double num, den;
	int x;

	// Allocate memory for psi, raising error if it fails
	psi = (double *)malloc((N+1)*sizeof(double));
	if(!psi){
		fprintf(stderr,"Error in function psi_init: couldn't allocate memory for array psi\n");
		exit(1);
	}

	/* Initialise caches with appropriate values */

	// Precompute some useful constants
	class_ratio = ((double)n)/N; class_ratio_bin = class_ratio*(1-class_ratio);

	// First compute the left side of "the W", i.e. the range [0,n]
	psi[0] = 0;
	for(x=1; x<=n; x++) {
		num = x*(1-class_ratio); num *= num;
		den = x*(1-((double)x)/N)*class_ratio_bin;
		psi[x] = num/den;
	}

	// Now, compute the minimum attainable p-values in the range [N-N_over_2,N]
	for(x=(n+1); x<=N_over_2; x++) {
		num = n*(1-((double)x)/N); num *= num;
		den = x*(1-((double)x)/N)*class_ratio_bin;
		psi[x] = num/den;
	}

	// Finally, since $\psi(x)$ is symmetric around N_over_2, complete the right half by symmetry
	for(x=(N_over_2+1); x<=N; x++) psi[x] = psi[N-x];
}

/* Free all allocated memory and give some output for debugging purposes */
void wy_end(){
	int j, idx_max;
	double Tth_corrected;

	// Sort p-values
	qsort(max_Tval,J,sizeof(double),doublecomp);
	// Tentative index to corrected significance threshold
	idx_max = (int)floor(alpha*J)-1; Tth_corrected = max_Tval[idx_max];
	// Check and correct (if necessary) boundary cases
	if(Tth_corrected==max_Tval[idx_max+1]){
		while(max_Tval[--idx_max]==Tth_corrected);
		Tth_corrected = max_Tval[idx_max];
	}

	// Print results
	fprintf(results_file,"RESULTS\n");
	fprintf(results_file,"\t Corrected significance threshold: %e\n",chdtrc(1,Tth_corrected));
	fprintf(results_file,"\t FWER at corrected significance threshold: %e\n",floor(idx_max+1)/J);
	fprintf(results_file,"\t Final LCM support: %d\n",LCM_th);
	fprintf(results_file,"\t Final P-value lower bound: %e\n",chdtrc(1,Tth));
	fprintf(results_file,"\t FWER at final P-value lower bound: %e\n",FWER);

	fprintf(minpvals_file,"MAXIMUM CHI2 STATISTICS (%d PERMUTATIONS)\n",J);
	for(j=0;j<(J-1);j++) fprintf(minpvals_file,"%e,",max_Tval[j]);
	fprintf(minpvals_file,"%e\n",max_Tval[J-1]);

	// Report the seed for reproducibility
	fprintf(minpvals_file,"\nRANDOM SEED USED FOR KEY GENERATION\n");
	fprintf(minpvals_file,"\t Seed = %u\n",(unsigned)seed);

	// Free allocated memory
	free(psi);
	free(max_Tval);
	free(a_cnt);

	// Close files
	fclose(results_file); fclose(minpvals_file);
}

/* --------------------------------CORE FAST WESTFALL-YOUNG PERMUTATION FUNCTIONS------------------------------------ */

/* Decrease the minimum p-value threshold one level
 * Main operations that need to be performed are:
 * 1) Figure out whether we have to shrink "the W" on the left side or the right side, that is, if the current region
 *    is Sigma_{k} = [sl1,sl2] U [N-sl2,N-sl1], we need to figure out if Sigma_{k+1} is of the form
 *    Sigma_{k+1} = [sl1+1,sl2] U [N-sl2,N-sl1-1] (shrink left side) or Sigma_{k+1} = [sl1,sl2-1] U [N-sl2+1,N-sl1-1]
 *    (shrink right side). This is done with help of a binary flag that remembers which of the two types of region
 *    change happened the last time the threshold was decreased.
 * 2) Update variables sl1, sl2 and Tth accordingly
 * 3) If sl1 has been modified, then the support of LCM has to be modified
 * 4) Since the tentative corrected significance threshold Tth has changed, the FWER needs to be recomputed
 * */
void decrease_threshold(){
	int j; //Loop iterator
	int false_positives; //Number of false positives (a false positive occurs if max_Tval[j] <= Tth)
	// Flag==1 means the last call to decrease_threshold() shrunk "the W" on the left side
	if(flag){
		sl1++; // Shrink Sigma_k on extremes of the W
		// Check what the new case will be
		if (psi[sl1] <= psi[sl2]) Tth = psi[sl1];
		else{ Tth = psi[sl2]; flag = 0; }
		//Update LCM minimum support
		LCM_th = sl1;
		//printf("\n\n\nTHRESHOLD CHANGE!!! NEW THRESHOLD=%d\n\n\n",LCM_th);
	}else{ // Flag==0 means the last call to decrease_threshold() shrunk "the W" on the right side
		sl2--; // Shrink Sigma_k on center of the W
		// Check what the new case will be
		if (psi[sl1] <= psi[sl2]){ Tth = psi[sl1]; flag = 1; }
		else Tth = psi[sl2];
		//No need to update LCM minimum support in this case, since sl1 remains the same
	}
	// Recompute FWER from scratch
	false_positives = 0;
	for(j=0; j<J; j++) false_positives += (max_Tval[j]>=Tth) ? 1 : 0;
	FWER = ((double)false_positives)/J;
}


/* -------------------FUNCTIONS TO PROCESS A NEWLY FOUND TESTABLE HYPOTHESIS-------------------------------------- */

/* This code contains 3 difference functions to process newly found hypotheses. All of them are virtually identical
 * and the only thing which differs is the way the function receives the list of observations (transactions) for
 * which the hypothesis has X=1.
 * LCM has a complex structure, with frequent itemsets being found at 4 different locations in the source code
 * and under 3 different circumstances. Therefore it was needed to introduce differentiation in the way the transaction
 * list is fed to the "solution processing functions" in order to keep the "transaction keeping" overhead minimal.
 *
 * To reuse this code for problems other than frequent itemset mining, the only thing that needs to be modified
 * is the line which computes the cell counts, for example, the following line in bm_process_solution:
 * 		for(i=0; i<current_trans.siz; i++) a += labels_perm[j][current_trans.list[i]];
 * 	There, current_trans.siz is the number of transactions for which the hypothesis has X=1, i.e. the margin x
 * 	of the 2x2 contingency table (note in this case it is redundant with the input argument x of the function)
 * 	Similarly, current_trans.list[i] with i ranging from 0 to (x-1) has the list of indices of the observations
 * 	for which X=1.
 * 	Simply changing those two parameters accordingly allows code reuse.
 * */

/* Process a solution involving the bitmap represented itemsets */
// x = frequency (i.e. number of occurrences) of newly found solution
void bm_process_solution(int x, int item, int *mask){
	int i, j;//Loop iterators
	double Tval, aux_chi, num_precomp, den_precomp; //Variable to hold p-values
	char *labels_perm_aux; //Auxiliary pointer

	// Sanity-check
	if (x != current_trans.siz) printf("Error: x = %d, current_trans.siz=%d\n",x,current_trans.siz);

	/* First, process the new hypothesis */

	// Maximum attainable Chi2 statistic for the hypothesis
	double psi_x = psi[x];
	// Check if the newly found solution is in the current testable region Sigma_k
	if(psi_x < Tth) return;

	#ifdef PROFILE_MINING
	ticp = measureTime();
	#endif

	// Precompute common parts of Chi2 statistic
	aux_chi = ((double)x)/N; num_precomp = -n*aux_chi; den_precomp = x*(1-aux_chi)*class_ratio_bin;
	n_pvalues_computed++; //Update profiling variable

	if(den_precomp>0){//If den_precomp==0, then pval=1 and the min_pvals will never change
		// Compute cell-counts for all J-permutations
		for(i=0; i<current_trans.siz; i++){
			labels_perm_aux = labels_perm[current_trans.list[i]];//Avoid recomputing labels_perm[current_trans.list[i]] all the time
			for(j=0; j<J; j++) a_cnt[j] += labels_perm_aux[j];
		}
		n_cellcounts_computed += J; //Update profiling variable
		effective_total_dataset_frq += x; // Update profiling variable

		// If not, compute permuted P-values for each permutation
		for(j=0; j<J; j++){
			// Obtain Chi2 test statistic and its associated p-value
			Tval = a_cnt[j]+num_precomp; Tval *= Tval; Tval /= den_precomp;
			// Sanity-check
			if(Tval < 0) printf("Negative Chi2 statistic detected in bm_process_solution!: j=%d, x=%d, a=%d, pval=%e.\n",j,x,a_cnt[j],Tval);
			a_cnt[j] = 0;
			// Check if the newly computed Chi2 statistic is larger than current mmaximum Chi2
			// statistic for the permutation
			if(Tval > max_Tval[j]){
				// Check if the increase in the current maximum Chi2 statistic for the j-th permutation
				// causes an increase in the FWER
				if( (Tval>=Tth) && (max_Tval[j]<Tth)) FWER += ((double)1)/J;
				max_Tval[j] = Tval;
			}
		}

		/* Finally, check if the FWER constraint is still satisfied, if not decrease threshold */
		while(FWER > alpha) {
			//printf("Threshold change BM\n");
			decrease_threshold();
			// Correct possible corruption of LCM data structures due to unexpected change in minimum support
			for(i=0; i<item; i++){
				//printf("Item %d, Frq %d, Current_th %d\n",i,LCM_Ofrq[i],LCM_th);
				if(LCM_Ofrq[i]==(LCM_th-1)){
					//printf("Bucket of item %d corrupted after TH change! Parent %d. Current Th%d.\n",i,item,LCM_th);
					LCM_BM_occurrence_delete(i);
					*mask &= ~BITMASK_1[i];
					//printf("Problem fixed!\n");
				}
			}
		}
	}
	#ifdef PROFILE_MINING
	tocp = measureTime(); time_minpval += tocp-ticp;
	#endif
}

/* Process a solution involving the most frequent item (item 0) */
// x = frequency (i.e. number of occurrences) of newly found solution
void process_solution0(int x){
	int i,j;//Loop iterators
	double Tval, aux_chi, num_precomp, den_precomp; //Variable to hold p-values
	char *labels_perm_aux; //Auxiliary pointer

	// Sanity-check
	if (x != bm_trans_list[1].siz) printf("Error: x = %d, bm_trans_list[1].siz=%d\n",x,bm_trans_list[1].siz);

	/* First, process the new hypothesis */

	// Maximum attainable Chi2 statistic for the hypothesis
	double psi_x = psi[x];
	// Check if the newly found solution is in the current testable region Sigma_k
	if(psi_x < Tth) return;

	#ifdef PROFILE_MINING
	ticp = measureTime();
	#endif

	// Precompute common parts of Chi2 statistic
	aux_chi = ((double)x)/N; num_precomp = -n*aux_chi; den_precomp = x*(1-aux_chi)*class_ratio_bin;
	n_pvalues_computed++; //Update profiling variable

	if(den_precomp>0){//If den_precomp==0, then pval=1 and the min_pvals will never change
		// Compute cell-counts for all permutations
		for(i=0; i<bm_trans_list[1].siz; i++){
			labels_perm_aux = labels_perm[bm_trans_list[1].list[i]]; //Avoid recomputing labels_perm[bm_trans_list[1].list[i]] all the time
			for(j=0; j<J; j++) a_cnt[j] += labels_perm_aux[j];
		}
		n_cellcounts_computed += J; //Update profiling variable
		effective_total_dataset_frq += x; // Update profiling variable

		// If not, compute permuted P-values for each permutation
		for(j=0; j<J; j++){
			// Obtain Chi2 test statistic and its associated p-value

			Tval = a_cnt[j]+num_precomp; Tval *= Tval; Tval /= den_precomp;
			// Sanity-check
			if(Tval < 0) printf("Negative Chi2 statistic detected in process_solution0!: j=%d, x=%d, a=%d, pval=%e.\n",j,x,a_cnt[j],Tval);
			a_cnt[j] = 0;
			// Check if the newly computed Chi2 statistic is larger than current mmaximum Chi2
			// statistic for the permutation
			if(Tval > max_Tval[j]){
				// Check if the increase in the current maximum Chi2 statistic for the j-th permutation
				// causes an increase in the FWER
				if( (Tval>=Tth) && (max_Tval[j]<Tth)) FWER += ((double)1)/J;
				max_Tval[j] = Tval;
			}
		}
		/* Finally, check if the FWER constraint is still satisfied, if not decrease threshold */
		while(FWER > alpha) {
			//printf("threshold change 0\n");
			decrease_threshold();
		}
	}
	#ifdef PROFILE_MINING
	tocp = measureTime(); time_minpval += tocp-ticp;
	#endif
}

/* Process a solution involving the array-list represented itemsets */
// x = frequency (i.e. number of occurrences) of newly found solution
// L = pointer to TRANS_LIST struct keeping track of merged transactions
// item = current node of the tree
void ary_process_solution(int x, TRANS_LIST *L, int item, int *mask){
	int j;//Loop iterator
	int aux; //Auxiliary counter
	int *t, *t_end, *ptr, *end_ptr; //Pointers for iterating on transaction list
	double Tval, aux_chi, num_precomp, den_precomp; //Variable to hold p-values
	char *labels_perm_aux; //Auxiliary pointer

	/* First, process the new hypothesis */

	// Maximum attainable Chi2 statistic for the hypothesis
	double psi_x = psi[x];
	// Check if the newly found solution is in the current testable region Sigma_k
	if(psi_x < Tth) return;

	#ifdef PROFILE_MINING
	ticp = measureTime();
	#endif

	// Sanity-check (this one is more complicated due to the way the transactions are stored)
	aux = 0;
	for(t=LCM_Os[item],t_end=LCM_Ot[item];t<t_end;t++){
		end_ptr = (*t == (L->siz2-1)) ? L->list + L->siz1 : L->ptr[*t+1];
		for(ptr = L->ptr[*t];ptr < end_ptr;ptr++) aux++;
	}
	if (x != aux) printf("Error: x = %d, trans_size=%d\n",x,aux);

	// Precompute common parts of Chi2 statistic
	aux_chi = ((double)x)/N; num_precomp = -n*aux_chi; den_precomp = x*(1-aux_chi)*class_ratio_bin;
	n_pvalues_computed++; //Update profiling variable

	if(den_precomp>0){//If den_precomp==0, then pval=1 and the min_pvals will never change
		// Compute cell-counts for all permutations
		for(t=LCM_Os[item],t_end=LCM_Ot[item];t<t_end;t++){
			end_ptr = (*t == (L->siz2-1)) ? L->list + L->siz1 : L->ptr[*t+1];
			for(ptr = L->ptr[*t];ptr < end_ptr;ptr++){
				labels_perm_aux = labels_perm[*ptr];
				for(j=0; j<J; j++) a_cnt[j] += labels_perm_aux[j];
			}
		}
		n_cellcounts_computed += J; //Update profiling variable
		effective_total_dataset_frq += x; // Update profiling variable

		// If not, compute permuted P-values for each permutation
		for(j=0; j<J; j++){
			// Obtain Chi2 test statistic and its associated p-value
			Tval = a_cnt[j]+num_precomp; Tval *= Tval; Tval /= den_precomp;
			// Sanity-check
			if(Tval < 0) printf("Negative Chi2 statistic detected in ary_process_solution!: j=%d, x=%d, a=%d, pval=%e.\n",j,x,a_cnt[j],Tval);
			a_cnt[j] = 0;
			// Check if the newly computed Chi2 statistic is larger than current mmaximum Chi2
			// statistic for the permutation
			if(Tval > max_Tval[j]){
				// Check if the increase in the current maximum Chi2 statistic for the j-th permutation
				// causes an increase in the FWER
				if( (Tval>=Tth) && (max_Tval[j]<Tth)) FWER += ((double)1)/J;
				max_Tval[j] = Tval;
			}
		}

		/* Finally, check if the FWER constraint is still satisfied, if not decrease threshold */
		while(FWER > alpha) {
			//printf("threshold change ary\n");
			decrease_threshold();
			// Correct possible corruption of LCM data structures due to unexpected change in minimum support
			for(j=0; j<LCM_BM_MAXITEM; j++){
				//printf("Item %d, Frq %d, Current_th %d\n",j,LCM_Ofrq[j],LCM_th);
				if(LCM_Ofrq[j]==(LCM_th-1)){
					//printf("Bucket of item %d corrupted after TH change! Parent %d. Current Th%d.\n",j,item,LCM_th);
					LCM_BM_occurrence_delete(j);
					*mask &= ~BITMASK_1[j];
					//printf("Problem fixed!\n");
				}
			}
		}
	}
	#ifdef PROFILE_MINING
	tocp = measureTime(); time_minpval += tocp-ticp;
	#endif

}

/* AUXILIARY FUNCTIONS */
// Comparison function used by quicksort implementation in C (decreasing sorting order)
int doublecomp(const void* elem1, const void* elem2){
    if(*(const double*)elem1 > *(const double*)elem2) return -1;
    else return 1;
}

#endif
