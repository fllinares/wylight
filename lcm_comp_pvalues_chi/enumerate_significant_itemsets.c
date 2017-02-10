#ifndef _enumerate_significant_itemsets_c_
#define _enumerate_significant_itemsets_c_

/* LIBRARY INCLUDES */
#include<math.h>

/* CODE DEPENDENCIES */
#include"time_keeping.c"
#include"var_declare.h"
#include"transaction_keeping.c"
#include"lcm_var.c"
#include "chi2.h"

/* CONSTANT DEFINES */
#define READ_BUF_SIZ 524288 //Size of the buffer to read chars from file

/* GLOBAL VARIABLES */
FILE* results_file;
// Number of observations, N; and midpoint of interval [0,N], floor(N/2)
int N, N_over_2;
// Number of observations in positive class
int n;

// Number of non-empty transaction
int Neff;
// Original vector of labels (dimension # non-empty transactions)
char *labels;

// Corrected significance threshold
double delta;

// Array with all values of minimum attainable P-value in [0,N] pre-computed
double *psi;
// And some constants which are useful to precompute
double class_ratio, class_ratio_bin;

// Output files
FILE *significant_itemsets_output_file, *pvalues_output_file;

/* FUNCTION DECLARATIONS */
void psi_init();
void get_N_n(char *);
void read_labels_file(char *, char*);
extern void LCMFREQ_output_itemset(int *);
extern double Chi2_sf(double,double); //Chi2 Survival function, from Dominik's EasyGWASCore
// Profiling variables
long long n_significant_patterns;

/* -------------------------------- INITIALISATION AND TERMINATION FUNCTIONS ----------------------------------------- */

/* Initialise the code
 * Input arguments are self-explanatory
 * */
void enum_sig_itemsets_init(double sig_th, char *labels_file){
	int j; //Loop variable
	char *labels_buffer;

	get_N_n(labels_file);

	// Store core constants
	N_over_2 = (N % 2) ? (N-1)/2 : N/2;//floor(N/2)
	delta = sig_th;

	// Allocate memory for the buffer containing the class labels, giving an error if it fails
	labels_buffer = (char *)malloc(N*sizeof(char));
	if(!labels_buffer){
		fprintf(stderr,"Error in function enum_sig_itemsets_init: couldn't allocate memory for array labels_buffer\n");
		exit(1);
	}

	/* Allocate memory for the vector of class labels, with labels of empty transactions removed */
	Neff = root_trans_list.siz1;
	labels = (char *)malloc(Neff*sizeof(char));
	if(!labels){
		fprintf(stderr,"Error in function enum_sig_itemsets_init: couldn't allocate memory for array labels\n");
		exit(1);
	}

	// Read file containing class labels and store them in array labels, taking care of removing labels
	// associated with empty transactions
	read_labels_file(labels_file,labels_buffer);
	// Ensure class 1 is the minority class
	if(n > (N/2)){
		for(j=0; j<N; j++) labels_buffer[j] = !labels_buffer[j];
		n = N-n;
	}
	for(j=0;j<Neff;j++) labels[j] = labels_buffer[non_empty_trans_idx[j]];
	free(labels_buffer);
	// The array containing the indices of all non-empty transactions is no longer needed
	free(non_empty_trans_idx);

	// Initialise cache for log(x!) and psi(x)
	psi_init();

	// Initialise profiling variables
	n_significant_patterns = 0;
}


/* Precompute minimum attainable P-values $\psi(x)$ for all x in [0,N] and store them in array psi */
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
	psi[0] = 1;
	for(x=1; x<=n; x++) {
		num = x*(1-class_ratio); num *= num;
		den = x*(1-((double)x)/N)*class_ratio_bin;
		psi[x] = Chi2_sf(num/den,1);
	}

	// Now, compute the minimum attainable p-values in the range [N-N_over_2,N]
	for(x=(n+1); x<=N_over_2; x++) {
		num = n*(1-((double)x)/N); num *= num;
		den = x*(1-((double)x)/N)*class_ratio_bin;
		psi[x] = Chi2_sf(num/den,1);
	}

	// Finally, since $\psi(x)$ is symmetric around N_over_2, complete the right half by symmetry
	for(x=(N_over_2+1); x<=N; x++) psi[x] = psi[N-x];
}

/* Free all allocated memory and give some output for debugging purposes */
void enum_sig_itemsets_end(){
	// Print results
	fprintf(results_file,"RESULTS\n");
	fprintf(results_file,"\t Corrected significance threshold: %e\n",delta);
	fprintf(results_file,"\t LCM support: %d\n",LCM_th);
	fprintf(results_file,"\t Number of significant patterns found: %lld\n",n_significant_patterns);

	// Free allocated memory
	free(psi);
	free(labels);

	// Close output files
	fclose(results_file);
	fclose(significant_itemsets_output_file);
	fclose(pvalues_output_file);
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
	int i,j;//Loop iterators
	double Tval, pval, aux_chi, num_precomp, den_precomp; //Variable to hold p-values
	int a; //Cell count of current itemset

	// Sanity-check
	if (x != current_trans.siz) printf("Error: x = %d, current_trans.siz=%d\n",x,current_trans.siz);

	// Minimum attainable P-value for the hypothesis
	double psi_x = psi[x];
	// Check if the newly found solution is in the current testable region Sigma_k
	if(psi_x > delta) return;

	// Precompute common parts of Chi2 statistic
	aux_chi = ((double)x)/N; num_precomp = -n*aux_chi; den_precomp = x*(1-aux_chi)*class_ratio_bin;
	if(den_precomp > 0){//If den_precomp==0, then pval=1 and the itemset cannot be significant
		// Compute the cell-count corresponding to the current itemset
		a = 0;
		for(i=0; i<current_trans.siz; i++) a += labels[current_trans.list[i]];
		// Obtain Chi2 test statistic and its associated p-value
		Tval = a+num_precomp; Tval *= Tval; Tval /= den_precomp;
		pval = Chi2_sf(Tval,1);
		// If p-value is significant, write current itemset and the corresponding p-value to the output files
		if(pval <= delta){
			n_significant_patterns++;
			fprintf(pvalues_output_file,"%d,%d,%.18e\n",a,x,pval);
			fprintf_current_itemset();
		}
	}

}

/* Process a solution involving the most frequent item (item 0) */
// x = frequency (i.e. number of occurrences) of newly found solution
void process_solution0(int x){
	int i,j;//Loop iterators
	double Tval, pval, aux_chi, num_precomp, den_precomp; //Variable to hold p-values
	int a; //Cell count of current itemset

	// Sanity-check
	if (x != bm_trans_list[1].siz) printf("Error: x = %d, bm_trans_list[1].siz=%d\n",x,bm_trans_list[1].siz);

	// Minimum attainable P-value for the hypothesis
	double psi_x = psi[x];
	// Check if the newly found solution is in the current testable region Sigma_k
	if(psi_x > delta) return;

	// Precompute common parts of Chi2 statistic
	aux_chi = ((double)x)/N; num_precomp = -n*aux_chi; den_precomp = x*(1-aux_chi)*class_ratio_bin;
	if(den_precomp>0){//If den_precomp==0, then pval=1 and the itemset cannot be significant
		// Compute the cell-count corresponding to the current itemset
		a = 0;
		for(i=0; i<bm_trans_list[1].siz; i++) a += labels[bm_trans_list[1].list[i]];
		// Obtain Chi2 test statistic and its associated p-value
		Tval = a+num_precomp; Tval *= Tval; Tval /= den_precomp;
		pval = Chi2_sf(Tval,1);
		// If p-value is significant, write current itemset and the corresponding p-value to the output files
		if(pval <= delta){
			n_significant_patterns++;
			fprintf(pvalues_output_file,"%d,%d,%.18e\n",a,x,pval);
			LCMFREQ_output_itemset(LCM_add.q+LCM_add.t);
		}
	}
}

/* Process a solution involving the array-list represented itemsets */
// x = frequency (i.e. number of occurrences) of newly found solution
// L = pointer to TRANS_LIST struct keeping track of merged transactions
// item = current node of the tree
void ary_process_solution(int x, TRANS_LIST *L, int item, int *mask){
	int j;//Loop iterator
	int aux; //Auxiliary counter
	int *t, *t_end, *ptr, *end_ptr; //Pointers for iterating on transaction list
	double Tval, pval, aux_chi, num_precomp, den_precomp; //Variable to hold p-values
	int a; //Cell count of current itemset

	/* First, process the new hypothesis */

	// Minimum attainable P-value for the hypothesis
	double psi_x = psi[x];
	// Check if the newly found solution is in the current testable region Sigma_k
	if(psi_x > delta) return;

	// Precompute common parts of Chi2 statistic
	aux_chi = ((double)x)/N; num_precomp = -n*aux_chi; den_precomp = x*(1-aux_chi)*class_ratio_bin;
	if(den_precomp>0){//If den_precomp==0, then pval=1 and the itemset cannot be significant
		// Compute the cell-count corresponding to the current itemset, plus sanity-check
		a = 0; aux = 0;
		for(t=LCM_Os[item],t_end=LCM_Ot[item];t<t_end;t++){
			end_ptr = (*t == (L->siz2-1)) ? L->list + L->siz1 : L->ptr[*t+1];
			for(ptr = L->ptr[*t];ptr < end_ptr;ptr++){
				a += labels[*ptr];
				aux++;
			}
		}
		if (x != aux) printf("Error: x = %d, trans_size=%d\n",x,aux);
		// Obtain Chi2 test statistic and its associated p-value
		Tval = a+num_precomp; Tval *= Tval; Tval /= den_precomp;
		pval = Chi2_sf(Tval,1);
		// If p-value is significant, write current itemset and the corresponding p-value to the output files
		if(pval <= delta){
			n_significant_patterns++;
			fprintf(pvalues_output_file,"%d,%d,%.18e\n",a,x,pval);
			fprintf_current_itemset();
		}
	}
}

/* Do a first scan of the file containing the class labels to compute the total number of observations, N,
 * and the total number of observations in the positive class, n
 * */
void get_N_n(char *labels_file){
	FILE *f_labels;//Stream with file containing class labels
	int n_read;//Number of chars read
	int i;// Iterator variable to be used in loops
	char char_to_int[256];//Array for converting chars to int fast
	char *read_buf, *read_buf_aux, *read_buf_end;//Buffer for reading from file and extra pointers for loops

	// Initialise both counters to 0 (the variables are defined as global variables in wy.c)
	N = 0; n = 0;

	//Try to open file, giving an error message if it fails
	if(!(f_labels = fopen(labels_file,"r"))){
		fprintf(stderr, "Error in function get_N_n when opening file %s\n",labels_file);
		exit(1);
	}

	//Try to allocate memory for the buffer, giving an error message if it fails
	read_buf = (char *)malloc(READ_BUF_SIZ*sizeof(char));
	if(!read_buf){
		fprintf(stderr,"Error in function read_labels_file: couldn't allocate memory for array read_buf\n");
		exit(1);
	}

	//Initialize the char to int converter
	for(i=0;i<256;i++) char_to_int[i] = 127;
	// We only care about the chars '0' and '1'. Everything else is mapped into the same "bucket"
	char_to_int['0'] = 0; char_to_int['1'] = 1;

	// Read the entire file
	while(1){
		// Try to read READ_BUF_SIZ chars from the file containing the class labels
		n_read = fread(read_buf,sizeof(char),READ_BUF_SIZ,f_labels);
		// If the number of chars read, n_read_ is smaller than READ_BUF_SIZ, either the file ended
		// or there was an error. Check if it was the latter
		if((n_read < READ_BUF_SIZ) && !feof(f_labels)){
			fprintf(stderr,"Error in function read_labels_file while reading the file %s\n",labels_file);
			exit(1);
		}
		// Process the n_read chars read from the file
		for(read_buf_aux=read_buf,read_buf_end=read_buf+n_read;read_buf_aux<read_buf_end;read_buf_aux++){
			//If the character is anything other than '0' or '1' go to process the next char
			if(char_to_int[*read_buf_aux] == 127) continue;
			N++;
			if(char_to_int[*read_buf_aux]) n++;
		}
		// Check if the file ended,. If yes, then exit the while loop
		if(feof(f_labels)) break;
	}

	//Close the file
	fclose(f_labels);

	//Free allocated memory
	free(read_buf);
}

void read_labels_file(char *labels_file, char *labels_buffer){
	FILE *f_labels;//Stream with file containing class labels
	int n_read;//Number of chars read
	int i;// Iterator variable to be used in loops
	char char_to_int[256];//Array for converting chars to int fast
	char *read_buf, *read_buf_aux, *read_buf_end;//Buffer for reading from file and extra pointers for loops
	char *labels_aux = labels_buffer;//Auxiliary pointer to array labels for increments

	//Try to open file, giving an error message if it fails
	if(!(f_labels = fopen(labels_file,"r"))){
		fprintf(stderr, "Error in function read_labels_file when opening file %s\n",labels_file);
		exit(1);
	}

	//Try to allocate memory for the buffer, giving an error message if it fails
	read_buf = (char *)malloc(READ_BUF_SIZ*sizeof(char));
	if(!read_buf){
		fprintf(stderr,"Error in function read_labels_file: couldn't allocate memory for array read_buf\n");
		exit(1);
	}

	//Initialize the char to int converter
	for(i=0;i<256;i++) char_to_int[i] = 127;
	// We only care about the chars '0' and '1'. Everything else is mapped into the same "bucket"
	char_to_int['0'] = 0; char_to_int['1'] = 1;

	// Read the entire file
	while(1){
		// Try to read READ_BUF_SIZ chars from the file containing the class labels
		n_read = fread(read_buf,sizeof(char),READ_BUF_SIZ,f_labels);
		// If the number of chars read, n_read_ is smaller than READ_BUF_SIZ, either the file ended
		// or there was an error. Check if it was the latter
		if((n_read < READ_BUF_SIZ) && !feof(f_labels)){
			fprintf(stderr,"Error in function read_labels_file while reading the file %s\n",labels_file);
			exit(1);
		}
		// Process the n_read chars read from the file
		for(read_buf_aux=read_buf,read_buf_end=read_buf+n_read;read_buf_aux<read_buf_end;read_buf_aux++){
			//If the character is anything other than '0' or '1' go to process the next char
			if(char_to_int[*read_buf_aux] == 127) continue;
			*labels_aux++ = char_to_int[*read_buf_aux];
		}
		// Check if the file ended,. If yes, then exit the while loop
		if(feof(f_labels)) break;
	}

	//Close the file
	fclose(f_labels);

	//Free allocated memory
	free(read_buf);
}

#endif
