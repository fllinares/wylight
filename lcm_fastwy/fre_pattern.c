#ifndef _fre_pattern_c_
#define _fre_pattern_c_


/* LIBRARY INCLUDES */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

/* CODE DEPENDENCIES */
#include"time_keeping.c"
#include"var_declare.h"

/* CONSTANT DEFINES */
#define READ_BUF_SIZ 524288 //Size of the buffer to read chars from file

/* GLOBAL VARIABLES */

// Filenames
char *trans_database_filename;
char *tmp_filename_lcm;
// Frequent patterns
int ***trans_per_frq;
int *frq_cnt;
char *frq_queried_lcm;
// Profiling variables
long long n_itemsets_in_memory;
long long n_transactions_in_memory;


void init_frequent_patterns(char *tmp_filename, char *trans_filename){
	int i;

	// Store filenames
	tmp_filename_lcm = tmp_filename;
	trans_database_filename = trans_filename;

	// Try to allocate the array of pointers which will hold all transaction lists
	trans_per_frq = (int ***)malloc((N+1)*sizeof(int **));
	if(!trans_per_frq){
		fprintf(stderr,"Error in function init_frequent_patterns: couldn't allocate memory for array trans_per_frq\n");
		exit(1);
	}
	for(i=0; i<=N; i++) trans_per_frq[i] = ((int **)0); //And initialise all pointers to NULL
	//Try to allocate memory for the frequency counter, giving an error message if it fails
	frq_cnt = (int *)malloc((N+1)*sizeof(int));
	if(!frq_cnt){
		fprintf(stderr,"Error in function init_frequent_patterns: couldn't allocate memory for array frq_cnt\n");
		exit(1);
	}
	for(i=0; i<=N; i++) frq_cnt[i] = 0; //And set all counters to 0
	//Try to allocate memory for the frequency queriedf flag, giving an error message if it fails
	frq_queried_lcm = (char *)malloc((N+1)*sizeof(char));
	if(!frq_queried_lcm){
		fprintf(stderr,"Error in function init_frequent_patterns: couldn't allocate memory for array frq_queried_lcm\n");
		exit(1);
	}
	for(i=0; i<=N; i++) frq_queried_lcm[i] = 0; //And set all flags to 0

	// Initialise profiling variables
	n_itemsets_in_memory = 0;
	n_transactions_in_memory = 0;
}

void end_frequent_patterns(){
	int i;
	// Free allocated memory
	for(i=0; i<=N; i++){
		if(trans_per_frq[i]){
			free(trans_per_frq[i][0]);
			free(trans_per_frq[i]);
		}
	}
	free(trans_per_frq);
	free(frq_cnt);
	free(frq_queried_lcm);
}

/* First pass through file calculates how much memory must be allocated */
void allocate_space_next_read(char *filename){
	FILE *f_tr;//Stream with file containing class labels
	int n_read;//Number of chars read
	int i,j;// Iterator variables to be used in loops
	int char_to_int[256];//Array for converting chars to int fast
	char c_aux, *read_buf, *read_buf_aux, *read_buf_end;//Buffer for reading from file and extra pointers for loops
	int frq; //Frequency of current itemset (accumulator)

	//Try to open file, giving an error message if it fails
	if(!(f_tr = fopen(filename,"r"))){
		fprintf(stderr, "Error in function allocate_space_next_read when opening file %s\n",filename);
		exit(1);
	}

	//Try to allocate memory for the buffer, giving an error message if it fails
	read_buf = (char *)malloc(READ_BUF_SIZ*sizeof(char));
	if(!read_buf){
		fprintf(stderr,"Error in function allocate_space_next_read: couldn't allocate memory for array read_buf\n");
		exit(1);
	}

	//Initialize the char to int converter
	for(i=0;i<256;i++) char_to_int[i] = 127;
	// We only care about the chars which represent digits. Everything else is mapped into the same "bucket"
	for(c_aux='0'; c_aux<='9'; c_aux++) char_to_int[c_aux] = c_aux-'0';
	char_to_int['\n'] = -1; char_to_int[','] = -2;

	// Read the entire file
	frq = 0;
	while(1){
		// Try to read READ_BUF_SIZ chars from the file containing the class labels
		n_read = fread(read_buf,sizeof(char),READ_BUF_SIZ,f_tr);
		// If the number of chars read, n_read_ is smaller than READ_BUF_SIZ, either the file ended
		// or there was an error. Check if it was the latter
		if((n_read < READ_BUF_SIZ) && !feof(f_tr)){
			fprintf(stderr,"Error in function allocate_space_next_read while reading the file %s\n",filename);
			exit(1);
		}
		// Process the n_read chars read from the file
		for(read_buf_aux=read_buf,read_buf_end=read_buf+n_read;read_buf_aux<read_buf_end;read_buf_aux++){
			// Ignore chars other than digits, commas or newlines
			if(char_to_int[*read_buf_aux] == 127) continue;
			// If the char is a comma or a newline, we have finished processing one transaction
			if((char_to_int[*read_buf_aux]==-2) || (char_to_int[*read_buf_aux]==-1)){
				frq++;
				// If it's a newline, move on to next itemset and reset current frequency counter
				if(char_to_int[*read_buf_aux]==-1){frq_cnt[frq]++; frq = 0;}
			}
		}
		// Check if the file ended,. If yes, then exit the while loop
		if(feof(f_tr)) break;
	}

	//Close file
	fclose(f_tr);

	//Free memory
	free(read_buf);

	// Now allocate space in trans_per_frq for the itemsets which will be read
	for(i=0; i<=N; i++){
		if(frq_cnt[i] && !(frq_queried_lcm[i])){
			trans_per_frq[i] = (int **)malloc(((long long)frq_cnt[i])*sizeof(int *));
			if(!trans_per_frq[i]){
				fprintf(stderr,"Error in function allocate_space_next_read: couldn't allocate memory for array trans_per_frq[%d]\n",i);
				exit(1);
			}
			trans_per_frq[i][0] = (int *)malloc(((long long)frq_cnt[i])*i*sizeof(int));
			if(!trans_per_frq[i][0]){
				fprintf(stderr,"Error in function allocate_space_next_read: couldn't allocate memory for array trans_per_frq[%d][0]\n",i);
				exit(1);
			}
			for(j=1; j<frq_cnt[i]; j++) trans_per_frq[i][j] = trans_per_frq[i][0] + ((long long)i)*j;
		}
	}


}

/* Second pass through file actually loads the transactions into the already allocated memory */
void load_transactions(char *filename){
	FILE *f_tr;//Stream with file containing class labels
	int n_read;//Number of chars read
	int i,j;// Iterator variable to be used in loops
	int char_to_int[256];//Array for converting chars to int fast
	char c_aux, *read_buf, *read_buf_aux, *read_buf_end;//Buffer for reading from file and extra pointers for loops
	int tr_idx, frq, *frq_cnt_tmp, *trans_buffer;

	//Try to open file, giving an error message if it fails
	if(!(f_tr = fopen(filename,"r"))){
		fprintf(stderr, "Error in function load_transactions when opening file %s\n",filename);
		exit(1);
	}

	//Try to allocate memory for the buffer, giving an error message if it fails
	read_buf = (char *)malloc(READ_BUF_SIZ*sizeof(char));
	if(!read_buf){
		fprintf(stderr,"Error in function load_transactions: couldn't allocate memory for array read_buf\n");
		exit(1);
	}

	//Try to allocate memory for the frequency counter, giving an error message if it fails
	frq_cnt_tmp = (int *)malloc((N+1)*sizeof(int));
	if(!frq_cnt_tmp){
		fprintf(stderr,"Error in function load_transactions: couldn't allocate memory for array frq_cnt_tmp\n");
		exit(1);
	}
	for(i=0; i<=N; i++) frq_cnt_tmp[i] = 0;

	// Buffer for current transaction list
	trans_buffer = (int *)malloc(N*sizeof(int));
	if(!trans_buffer){
		fprintf(stderr,"Error in function load_transactions: couldn't allocate memory for array trans_buffer\n");
		exit(1);
	}

	//Initialize the char to int converter
	for(i=0;i<256;i++) char_to_int[i] = 127;
	// We only care about the chars which represent digits. Everything else is mapped into the same "bucket"
	for(c_aux='0'; c_aux<='9'; c_aux++) char_to_int[c_aux] = c_aux-'0';
	char_to_int['\n'] = -1; char_to_int[','] = -2;

	// Initialise other variables
	tr_idx = 0; frq = 0;

	// Read the entire file
	while(1){
		// Try to read READ_BUF_SIZ chars from the file containing the class labels
		n_read = fread(read_buf,sizeof(char),READ_BUF_SIZ,f_tr);
		// If the number of chars read, n_read_ is smaller than READ_BUF_SIZ, either the file ended
		// or there was an error. Check if it was the latter
		if((n_read < READ_BUF_SIZ) && !feof(f_tr)){
			fprintf(stderr,"Error in function load_transactions while reading the file %s\n",filename);
			exit(1);
		}
		// Process the n_read chars read from the file
		for(read_buf_aux=read_buf,read_buf_end=read_buf+n_read;read_buf_aux<read_buf_end;read_buf_aux++){
			if(char_to_int[*read_buf_aux] == 127) continue;
			if((char_to_int[*read_buf_aux]==-2) || (char_to_int[*read_buf_aux]==-1)){
				trans_buffer[frq] = tr_idx; //Append tr_idx to current transaction list
				//printf("tr_idx=%d\n",tr_idx);
				tr_idx = 0; //Reset current transaction
				frq++; //Increase frequency counter
				if(char_to_int[*read_buf_aux]==-1){//Newline
					for(j=0; j<frq; j++) trans_per_frq[frq][frq_cnt_tmp[frq]][j] = trans_buffer[j];
					frq_cnt_tmp[frq]++;
					frq = 0;
				}
				continue;
			}
			tr_idx = 10*tr_idx + char_to_int[*read_buf_aux];
		}
		// Check if the file ended,. If yes, then exit the while loop
		if(feof(f_tr)) break;
	}

	//Close file
	fclose(f_tr);

	//Free memory
	free(read_buf);
	free(frq_cnt_tmp);
	free(trans_buffer);
}

void load_patterns_first(int support){
	int i;
	char command_buf[1024]; //Buffer to store command string to call LCM
	printf("Load all patterns with support >= %d\n",support);
	// Call LCM
	printf("%d, %d\n",support,N);
	sprintf(command_buf, "../lcm_wy_terada_cache/fim_closed %d %d %s %s",support,N,tmp_filename_lcm,trans_database_filename);
	tic = measureTimeChildren();
	system(command_buf);
	toc = measureTimeChildren(); time_mining += toc-tic;
	// Store transactions in memory
	tic = measureTime();
	allocate_space_next_read(tmp_filename_lcm);
	load_transactions(tmp_filename_lcm);
	toc = measureTime(); time_fileio += toc-tic;
	// Delete temporary file
	tic = measureTimeChildren();
	sprintf(command_buf,"rm %s",tmp_filename_lcm);
	system(command_buf);
	toc = measureTimeChildren(); time_fileio += toc-tic;
	// Update frequency queried flags and profile variables
	for(i=support; i<=N; i++) {frq_queried_lcm[i] = 1; n_itemsets_in_memory += frq_cnt[i]; n_transactions_in_memory += ((long long)frq_cnt[i])*i;}
	printf("Number of transaction lists in memory so far: %lld\n",n_itemsets_in_memory);
	printf("Number of transactions in memory so far: %lld\n",n_transactions_in_memory);
	printf("Finished loading patterns!\n");
}

void load_patterns(int support){
	char command_buf[1024]; //Buffer to store command string to call LCM
	printf("Loading patterns with support = %d\n",support);
	// Call LCM
	sprintf(command_buf, "../lcm_wy_terada_cache/fim_closed %d %d %s %s",support,support,tmp_filename_lcm,trans_database_filename);
	tic = measureTimeChildren();
	system(command_buf);
	toc = measureTimeChildren(); time_mining += toc-tic;
	// Store transactions in memory
	tic = measureTime();
	allocate_space_next_read(tmp_filename_lcm);
	load_transactions(tmp_filename_lcm);
	toc = measureTime(); time_fileio += toc-tic;
	// Delete temporary file
	sprintf(command_buf,"rm %s",tmp_filename_lcm);
	tic = measureTimeChildren();
	system(command_buf);
	toc = measureTimeChildren(); time_fileio += toc-tic;
	// Update frequency queried flags and profile variables
	frq_queried_lcm[support] = 1; n_itemsets_in_memory += frq_cnt[support]; n_transactions_in_memory += ((long long)frq_cnt[support])*support;
	printf("Number of transaction lists in memory so far: %lld\n",n_itemsets_in_memory);
	printf("Number of transactions in memory so far: %lld\n",n_transactions_in_memory);
	printf("Finished loading patterns!\n");
}

#endif
