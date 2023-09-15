/*
    CS342 - Programming Assignment1, 2023
    Author: Minos Stavrakakis - csd4120

*/

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <pthread.h>
#include <string.h>
#include <time.h>
#include "cvector.h"

#include <time.h>
#include <sys/time.h>
#include <stdint.h>

uint64_t get_posix_clock_time ()
{
    struct timespec ts;

    if (clock_gettime (CLOCK_MONOTONIC, &ts) == 0)
        return (uint64_t) (ts.tv_sec * 1000000 + ts.tv_nsec / 1000);
    else
        return 0;
}

#define TRUE 1
#define FALSE 0

#define DAMPING_FACTOR 0.85
#define BASE_RANK 1.0 - DAMPING_FACTOR
#define NO_ITERATIONS 500

typedef struct node_s {
    int64_t id; // = -1; @ init
    double_t rank; // = 1.0; @ init
    double_t rank_to_give; // = -1.0; @ init
    cvector_vector_type(int64_t) vec_nbor_in; // = NULL; @ init
    uint16_t v_in_sz; // = 0; @ init
    cvector_vector_type(int64_t) vec_nbor_out; // = NULL; @ init
    uint16_t v_out_sz; // = 0; @ init
} node;

typedef enum graph_t {
    UNDIRECTED  =522,
    DIRECTED    =723
} graph_t;

typedef struct thread_args {
	int16_t given_id; // thread id
	int16_t total_threads;
    int64_t BGN;
    int64_t END;
} thread_args;

/* GLOBAL & Thread Viewable Variables */
node *node_arr = NULL;

pthread_barrier_t thread_barrier;


void init_node_arr(node *node_arr, int64_t no_nodes) {
    
    for (int64_t i=0; i<no_nodes; ++i) {
        node_arr[i].id = -1;
        node_arr[i].rank = 1.0;
        node_arr[i].rank_to_give = -1.0;
        node_arr[i].vec_nbor_in = NULL;
        node_arr[i].vec_nbor_out = NULL;
        node_arr[i].v_in_sz = 0;
        node_arr[i].v_out_sz = 0;
    }

}

void free_node_arr(node *node_arr, int64_t no_nodes) {
    
    for (int64_t i=0; i<no_nodes; ++i) {

        cvector_free(node_arr[i].vec_nbor_in);
        cvector_free(node_arr[i].vec_nbor_out);
    }
    free(node_arr);
}

void print_node_arr(node *node_arr, int64_t no_nodes, FILE *fp) {
    int64_t i;
    size_t j;

    for (i=0; i<no_nodes; ++i) {

        fprintf(fp, "Id: %ld  \nIncoming: [ ", node_arr[i].id);                                  
        for (j=0; j < cvector_size(node_arr[i].vec_nbor_in); ++j) {
            if (j!=0) {
                fprintf(fp, ", ");
                if (j%10 == 0) fprintf(fp, "\n\t        ");
            }  
            fprintf(fp, "%ld",node_arr[i].vec_nbor_in[j]);
        }
        fprintf(fp, " ] v_in_sz: %u\n", node_arr[i].v_in_sz);

        fprintf(fp, "Outgoing: [ ");
        for (j=0; j < cvector_size(node_arr[i].vec_nbor_out); ++j) {
            if (j!=0) {
                fprintf(fp, ", ");
                if (j%10 == 0) fprintf(fp, "\n\t        ");
            }   
            fprintf(fp, "%ld",node_arr[i].vec_nbor_out[j]);
        }
        fprintf(fp, " ] v_out_sz: %u\n\n", node_arr[i].v_out_sz);
    }

}

void print_ranks(node *node_arr, int64_t no_nodes, FILE *fp) {
    int64_t i;
    for (i=0; i < no_nodes; ++i) {
        fprintf(fp, "%ld,%lf\n", node_arr[i].id, node_arr[i].rank);
    }
}


// Thread Func
void *
page_rank_thrd(void *myargs) {

    thread_args *t_arg = (thread_args*) myargs;
	// Transfer args to local vars for eou
    int64_t t_BGN = t_arg->BGN;
    int64_t t_END = t_arg->END;

    double_t sum_from_in_nbors;
    uint16_t vi;

    /* Calculate rank_to_give for Iteration #1 */
    for (int64_t j = t_BGN; j <= t_END; ++j) {
            // node_arr[j].rank_old = node_arr[j].rank_new;
            node_arr[j].rank_to_give = node_arr[j].rank / ((double_t) node_arr[j].v_out_sz);
        }

    /* BARRIER */
    // make sure rank_to_give have initialized
    pthread_barrier_wait(&thread_barrier); 

    /* WHILE BEGIN */      // NO_ITERATIONS 500
    for (size_t it=0; it < NO_ITERATIONS; ++it) {
        
        for (int64_t i = t_BGN; i <= t_END; ++i) {
        
            sum_from_in_nbors=0.0;

            for (vi=0; vi < node_arr[i].v_in_sz; ++vi) {

                sum_from_in_nbors += node_arr[node_arr[i].vec_nbor_in[vi]].rank_to_give;
            }

            /* Calculate Rank */
            node_arr[i].rank = BASE_RANK + ( DAMPING_FACTOR * sum_from_in_nbors);
            /* Calculate Rank to disperse to nbors */
            node_arr[i].rank_to_give = node_arr[i].rank / ((double_t) node_arr[i].v_out_sz);
        }

        /* BARRIER */
        // make sure all new PR's finished calculating
	    pthread_barrier_wait(&thread_barrier);

    /* WHILE END */
    }

    return NULL;
}




int main(int argc, char const *argv[])
{

    int64_t no_nodes;
    graph_t graph_type;
    int no_threads=4, main_ret=0;
    uint8_t data_stats = FALSE;
    uint8_t time_stats = FALSE;

    char fname_delimit[]="/\\";
    char fname[100], *raw_fname, *fn_tok; // fname should be mallcd

    // get testfile name
    if (argc>1){
        raw_fname = strdup(argv[1]);
        fn_tok = strtok(raw_fname, fname_delimit);
    
        while (fn_tok != NULL) {
            strcpy(fname,fn_tok);
            fn_tok = strtok(NULL, fname_delimit);
        }
        //printf("filename: %s\n", fname);
    }
    else {
        fprintf(stderr, "Please specify an input file.\n");
        return 1;
    }


    // Calculate flags / # of threads
    for (uint8_t i=2; i < argc; ++i) {
        if (strcmp(argv[i], "-datastats") == 0) {
            data_stats = TRUE;
        }
        else if (strcmp(argv[i], "-timestats") == 0) {
            time_stats = TRUE;
        }
        else if (strcmp(argv[i], "-threads") == 0) {
            if (argc > i+1) { 
                no_threads = atoi(argv[i+1]);
                if (no_threads > 4) no_threads=4;
                // else if (no_threads == 3) no_threads=3;
                else if (no_threads < 1) no_threads=1;
                else
                 ;
            }
            else no_threads=4;
            
        }
        else ;
    }   // For Debug
    // printf("# threds: %d\n", no_threads);


    // Determine node_arr size based on input file 
    if (strcmp(fname, "Email-Enron.txt") == 0) {
        no_nodes = 36692;
        graph_type = DIRECTED; // UNDIRECTED! though treated as DIRECTED due to datafile format
    }
    else if (strcmp(fname, "facebook_combined.txt") == 0) {
        no_nodes = 4039;
        graph_type = UNDIRECTED;
    }
    else if (strcmp(fname, "p2p-Gnutella24.txt") == 0) {
        no_nodes = 26518;
        graph_type = DIRECTED;
    }
    // printf("~ filename: %s\n~ no_nodes: %ld\n", fname, no_nodes);


    //~~ Memory allocation for node_arr ~~
    node_arr = (node*)malloc(no_nodes * sizeof(node));
    if (!node_arr) {
        fprintf(stderr, "ERROR: Malloc Failure.\n");
		return 1;
    }
    // initialize node_arr
    init_node_arr(node_arr, no_nodes);


    // --------------------
    // parse file
    FILE *tfp;
    char *txt_line = NULL;
    size_t line_len = 0;
    char line_buf[300] = { '\0' };
    char delimit[]=" \t\r\n\v\f";

    if (!(tfp = fopen(argv[1], "r+"))) {
		fprintf(stderr, "Failed to read input file: %s\n", fname);
		return 1;
	}

    char *token;
    int64_t node_x, nbor_x;
    while (getline(&txt_line, &line_len, tfp)!= -1) {

        strcpy(line_buf, txt_line);
        if(line_buf[0] == '#'){
            ;
        }
        else {

            if(graph_type == UNDIRECTED){
                // node
                token = strtok(line_buf, delimit);
                node_x = atol(token);

                if (node_arr[node_x].id == -1) {
                    node_arr[node_x].id = node_x;
                }

                // nbor
                token = strtok(NULL, delimit);
                nbor_x = atol(token);
                cvector_push_back(node_arr[node_x].vec_nbor_in, nbor_x);
                ++node_arr[node_x].v_in_sz;
                cvector_push_back(node_arr[node_x].vec_nbor_out, nbor_x);
                ++node_arr[node_x].v_out_sz;

                if (node_arr[nbor_x].id == -1) {
                    node_arr[nbor_x].id = nbor_x;
                }
                cvector_push_back(node_arr[nbor_x].vec_nbor_in, node_x);
                ++node_arr[nbor_x].v_in_sz;
                cvector_push_back(node_arr[nbor_x].vec_nbor_out, node_x);
                ++node_arr[nbor_x].v_out_sz;

            }
            else if (graph_type == DIRECTED) {
                // node
                token = strtok(line_buf, delimit);
                node_x = atol(token);

                if (node_arr[node_x].id == -1) {
                    node_arr[node_x].id = node_x;
                }

                // nbor
                token = strtok(NULL, delimit);
                nbor_x = atol(token);
                cvector_push_back(node_arr[node_x].vec_nbor_out, nbor_x);
                ++node_arr[node_x].v_out_sz;

                if (node_arr[nbor_x].id == -1) {
                    node_arr[nbor_x].id = nbor_x;
                }
                cvector_push_back(node_arr[nbor_x].vec_nbor_in, node_x);
                ++node_arr[nbor_x].v_in_sz;
                
            }
            else {
                fprintf(stderr, "ERROR: Unspecified Graph Type.\n");
                main_ret = -1;
            }
        }
    }
    fclose(tfp);


    if (data_stats) {
        FILE *dtfp;
        
        if (!(dtfp = fopen("data_stats.txt", "w+"))) {
			fprintf(stderr, "Failed to create/write data stats file.\n");
			return 1;
		}
        print_node_arr(node_arr, no_nodes, dtfp);
        fclose(dtfp);
    }
    // print_node_arr(node_arr, no_nodes, stdout);

    /* Time Measurement Vars */

    /* 2nd Way of Measurement */
    // uint64_t prev_time_value, time_value;
    // uint64_t time_diff;

    /* 3rd Way of Measurement */
    struct timespec start, finish;
    double elapsed;
    
    
    /* BEGIN TIMER */

    /* 1st Way of Measurement */    // NOT Wall Time
    // clock_t start_tm = clock();
    
    /* 2nd Way of Measurement */
    // prev_time_value = get_posix_clock_time();

    /* 3rd Way of Measurement */
    clock_gettime(CLOCK_MONOTONIC, &start);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\'
    //             Thread Code

    pthread_t my_threads[no_threads];
	thread_args arg[no_threads];

    
    // Distribute workload among threads
    int64_t tmp_cnt=0;
    int64_t arr_prtn = no_nodes / no_threads;
    int16_t i;

    pthread_barrier_init (&thread_barrier,
	 					NULL, no_threads);

    for (i=0; i<no_threads; ++i) {
    	arg[i].given_id = i;
    	arg[i].total_threads = no_threads;
        arg[i].BGN = tmp_cnt;
        tmp_cnt += arr_prtn;
        arg[i].END = (i == (no_threads-1)) ? ((tmp_cnt-1)+(no_nodes % no_threads)) : (tmp_cnt-1);
    	
        pthread_create (&(my_threads[i]), NULL, 
    	page_rank_thrd, (void*)(&arg[i]));
        
        // For Debug
        // printf("thread_%d: BGN=%ld END=%ld\n", i, arg[i].BGN, arg[i].END);
    }

    for (i=0; i<no_threads; ++i) {
		pthread_join((my_threads[i]), NULL);
	}


    /* TIMER END */ 

    /* 1st Way of Measurement */    // NOT Wall Time
    // clock_t end_tm = clock();
    // printf("Time lapsed: %lf\n", (double) (end_tm - start_tm) / CLOCKS_PER_SEC);
    
    /* 2nd Way of Measurement */
    // time_value = get_posix_clock_time();
    // time_diff = time_value - prev_time_value;
    // printf("Time lapsed: %lfsec\n", (double) time_diff / CLOCKS_PER_SEC);

    /* 3rd Way of Measurement */
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    
    /* Print Time stats */
    if (time_stats) {
        FILE *tmfp;
        
        if (!(tmfp = fopen("time_stats.txt", "a+"))) {
			fprintf(stderr, "Failed to create/write data stats file.\n");
			return 1;
		}
        fprintf(tmfp,"Time lapsed: %lfsec\n", elapsed);
        fclose(tmfp);
    }
    // printf("Time lapsed: %lfsec\n", elapsed);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\'

    // Print Results
    FILE *rfp; 
    if (!(rfp = fopen("pagerank.csv", "w+"))) {
		fprintf(stderr, "Failed to create/write rankings output file.\n");
		return 1;
	}
    print_ranks(node_arr, no_nodes, rfp);
    fclose(rfp);
    // print_ranks(node_arr, no_nodes, stdout);


    // Free Section
    free(raw_fname);
    //free(fname);
    free(node_arr);

    return main_ret;
}