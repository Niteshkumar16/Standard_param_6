#ifndef _NRUNPLC_H_
#define _NRUNPLC_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "clik.h"
#include "param_config_single.h"

extern Config global_config;

// --- Thread-Local Storage for Clik Objects ---
typedef struct {
    clik_object *camspec;
    clik_object *commander; 
    clik_object *lowlike;
    int initialized;
    int n_nuis_total;
} ClikCache;

static ClikCache* get_clik_cache() {
    static ClikCache cache = {0}; 
    return &cache;
}

// --- Initialize Clik Objects ---
int initialize_clik_objects(error **err) {
    ClikCache *cache = get_clik_cache();
    if (cache->initialized) return 1;

    if(global_config.likelihood_count == 0) {
        // Fallback or Error
        return 0;
    }

    // DEBUG PRINT: Show paths being initialized
    for(int i=0; i<global_config.likelihood_count && i<3; i++) {
        printf("Initializing Likelihood %d: %s\n", i+1, global_config.likelihood_paths[i]);
    }

    // Init up to 3 likelihoods
    if(global_config.likelihood_count > 0) cache->camspec = clik_init(global_config.likelihood_paths[0], err);
    if(global_config.likelihood_count > 1) cache->commander = clik_init(global_config.likelihood_paths[1], err);
    if(global_config.likelihood_count > 2) cache->lowlike = clik_init(global_config.likelihood_paths[2], err);

    if (isError(*err)) return 0;

    // Calculate total nuisance parameters expected by the first likelihood (usually CAMspec)
    if(cache->camspec) {
        parname *nuis_names;
        cache->n_nuis_total = clik_get_extra_parameter_names(cache->camspec, &nuis_names, err);
        if(!isError(*err)) free(nuis_names);
    }

    cache->initialized = 1;
    return 1;
}

// --- Helper: Get Dimensions Correctly ---
int get_clik_dimension(clik_object *clikid, error **err) {
    if (!clikid) return 0;
    int has_cl[6], lmax[6];
    clik_get_has_cl(clikid, has_cl, err);
    clik_get_lmax(clikid, lmax, err);
    
    parname *names;
    int n_nuis = clik_get_extra_parameter_names(clikid, &names, err);
    if(n_nuis > 0) free(names);
    
    int dim = n_nuis;
    for(int i=0; i<6; i++) {
        if(has_cl[i] && lmax[i] >= 0) dim += (lmax[i] + 1);
    }
    return dim;
}

// --- Helper: Fill Vector ---
int fill_clik_vector_generic(clik_object *clikid, double *clvec, 
                             const double *TT, const double *TE, 
                             const double *EE, const double *BB,
                             const double *nuisance_vector, int n_nuis_provided, error **err) {
    if (!clikid || !clvec) return 0;
    
    int has_cl[6], lmax[6];
    int offset = 0;
    
    clik_get_has_cl(clikid, has_cl, err);
    clik_get_lmax(clikid, lmax, err);
    
    // 1. Fill Spectra
    for (int i = 0; i < 6; i++) {
        if (has_cl[i] && lmax[i] >= 0) {
            int length = lmax[i] + 1;
            const double *source = NULL;
            switch(i) {
                case 0: source = TT; break; // TT
                case 1: source = EE; break; // EE  
                case 2: source = BB; break; // BB
                case 3: source = TE; break; // TE
            }
            if (source) {
                // Limit copy to 2501 to prevent buffer overflow if lmax > 2500
                int copy_len = (length <= 2501) ? length : 2501;
                memcpy(clvec + offset, source, copy_len * sizeof(double));
            } else {
                memset(clvec + offset, 0, length * sizeof(double));
            }
            offset += length;
        }
    }
    
    // 2. Fill Nuisance
    parname *names;
    int n_expected = clik_get_extra_parameter_names(clikid, &names, err);
    if(n_expected > 0) {
        free(names);
        if (n_nuis_provided >= n_expected) {
            memcpy(clvec + offset, nuisance_vector, n_expected * sizeof(double));
        } else {
            memset(clvec + offset, 0, n_expected * sizeof(double));
        }
    }
    return 1;
}

// --- Main Likelihood Function ---
double run_plc(int rank, double *task, double *cl_tt_in, double *cl_te_in, double *cl_ee_in, double *cl_bb_in) {
    
    error *myerr = initError();
    if (!initialize_clik_objects(&myerr)) return 1e30;
    ClikCache *cache = get_clik_cache();

    // 1. Prepare Spectra
    double TT[2601] = {0}, TE[2601] = {0}, EE[2601] = {0}, BB[2601] = {0};
    for (int l = 2; l <= 2600; l++) {
        double norm = 2.0 * 3.14159265359 / (l * (l + 1.0));
        TT[l] = cl_tt_in[l] * norm;
        TE[l] = cl_te_in[l] * norm; 
        EE[l] = cl_ee_in[l] * norm;
        BB[l] = cl_bb_in[l] * norm;
    }

    // 2. Build Nuisance Vector
    int n_nuis_alloc = (cache->n_nuis_total > 0) ? cache->n_nuis_total : 1;
    double *dynamic_nuisance = (double*)calloc(n_nuis_alloc, sizeof(double));
    
    int task_index = 0;
    for(int i = 0; i < global_config.param_count; i++) {
        ParameterConfig* p = &global_config.params[i];
        double val = p->is_estimated ? task[task_index++] : p->lower_bound;
        
        if (p->usage == USAGE_NUISANCE && p->target_index >= 0 && p->target_index < n_nuis_alloc) {
            dynamic_nuisance[p->target_index] = val;
        }
    }

    // 3. Compute Likelihoods
    double loglike_sum = 0.0;
    clik_object* liks[3] = {cache->camspec, cache->commander, cache->lowlike};
    
    // Array to store individual likelihood results for printing: [0]=CAMspec, [1]=Commander, [2]=LowLike
    double individual_likes[3] = {0.0, 0.0, 0.0}; 

    for(int k=0; k<3; k++) {
        if(liks[k]) {
            int ndim = get_clik_dimension(liks[k], &myerr); 
            if(ndim > 0) {
                double *clvec = (double*)calloc(ndim, sizeof(double));
                if (fill_clik_vector_generic(liks[k], clvec, TT, TE, EE, BB, dynamic_nuisance, n_nuis_alloc, &myerr)) {
                    double res = clik_compute(liks[k], clvec, &myerr);
                    if (!isError(myerr)) {
                        loglike_sum += res;
                        individual_likes[k] = res; // Store for debug print
                    }
                }
                free(clvec);
            }
        }
    }

    free(dynamic_nuisance);

    if (isError(myerr)) return 1e30;

    // Calculate total chi2
    double total = -2.0 * loglike_sum;

    // DEBUG PRINT: Final results for this step
    printf("Rank %d: likelihoods = (%f, %f, %f) total = %f\n", rank, individual_likes[0], individual_likes[1], individual_likes[2], total);

    return total; 
}

void cleanup_clik_objects() {
    ClikCache *cache = get_clik_cache();
    if (cache->camspec) clik_cleanup(&cache->camspec);
    if (cache->commander) clik_cleanup(&cache->commander);
    if (cache->lowlike) clik_cleanup(&cache->lowlike);
    cache->initialized = 0;
}

#endif