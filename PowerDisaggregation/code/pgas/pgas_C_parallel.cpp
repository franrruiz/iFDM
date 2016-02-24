#include "pgas_C_parallel.h"

#define input_Nt prhs[0]
#define input_Nr prhs[1]
#define input_NPF prhs[2]
#define input_NPG prhs[3]
#define input_M prhs[4]
#define input_T prhs[5]
#define input_Q prhs[6]
#define input_Y prhs[7]
#define input_flagPG prhs[8]
#define input_Z prhs[9]
#define input_sy2 prhs[10]
#define input_ptrans prhs[11]
#define input_Xt prhs[12]
#define input_P prhs[13]

#define output_XPG plhs[0]

/* INPUTS:
 *
 *  0: Nt,          #devices (or users)
 *  1: Nr,          #dimensions
 *  2: NPF,         #particles for SMC
 *  3: NPG,         #particles for PGAS
 *  4: M,           #desired samples
 *  5: T,           #time steps
 *  6: Q,           number of states (excluding 0)
 *  7: Y,           observations (Nr x T)
 *  8: flagPG,      flag to choose SMC/PGAS (SMC ignores Z)
 *  9: Z,           particle (integers) we condition on (Nt x T)
 * 10: sy2,         noise variance
 * 11: ptrans,      transition probabilities (Q+1 x Q+1 x Nt)
 * 12: Xt,          matrix to store the particles. It should be initialized as int16 (maxNt x N x T)
 * 13: P,  			unique power values (Q+1 x Nt)
 *
 */

/* OUTPUTS:
 *
 *  0: XPG,     MCMC samples of trajectories (Nt x M x T)
 *
 */

void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] ) {
    
    /**************************** Read inputs ****************************/
    int Nt = mxGetScalar(input_Nt);
    int Nr = mxGetScalar(input_Nr);
    int NPF = mxGetScalar(input_NPF);
    int NPG = mxGetScalar(input_NPG);
    int M = mxGetScalar(input_M);
    int T = mxGetScalar(input_T);
    int Q = mxGetScalar(input_Q);
    double sy2 = mxGetScalar(input_sy2);
    int flagPG = mxGetScalar(input_flagPG);
    
    double *Y = mxGetPr(input_Y);
    double *ptrans = mxGetPr(input_ptrans);
    double *P = mxGetPr(input_P);
    
    double *Z = mxGetPr(input_Z);
    int16_T *Xt = (int16_T*)mxGetData(input_Xt);

    /******************** Allocate memory for output *********************/
    // Note: XPG(:,m,t) is a vector that describes the trasmitted symbols
    // at t, according to the m-th sample of the sequence of symbols
    mwSize *dims = (mwSize*)calloc(3,sizeof(mwSize));
    dims[0] = Nt;
    dims[1] = M;
    dims[2] = T;
    output_XPG = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
    double *XPG = mxGetPr(output_XPG);
    free(dims);
    
    /************** Allocate memory for auxiliary variables **************/
    double *log_ptrans = (double*)calloc(Nt*(Q+1)*(Q+1),sizeof(double));
    int *xc = (int*)calloc(Nt*T,sizeof(int));
    
    /***************************** Main body *****************************/
    // Set the seed for each thread
    gsl_rng *semillaGlobal = gsl_rng_alloc(gsl_rng_taus2);
    time_t clckG = time(NULL);
    gsl_rng_set(semillaGlobal,clckG);
    gsl_rng **semilla = (gsl_rng**)calloc(omp_get_max_threads(),sizeof(gsl_rng*));
    for(int np=0; np<omp_get_max_threads(); np++) {
        semilla[np] = gsl_rng_alloc(gsl_rng_taus2);
        time_t clck = time(NULL);
        unsigned long int auxSeed = 0;
        if(np>0) {
            auxSeed = gsl_rng_uniform_int(semilla[np-1],100000);
        }
        gsl_rng_set(semilla[np],clck+auxSeed);
    }
    
    // Number of particles
    int N;
    
    // Compute log(ptrans)
    for(int nt=0; nt<Nt; nt++) {
	    double aux;
    	for(int q1=0; q1<Q+1; q1++) {
	    	for(int q2=0; q2<Q+1; q2++) {
	    		aux = my_log(getdouble_3D(ptrans,q1,q2,nt,Q+1,Q+1,Nt));
	    		setdouble_3D(aux,log_ptrans,q1,q2,nt,Q+1,Q+1,Nt);
	    	}
    	}
    }
    
    // If Z (initial particle) is provided, copy Z to xc
    if(flagPG) {
        // Copy Z as it is
        for(int nt=0; nt<Nt; nt++) {
            for(int t=0; t<T; t++) {
                setint_2D(getint_2D(Z,nt,t,Nt,T),xc,nt,t,Nt,T);
            }
        }
    }
    
    // Main loop. For each desired sample m
    for(int m=0; m<M; m++) {
        
        // If first iteration, set number of particles to NPF. Else, use NPG.
        if((m==0) && (!flagPG)) {
            N = NPF;
        } else {
            N = NPG;
        }
        
        // Allocate memory for some important variables
        int *a_ind = (int*)calloc(N*T,sizeof(int));  // Stores the ancestor indices
        double **W = (double**)calloc(T,sizeof(double*));  // Stores the particles weights
        double **logW = (double**)calloc(T,sizeof(double*));  // Stores the particles log-weights
        for(int t=0; t<T; t++) {
            W[t] = (double*)calloc(N,sizeof(double));
            logW[t] = (double*)calloc(N,sizeof(double));
        }
        
        // Loop over time
        for(int t=0; t<T; t++) {
            // At time t=0 we propagate particles from the inactive state
            if(t==0) {
                // Propagate particles
                #pragma omp parallel for
                for(int n=0; n<N; n++) {
                	double *aux_vec = new double[Q+1];
                	int aux_int;
                    for(int nt=0; nt<Nt; nt++) {
                    	for(int q=0; q<Q+1; q++) {
                    		aux_vec[q] = getdouble_3D(ptrans,0,q,nt,Q+1,Q+1,Nt);
                    	}
                    	gsl_ran_discrete_t *aux_prob = gsl_ran_discrete_preproc(Q+1,aux_vec);
                    	aux_int = gsl_ran_discrete(semilla[omp_get_thread_num()],aux_prob);
                    	setint_Xt(aux_int,Xt,nt,n,t,Nt,N,T);
                    	gsl_ran_discrete_free(aux_prob);
                    }
                    delete [] aux_vec;
                }
                // If m>0, the N-th particle should be xc, i.e., Xt(:,N,t)=xc(:,t)
                if((m!=0)||(flagPG)) {
                    for(int nt=0; nt<Nt; nt++) {
                        setint_Xt(getint_2D(xc,nt,t,Nt,T),Xt,nt,N-1,t,Nt,N,T);
                    }
                }
            // For a general time instant t (t!=0)
            } else {
                // First, resampling the particles at t-1 (saves it in a_ind)
                gsl_ran_discrete_t *auxTable = gsl_ran_discrete_preproc(N,W[t-1]);
                #pragma omp parallel for shared(auxTable)
                for(int n=0; n<N; n++) {
                    setint_2D(gsl_ran_discrete(semilla[omp_get_thread_num()],auxTable),a_ind,n,t,N,T);
                }
                gsl_ran_discrete_free(auxTable);
                
                // Second, propagate the particles from t-1 to t
                #pragma omp parallel for
                for(int n=0; n<N; n++) {
                    int prev;
                	double *aux_vec = new double[Q+1];
                    for(int nt=0; nt<Nt; nt++) {
                        // Obtain the state at t-1
                        prev = getint_Xt(Xt,nt,getint_2D(a_ind,n,t,N,T),t-1,Nt,N,T);
                        // Create vector of transition probabilities from previous state
	                	int aux_int;
                    	for(int q=0; q<Q+1; q++) {
                    		aux_vec[q] = getdouble_3D(ptrans,prev,q,nt,Q+1,Q+1,Nt);
                    	}
                    	gsl_ran_discrete_t *aux_prob = gsl_ran_discrete_preproc(Q+1,aux_vec);
                    	// Propagate particle
                    	aux_int = gsl_ran_discrete(semilla[omp_get_thread_num()],aux_prob);
                    	setint_Xt(aux_int,Xt,nt,n,t,Nt,N,T);
                    	gsl_ran_discrete_free(aux_prob);
                    }
                    delete [] aux_vec;                        
                }
                
                // For PGAS, the ancestor index for the N-th particle is sampled differently
                if((m!=0)||(flagPG)) {
                    // Set the N-th particle to xc, i.e., Xt(:,N,t)=xc(:,t)
                    for(int nt=0; nt<Nt; nt++) {
                        setint_Xt(getint_2D(xc,nt,t,Nt,T),Xt,nt,N-1,t,Nt,N,T);
                    }
                    
                    // Factor 1 in Eq. 5: Likelihood. It is equal to 0 because there's no memory (L=1)
                    
                    // Factor 2 in Eq. 5: Transition probabilities from t-1
                    // to the particle that we condition on at time t
                    double *logWZ = computeFactor2(Nt,N,t,T,Q,log_ptrans,xc,Xt);

                    // Factor 3 are just the (log-)weights for time t-1
                    
                    // Add together the three contributions in variable w_a
                    double *logw_a = (double*)calloc(N,sizeof(double));
                    double *w_a = (double*)calloc(N,sizeof(double));
                    double max_w_a = -INFINITY;
                    
                    for(int n=0; n<N; n++) {
                        logw_a[n] = logWZ[n]+logW[t-1][n];
                        // Keep track of the maximum of log(w_a)
                        max_w_a = (max_w_a>logw_a[n])?max_w_a:logw_a[n];
                    }
                    #pragma omp parallel for
                    for(int n=0; n<N; n++) {
                        // Compute the weights (they are NOT normalized because
                        // gsl_ran_discrete doesn't need them to be normalized)
                        w_a[n] = my_exp(logw_a[n]-max_w_a);
                    }
                    
                    // Sample the N-th ancestor and store it in a_ind(N,t)
                    gsl_ran_discrete_t *auxTable = gsl_ran_discrete_preproc(N,w_a);
                    setint_2D(gsl_ran_discrete(semillaGlobal,auxTable),a_ind,N-1,t,N,T);
                    gsl_ran_discrete_free(auxTable);
                    
                    // Free memory
                    free(logWZ);
                    free(w_a);
                    free(logw_a);
                }
            }
            
            // Compute the importance weights for each particle, W(:,t)
            double *sumr = (double*)calloc(N*Nr,sizeof(double));
            double maxLogW = -INFINITY;
            #pragma omp parallel for shared(maxLogW)
            for(int n=0; n<N; n++) {
                // For each particle, initialize some variables
                int nCur = n;
                int r = 0;
                double aux1;
                for(int nr=0; nr<Nr; nr++) {
                    sumr[N*nr+n] = 0;
                }
                // Contribution of particle n to Y(:,t)
                while((r<1)&&(t-r>=0)) {
                    for(int nr=0; nr<Nr; nr++) {
                        for(int nt=0; nt<Nt; nt++) {
                        	int qq = getint_forC_Xt(Xt,nt,nCur,t-r,Nt,N,T,NULL,0);
                            sumr[N*nr+n] += getdouble_2D(P,qq,nt,Q+1,Nt);
                        }
                    }
                    nCur = getint_2D(a_ind,nCur,t-r,N,T);
                    r++;
                }
                // Compute the log-likelihood of particle n
                logW[t][n] = 0;
                for(int nr=0; nr<Nr; nr++) {
                    aux1 = getdouble_2D(Y,nr,t,Nr,T)-sumr[N*nr+n];
                    logW[t][n] -= aux1*aux1;
                }
                logW[t][n] /= sy2;
                // Keep track of the maximum log-weight
                #pragma omp critical
                {   
                    maxLogW = (maxLogW>logW[t][n])?maxLogW:logW[t][n];
                }
            }
            
            // Take the exponential of the weights
            #pragma omp parallel for
            for(int n=0; n<N; n++) {
                W[t][n] = my_exp(logW[t][n]-maxLogW);
            }
            
            // Free sumr
            free(sumr);
            
        } // end of for loop in t
        
        // Select one trajectory at random (according to W(:,T))
        gsl_ran_discrete_t *auxTable = gsl_ran_discrete_preproc(N,W[T-1]);
        int nChosen = gsl_ran_discrete(semillaGlobal,auxTable);
        gsl_ran_discrete_free(auxTable);
        
        // Construct the whole trajectory for the chosen particle, based on Xt and a_ind.
        // Save this trajectory into XPG(:,m,:) and xc(:,:)
        for(int t=T-1; t>=0; t--) {
            for(int nt=0; nt<Nt; nt++) {
                setdouble_3D(getint_forC_Xt(Xt,nt,nChosen,t,Nt,N,T,NULL,0),XPG,nt,m,t,Nt,M,T);
                setint_2D(getint_Xt(Xt,nt,nChosen,t,Nt,N,T),xc,nt,t,Nt,T);
            }
            nChosen = getint_2D(a_ind,nChosen,t,N,T);
        }
        
        // Free a_ind, W and logW
        free(a_ind);
        for(int t=0; t<T; t++) {
            free(W[t]);
            free(logW[t]);
        }
        free(W);
        free(logW);
    }
    
    /**************************** Free memory ****************************/
    gsl_rng_free(semillaGlobal);
    for(int np=0; np<omp_get_max_threads(); np++) {
        gsl_rng_free(semilla[np]);
    }
    free(semilla);
    free(log_ptrans);
    free(xc);
}



/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/************************** Auxiliary functions **************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

double *computeFactor2(int Nt, int N, int t, int T, int Q, double *log_ptrans, int *xc, int16_T *Xt) {
    // Allocate memory for result
    double *logWZ = (double*)calloc(N,sizeof(double));
    
    // For each particle, sum the log-transition probability of each transmitter
    #pragma omp parallel for
    for(int n=0; n<N; n++) {
        int prev;
        int curr;
        for(int nt=0; nt<Nt; nt++) {
            prev = getint_Xt(Xt,nt,n,t-1,Nt,N,T);
            curr = getint_2D(xc,nt,t,Nt,T);
            logWZ[n] += getdouble_3D(log_ptrans,prev,curr,nt,Q+1,Q+1,Nt);
        }
    }
    
    // Return result
    return logWZ;
}


/************************ Misc. Utility functions ************************/

inline double my_exp(double val) {
    return (val<-700.0?gsl_sf_exp(-700.0):gsl_sf_exp(val));
}

inline double my_log(double val) {
    return (val<1e-300?gsl_sf_log(1e-300):gsl_sf_log(val));
}

/*************************** Get/Set functions ***************************/

inline double getdouble_2D(double *x, int n1, int n2, int N1, int N2) {
    return x[(unsigned long long)N1*n2+n1];
}

inline double getdouble_3D(double *x, int n1, int n2, int n3, int N1, int N2, int N3) {
    return x[(unsigned long long)N1*N2*n3+(unsigned long long)N1*n2+n1];
}

inline void setdouble_2D(double val, double *x, int n1, int n2, int N1, int N2) {
    x[(unsigned long long)N1*n2+n1] = val;
}

inline void setdouble_3D(double val, double *x, int n1, int n2, int n3, int N1, int N2, int N3) {
    x[(unsigned long long)N1*N2*n3+(unsigned long long)N1*n2+n1] = val;
}

inline int getint_2D(double *x, int n1, int n2, int N1, int N2) {
    return (int)(x[(unsigned long long)N1*n2+n1]);
}

inline int getint_2D(int *x, int n1, int n2, int N1, int N2) {
    return x[(unsigned long long)N1*n2+n1];
}

inline int getint_forC_2D(int *x, int n1, int n2, int N1, int N2, double *head, int offState) {
    int aux = x[(unsigned long long)N1*n2+n1];
    if(aux>=0) {
        return aux;
    } else if(aux==offState) {
        return 0;
    } else {
        return (int)head[-aux-1];
    }
}

inline int getint_3D(double *x, int n1, int n2, int n3, int N1, int N2, int N3) {
    return (int)(x[(unsigned long long)N1*N2*n3+(unsigned long long)N1*n2+n1]);
}

inline int getint_3D(int *x, int n1, int n2, int n3, int N1, int N2, int N3) {
    return x[(unsigned long long)N1*N2*n3+(unsigned long long)N1*n2+n1];
}

inline int getint_Xt(int16_T *x, int n1, int n2, int n3, int N1, int N2, int N3) {
    return x[(unsigned long long)N1*N2*n3+(unsigned long long)N1*n2+n1];
}

inline int getint_forC_Xt(int16_T *x, int n1, int n2, int n3, int N1, int N2, int N3, double *head, int offState) {
    int16_T aux = x[(unsigned long long)N1*N2*n3+(unsigned long long)N1*n2+n1];
    if(aux>=0) {
        return aux;
    } else if(aux==offState) {
        return 0;
    } else {
        return (int)head[-aux-1];
    }
}

inline void setint_2D(int val, int *x, int n1, int n2, int N1, int N2) {
    x[(unsigned long long)N1*n2+n1] = val;
}

inline void setint_3D(int val, int *x, int n1, int n2, int n3, int N1, int N2, int N3) {
    x[(unsigned long long)N1*N2*n3+(unsigned long long)N1*n2+n1] = val;
}

inline void setint_Xt(int val, int16_T *x, int n1, int n2, int n3, int N1, int N2, int N3) {
    x[(unsigned long long)N1*N2*n3+(unsigned long long)N1*n2+n1] = (int16_T)val;
}

inline bool getbool_2D(bool *x, int n1, int n2, int N1, int N2) {
    return x[(unsigned long long)N1*n2+n1];
}

inline void setbool_2D(int val, bool *x, int n1, int n2, int N1, int N2) {
    x[(unsigned long long)N1*n2+n1] = val;
}

/*************************************************************************/

