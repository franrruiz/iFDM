#include "pgas_C.h"

#define input_Nt prhs[0]
#define input_Nr prhs[1]
#define input_NPF prhs[2]
#define input_NPG prhs[3]
#define input_M prhs[4]
#define input_T prhs[5]
#define input_L prhs[6]
#define input_Q prhs[7]
#define input_C prhs[8]
#define input_Y prhs[9]
#define input_flagPG prhs[10]
#define input_Z prhs[11]
#define input_H prhs[12]
#define input_sy2 prhs[13]
#define input_am prhs[14]
#define input_bm prhs[15]
#define input_head prhs[16]
#define input_lHead prhs[17]
#define input_Xt prhs[18]
#define input_onOffModel prhs[19]

#define output_XPG plhs[0]

/* INPUTS:
 *
 *  0: Nt,          #transmitters
 *  1: Nr,          #receivers
 *  2: NPF,         #particles for SMC
 *  3: NPG,         #particles for PGAS
 *  4: M,           #desired samples
 *  5: T,           #time steps
 *  6: L,           memory
 *  7: Q,           constellation order, i.e., |A| (including 0)
 *  8: C,           constellation (1 x Q). The first element should be 0
 *  9: Y,           observations (Nt x T)
 * 10: flagPG,      flag to choose SMC/PGAS (SMC ignores Z)
 * 11: Z,           particle (integers) we condition on (Nt x T)
 * 12: H,           channel (Nr x Nt x L)
 * 13: sy2,         noise variance
 * 14: am,          transition probability from 0 to 0 (1 x Nt)
 * 15: bm,          transition probability from active to 0 (1 x Nt)
 * 16: head,        header symbols (integers) (1 x lHead)
 * 17: lHead,       header length
 * 18: Xt,          matrix to store the particles. It should be initialized as int16 (maxNt x N x T)
 * 19: onOffModel,  if true (!=0), a transmitter that becomes inactive cannot start transmitting again
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
    int L = mxGetScalar(input_L);
    int Q = mxGetScalar(input_Q);
    double sy2 = mxGetScalar(input_sy2);
    int flagPG = mxGetScalar(input_flagPG);
    int lHead = mxGetScalar(input_lHead);
    int onOffModel = mxGetScalar(input_onOffModel);
    int offState = -lHead-1;
    
    bool flagComplex = true;
    mxComplexity myComplex = mxCOMPLEX;
    double *Cr = mxGetPr(input_C);
    double *Ci = mxGetPi(input_C);
    double *Yr = mxGetPr(input_Y);
    double *Yi = mxGetPi(input_Y);
    double *Hr = mxGetPr(input_H);
    double *Hi = mxGetPi(input_H);
    if(NULL==Ci) {
        flagComplex = false;
        myComplex = mxREAL;
    }
    
    double *Z = mxGetPr(input_Z);
    double *am = mxGetPr(input_am);
    double *bm = mxGetPr(input_bm);
    double *head = mxGetPr(input_head);
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
    double *logam = (double*)calloc(Nt,sizeof(double));
    double *log1_am = (double*)calloc(Nt,sizeof(double));
    double *logbm = (double*)calloc(Nt,sizeof(double));
    double *log1_bm = (double*)calloc(Nt,sizeof(double));
    int *xc = (int*)calloc(Nt*T,sizeof(int));
    
    /***************************** Main body *****************************/
    // Set the seed
    srand48(time(NULL));
    gsl_rng *semilla = gsl_rng_alloc(gsl_rng_taus);
    time_t clck = time(NULL);
    gsl_rng_set(semilla,clck);
    
    // Number of particles
    int N;
    
    // Compute log(am), log(bm), log(1-am), log(1-bm) and log(Q-1)
    double logQ_1 = my_log(Q-1.0);
    for(int nt=0; nt<Nt; nt++) {
        logam[nt] = my_log(am[nt]);
        logbm[nt] = my_log(bm[nt]);
        log1_am[nt] = my_log(1.0-am[nt]);
        log1_bm[nt] = my_log(1.0-bm[nt]);
    }
    
    // If Z (initial particle) is provided, copy Z to xc
    if(flagPG && !onOffModel) {
        // If there is not header, copy Z as it is
        if(lHead==0) {
            for(int nt=0; nt<Nt; nt++) {
                for(int t=0; t<T; t++) {
                    setint_2D(getint_2D(Z,nt,t,Nt,T),xc,nt,t,Nt,T);
                }
            }
        // Else, copy but use negative integers to denote header symbols
        // (first header symbol is encoded as -1, second as -2, etc.)
        } else {
            for(int nt=0; nt<Nt; nt++) {
                int t = 0;
                int countH = 0;
                if(getint_2D(Z,nt,t,Nt,T)==0) {
                    setint_2D(0,xc,nt,t,Nt,T);
                } else {
                    countH++;
                    setint_2D(-countH,xc,nt,t,Nt,T);
                }
                for(t=1; t<T; t++) {
                    // If it was inactive at t-1
                    if(getint_2D(Z,nt,t-1,Nt,T)==0) {
                        // and it is still inactive
                        if(getint_2D(Z,nt,t,Nt,T)==0) {
                            setint_2D(0,xc,nt,t,Nt,T);
                        // and it is active now
                        } else {
                            countH++;
                            setint_2D(-countH,xc,nt,t,Nt,T);
                        }
                    // If it was active at t-1
                    } else {
                        // and it is inactive now
                        if(getint_2D(Z,nt,t,Nt,T)==0) {
                            countH = 0;
                            setint_2D(0,xc,nt,t,Nt,T);
                        // and it is still active
                        } else {
                            countH++;
                            if(countH>lHead) {
                                setint_2D(getint_2D(Z,nt,t,Nt,T),xc,nt,t,Nt,T);
                            } else {
                                setint_2D(-countH,xc,nt,t,Nt,T);
                            }
                        }
                    }
                }
            }
        }
    // If Z (initial particle) is provided and the on/off model is considered,
    // copy Z but includes symbol 'off' (encoded as offState=-lHead-1)
    } else if(flagPG && onOffModel) {
        int aux;
        // If there is not header, copy Z as it is
        if(lHead==0) {
            for(int nt=0; nt<Nt; nt++) {
                bool flagAct = false;
                bool flagDeact = false;
                for(int t=0; t<T; t++) {
                    aux = getint_2D(Z,nt,t,Nt,T);
                    if(aux>0) {
                        flagAct = true;
                    } else if(aux==0 && flagAct) {
                        flagDeact = true;
                    }
                    setint_2D(flagDeact?offState:aux,xc,nt,t,Nt,T);
                }
            }
        // Else, copy but use negative integers to denote header symbols
        // (first header symbol is encoded as -1, second as -2, etc.)
        } else {
            for(int nt=0; nt<Nt; nt++) {
                bool flagDeact = false;
                int t = 0;
                int countH = 0;
                if(getint_2D(Z,nt,t,Nt,T)==0) {
                    setint_2D(0,xc,nt,t,Nt,T);
                } else {
                    countH++;
                    setint_2D(-countH,xc,nt,t,Nt,T);
                }
                for(t=1; t<T; t++) {
                    // If it was inactive at t-1
                    if(getint_2D(Z,nt,t-1,Nt,T)==0) {
                        // and it is still inactive
                        if(getint_2D(Z,nt,t,Nt,T)==0) {
                            setint_2D(flagDeact?offState:0,xc,nt,t,Nt,T);
                        // and it is active now
                        } else {
                            countH++;
                            setint_2D(-countH,xc,nt,t,Nt,T);
                        }
                    // If it was active at t-1
                    } else {
                        // and it is inactive now
                        if(getint_2D(Z,nt,t,Nt,T)==0) {
                            countH = 0;
                            setint_2D(offState,xc,nt,t,Nt,T);
                            flagDeact = true;
                        // and it is still active
                        } else {
                            countH++;
                            if(countH>lHead) {
                                setint_2D(getint_2D(Z,nt,t,Nt,T),xc,nt,t,Nt,T);
                            } else {
                                setint_2D(-countH,xc,nt,t,Nt,T);
                            }
                        }
                    }
                }
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
                for(int n=0; n<N; n++) {
                    for(int nt=0; nt<Nt; nt++) {
                        // Keep inactive
                        if(drand48()<am[nt]) {
                            setint_Xt(0,Xt,nt,n,t,Nt,N,T);
                        // Activate
                        } else {
                            if(lHead==0) {
                                setint_Xt(1+(int)((Q-1.0)*drand48()),Xt,nt,n,t,Nt,N,T);
                            } else {
                                setint_Xt(-1,Xt,nt,n,t,Nt,N,T);
                            }
                        }
                    }
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
                for(int n=0; n<N; n++) {
                    setint_2D(gsl_ran_discrete(semilla,auxTable),a_ind,n,t,N,T);
                }
                gsl_ran_discrete_free(auxTable);
                
                // Second, propagate the particles from t-1 to t
                int prev;
                for(int n=0; n<N; n++) {
                    for(int nt=0; nt<Nt; nt++) {
                        // Obtain the state at t-1
                        prev = getint_Xt(Xt,nt,getint_2D(a_ind,n,t,N,T),t-1,Nt,N,T);
                        // If it was inactive
                        if(prev==0) {
                            // Keep inactive
                            if(drand48()<am[nt]) {
                                setint_Xt(0,Xt,nt,n,t,Nt,N,T);
                            // Activate
                            } else {
                                if(lHead==0) {
                                    setint_Xt(1+(int)((Q-1.0)*drand48()),Xt,nt,n,t,Nt,N,T);
                                } else {
                                    setint_Xt(-1,Xt,nt,n,t,Nt,N,T);
                                }
                            }
                        // If it was active
                        } else if(prev>offState) {
                            // Deactivate
                            if(drand48()<bm[nt]) {
                                setint_Xt(onOffModel?offState:0,Xt,nt,n,t,Nt,N,T);
                            // Keep active
                            } else {
                                if(prev>0 || -prev>=lHead) {
                                    setint_Xt(1+(int)((Q-1.0)*drand48()),Xt,nt,n,t,Nt,N,T);
                                } else {
                                    setint_Xt(prev-1,Xt,nt,n,t,Nt,N,T);
                                }
                            }
                        }
                        // If it was in the off state
                        else {
                            setint_Xt(offState,Xt,nt,n,t,Nt,N,T);
                        }
                    }
                }
                
                // For PGAS, the ancestor index for the N-th particle is
                // sampled differently
                if((m!=0)||(flagPG)) {
                    
                    // Set the N-th particle to xc, i.e., Xt(:,N,t)=xc(:,t)
                    for(int nt=0; nt<Nt; nt++) {
                        setint_Xt(getint_2D(xc,nt,t,Nt,T),Xt,nt,N-1,t,Nt,N,T);
                    }
                    
                    // Factor 1 in Eq. 5: Likelihood
                    double *logWY;
                    if(L>1) {
                        logWY = computeFactor1(Nt,Nr,N,t,T,L,Q,Cr,Ci,Yr,Yi,a_ind,xc,Xt,Hr,Hi,flagComplex,myComplex,sy2,head,offState);
                    } else {
                        logWY = (double*)calloc(N,sizeof(double));
                    }
                    
                    // Factor 2 in Eq. 5: Transition probabilities from t-1
                    // to the particle that we condition on at time t
                    double *logWZ = computeFactor2(Nt,N,t,T,Q,logQ_1,xc,Xt,logam,logbm,log1_am,log1_bm,lHead,onOffModel,offState);
                    
                    // Factor 3 are just the (log-)weights for time t-1
                    
                    // Add together the three contributions in variable w_a
                    double *logw_a = (double*)calloc(N,sizeof(double));
                    double *w_a = (double*)calloc(N,sizeof(double));
                    double max_w_a = -INFINITY;
                    for(int n=0; n<N; n++) {
                        logw_a[n] = logWY[n]+logWZ[n]+logW[t-1][n];
                        // Keep track of the maximum of log(w_a)
                        max_w_a = ((max_w_a>logw_a[n])?max_w_a:logw_a[n]);
                    }
                    for(int n=0; n<N; n++) {
                        // Compute the weights (they are NOT normalized because
                        // gsl_ran_discrete doesn't need them to be normalized)
                        w_a[n] = my_exp(logw_a[n]-max_w_a);
                    }
                    
                    // Sample the N-th ancestor and store it in a_ind(N,t)
                    gsl_ran_discrete_t *auxTable = gsl_ran_discrete_preproc(N,w_a);
                    setint_2D(gsl_ran_discrete(semilla,auxTable),a_ind,N-1,t,N,T);
                    gsl_ran_discrete_free(auxTable);
                    
                    // Free logWY and logWZ
                    free(logWY);
                    free(logWZ);
                    free(w_a);
                    free(logw_a);
                }
            }
            
            // Compute the importance weights for each particle, W(:,t)
            int nCur;
            int r;
            double *sumr = (double*)calloc(Nr,sizeof(double));
            double *sumi = (double*)calloc(Nr,sizeof(double));
            double aux1, aux2;
            double maxLogW = -INFINITY;
            if(flagComplex) {
                for(int n=0; n<N; n++) {
                    // For each particle, initialize some variables
                    nCur = n;
                    r=0;
                    for(int nr=0; nr<Nr; nr++) {
                        sumr[nr] = 0;
                        sumi[nr] = 0;
                    }
                    // Contribution of particle n to Y(:,t)
                    while((r<L)&&(t-r>=0)) {
                        for(int nr=0; nr<Nr; nr++) {
                            for(int nt=0; nt<Nt; nt++) {
                                sumr[nr] += getdouble_3D(Hr,nr,nt,r,Nr,Nt,L)*Cr[getint_forC_Xt(Xt,nt,nCur,t-r,Nt,N,T,head,offState)] \
                                            -getdouble_3D(Hi,nr,nt,r,Nr,Nt,L)*Ci[getint_forC_Xt(Xt,nt,nCur,t-r,Nt,N,T,head,offState)];
                                sumi[nr] += getdouble_3D(Hr,nr,nt,r,Nr,Nt,L)*Ci[getint_forC_Xt(Xt,nt,nCur,t-r,Nt,N,T,head,offState)] \
                                            +getdouble_3D(Hi,nr,nt,r,Nr,Nt,L)*Cr[getint_forC_Xt(Xt,nt,nCur,t-r,Nt,N,T,head,offState)];
                            }
                        }
                        nCur = getint_2D(a_ind,nCur,t-r,N,T);
                        r++;
                    }
                    // Compute the log-likelihood of particle n
                    logW[t][n] = 0;
                    for(int nr=0; nr<Nr; nr++) {
                        aux1 = getdouble_2D(Yr,nr,t,Nr,T)-sumr[nr];
                        aux2 = getdouble_2D(Yi,nr,t,Nr,T)-sumi[nr];
                        logW[t][n] -= aux1*aux1+aux2*aux2;
                    }
                    logW[t][n] /= sy2;
                    // Keep track of the maximum log-weight
                    maxLogW = ((maxLogW>logW[t][n])?maxLogW:logW[t][n]);
                }
            // For non-complex constellations
            } else {
                for(int n=0; n<N; n++) {
                    // For each particle, initialize some variables
                    nCur = n;
                    r=0;
                    for(int nr=0; nr<Nr; nr++) {
                        sumr[nr] = 0;
                    }
                    // Contribution of particle n to Y(:,t)
                    while((r<L)&&(t-r>=0)) {
                        for(int nr=0; nr<Nr; nr++) {
                            for(int nt=0; nt<Nt; nt++) {
                                sumr[nr] += getdouble_3D(Hr,nr,nt,r,Nr,Nt,L)*Cr[getint_forC_Xt(Xt,nt,nCur,t-r,Nt,N,T,head,offState)];
                            }
                        }
                        nCur = getint_2D(a_ind,nCur,t-r,N,T);
                        r++;
                    }
                    // Compute the log-likelihood of particle n
                    logW[t][n] = 0;
                    for(int nr=0; nr<Nr; nr++) {
                        aux1 = getdouble_2D(Yr,nr,t,Nr,T)-sumr[nr];
                        logW[t][n] -= aux1*aux1;
                    }
                    logW[t][n] /= sy2;
                    // Keep track of the maximum log-weight
                    maxLogW = ((maxLogW>logW[t][n])?maxLogW:logW[t][n]);
                }
            }
            
            // De-normalize the weights
            for(int n=0; n<N; n++) {
                W[t][n] = my_exp(logW[t][n]-maxLogW);
            }
            
            // Free sumr, sumi
            free(sumr);
            free(sumi);
            
        } // end of for loop in t
        
        // Select one trajectory at random (according to W(:,T))
        gsl_ran_discrete_t *auxTable = gsl_ran_discrete_preproc(N,W[T-1]);
        int nChosen = gsl_ran_discrete(semilla,auxTable);
        gsl_ran_discrete_free(auxTable);
        
        // Construct the whole trajectory for the chosen particle, based on Xt and a_ind.
        // Save this trajectory into XPG(:,m,:) and xc(:,:)
        for(int t=T-1; t>=0; t--) {
            for(int nt=0; nt<Nt; nt++) {
                setdouble_3D(getint_forC_Xt(Xt,nt,nChosen,t,Nt,N,T,head,offState),XPG,nt,m,t,Nt,M,T);
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
    gsl_rng_free(semilla);
    free(logam);
    free(log1_am);
    free(logbm);
    free(log1_bm);
    free(xc);
}



/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/************************** Auxiliary functions **************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

double *computeFactor1(int Nt, int Nr, int N, int t, int T, int L, int Q, double *Cr, double *Ci, double *Yr, double *Yi, int *a_ind, int *xc, int16_T *Xt, double *Hr, double *Hi, bool flagComplex, mxComplexity myComplex, double sy2, double *head, int offState) {
    double *logWY = (double*)calloc(N,sizeof(double));
    
    // Allocate memory for auxiliary variables
    mxArray *Ypred_mx = mxCreateDoubleMatrix(Nr,1,myComplex);
    double *Ypredr = mxGetPr(Ypred_mx);
    double *Ypredi = mxGetPi(Ypred_mx);
    
    mxArray **ancContrib_mx = (mxArray**)calloc(L-1,sizeof(mxArray*));
    double **ancContribr = (double**)calloc(L-1,sizeof(double*));
    double **ancContribi = (double**)calloc(L-1,sizeof(double*));
    for(int ll=0; ll<L-1; ll++) {
        ancContrib_mx[ll] = mxCreateDoubleMatrix(Nr,N,myComplex);
        ancContribr[ll] = mxGetPr(ancContrib_mx[ll]);
        ancContribi[ll] = mxGetPi(ancContrib_mx[ll]);
    }
    
    bool *Act_ancB = (bool*)calloc(N*(L-1),sizeof(bool));
    int *Act_anc = (int*)calloc(N*(L-1),sizeof(int));

    // Create auxiliary variable Act_ancB
    for(int n=0; n<N; n++) {
        setbool_2D(true,Act_ancB,n,L-2,N,L-1);   // Act_ancB(1:N,L-1)=true
    }
    int rr = L-2;
    while((rr>=1)&&(t+rr-L+1>=0)) {
        for(int n=0; n<N; n++) {
            if(getbool_2D(Act_ancB,n,rr,N,L-1)) {
                setbool_2D(true,Act_ancB,getint_2D(a_ind,n,t+rr-L+1,N,T),rr-1,N,L-1);
            }
        }
        rr--;
    }
    
    // Create auxiliary variable Act_anc
    int count;
    for(int rr=0; rr<L-1; rr++) {
        count = 0;
        for(int n=0; n<N; n++) {
            if(getbool_2D(Act_ancB,n,rr,N,L-1)) {
                setint_2D(n+1,Act_anc,count,rr,N,L-1);
                count++;
            }
        }
    }
    
    // We increase t. Note: This should not be necessary (of course), but it
    // comes from the fact that the code below was initially designed to work
    // under a Matlab call. Since the indexes in Matlab and the indexes in C
    // differ by 1, we need to artificially increase t now in order for the
    // code below to be correct.
    t++;

    // LOOP: for tau=t:min(t+L-2,T)
    int maxtau = (t+L-2>T?T:t+L-2);
    if(flagComplex) {
        for(int tau=t; tau<=maxtau; tau++) {
            //Ypred = Y(:,tau);
            for(int nr=0; nr<Nr; nr++) {
                Ypredr[nr] = getdouble_2D(Yr,nr,tau-1,Nr,T);
                Ypredi[nr] = getdouble_2D(Yi,nr,tau-1,Nr,T);
            }
            
            // Remove the contribution of the fixed particle to instant tau
            for(int r=0; r<=tau-t; r++) {
                for(int nr=0; nr<Nr; nr++) {
                    for(int nt=0; nt<Nt; nt++) {
                        Ypredr[nr] -= getdouble_3D(Hr,nr,nt,r,Nr,Nt,L)*Cr[getint_forC_2D(xc,nt,tau-r-1,Nt,T,head,offState)] \
                                      -getdouble_3D(Hi,nr,nt,r,Nr,Nt,L)*Ci[getint_forC_2D(xc,nt,tau-r-1,Nt,T,head,offState)];
                        Ypredi[nr] -= getdouble_3D(Hr,nr,nt,r,Nr,Nt,L)*Ci[getint_forC_2D(xc,nt,tau-r-1,Nt,T,head,offState)] \
                                      +getdouble_3D(Hi,nr,nt,r,Nr,Nt,L)*Cr[getint_forC_2D(xc,nt,tau-r-1,Nt,T,head,offState)];
                    }
                }
            }
            
            // Compute the contribution of ancestors to instant tau
            int minr = tau-t+1;
            int maxr = (L>tau?tau-1:L-1);
            for(int r=maxr; r>=minr; r--) {
                int indaux = tau-t+L-r-1;
                int ii = 0;
                while((ii<N)&&(getint_2D(Act_anc,ii,indaux,N,L-1)>0)) {
                    int n = getint_2D(Act_anc,ii,indaux,N,L-1)-1;
                    for(int nr=0; nr<Nr; nr++) {
                        double valr = 0;
                        double vali = 0;
                        for(int nt=0; nt<Nt; nt++) {
                            valr += getdouble_3D(Hr,nr,nt,r,Nr,Nt,L)*Cr[getint_forC_Xt(Xt,nt,n,tau-r-1,Nt,N,T,head,offState)] \
                                    -getdouble_3D(Hi,nr,nt,r,Nr,Nt,L)*Ci[getint_forC_Xt(Xt,nt,n,tau-r-1,Nt,N,T,head,offState)];
                            vali += getdouble_3D(Hr,nr,nt,r,Nr,Nt,L)*Ci[getint_forC_Xt(Xt,nt,n,tau-r-1,Nt,N,T,head,offState)] \
                                    +getdouble_3D(Hi,nr,nt,r,Nr,Nt,L)*Cr[getint_forC_Xt(Xt,nt,n,tau-r-1,Nt,N,T,head,offState)];
                        }
                        setdouble_2D(valr,ancContribr[indaux],nr,n,Nr,N);
                        setdouble_2D(vali,ancContribi[indaux],nr,n,Nr,N);
                    }
                    if(r<=maxr-1) {
                        int nAnt = getint_2D(a_ind,n,tau-r-1,N,T);
                        for(int nr=0; nr<Nr; nr++) {
                            double valr = getdouble_2D(ancContribr[indaux],nr,n,Nr,N)+getdouble_2D(ancContribr[indaux-1],nr,nAnt,Nr,N);
                            double vali = getdouble_2D(ancContribi[indaux],nr,n,Nr,N)+getdouble_2D(ancContribi[indaux-1],nr,nAnt,Nr,N);
                            setdouble_2D(valr,ancContribr[indaux],nr,n,Nr,N);
                            setdouble_2D(vali,ancContribi[indaux],nr,n,Nr,N);
                        }
                    }
                    ii++;
                }
            }
            
            // For each particle, remove the contribution of ancestors to instant
            // tau and compute the likelihood (add it to logWY)
            for(int n=0; n<N; n++) {
                double val = 0;
                double aux1;
                double aux2;
                for(int nr=0; nr<Nr; nr++) {
                    aux1 = Ypredr[nr] - getdouble_2D(ancContribr[L-2],nr,n,Nr,N);
                    aux2 = Ypredi[nr] - getdouble_2D(ancContribi[L-2],nr,n,Nr,N);
                    val += aux1*aux1+aux2*aux2;
                }
                logWY[n] -= val/sy2;
            }
        }
    }
    else {
        for(int tau=t; tau<=maxtau; tau++) {
            //Ypred = Y(:,tau);
            for(int nr=0; nr<Nr; nr++) {
                Ypredr[nr] = getdouble_2D(Yr,nr,tau-1,Nr,T);
            }
            
            // Remove the contribution of the fixed particle to instant tau
            for(int r=0; r<=tau-t; r++) {
                for(int nr=0; nr<Nr; nr++) {
                    for(int nt=0; nt<Nt; nt++) {
                        Ypredr[nr] -= getdouble_3D(Hr,nr,nt,r,Nr,Nt,L)*Cr[getint_forC_2D(xc,nt,tau-r-1,Nt,T,head,offState)];
                    }
                }
            }
            
            // Compute the contribution of ancestors to instant tau
            int minr = tau-t+1;
            int maxr = (L>tau?tau-1:L-1);
            for(int r=maxr; r>=minr; r--) {
                int indaux = tau-t+L-r-1;
                int ii = 0;
                while((ii<N)&&(getint_2D(Act_anc,ii,indaux,N,L-1)>0)) {
                    int n = getint_2D(Act_anc,ii,indaux,N,L-1)-1;
                    for(int nr=0; nr<Nr; nr++) {
                        double valr = 0;
                        for(int nt=0; nt<Nt; nt++) {
                            valr += getdouble_3D(Hr,nr,nt,r,Nr,Nt,L)*Cr[getint_forC_Xt(Xt,nt,n,tau-r-1,Nt,N,T,head,offState)];
                        }
                        setdouble_2D(valr,ancContribr[indaux],nr,n,Nr,N);
                    }
                    if(r<=maxr-1) {
                        int nAnt = getint_2D(a_ind,n,tau-r-1,N,T);
                        for(int nr=0; nr<Nr; nr++) {
                            double valr = getdouble_2D(ancContribr[indaux],nr,n,Nr,N)+getdouble_2D(ancContribr[indaux-1],nr,nAnt,Nr,N);
                            setdouble_2D(valr,ancContribr[indaux],nr,n,Nr,N);
                        }
                    }
                    ii++;
                }
            }
            
            // For each particle, remove the contribution of ancestors to instant
            // tau and compute the likelihood (add it to logWY)
            for(int n=0; n<N; n++) {
                double val = 0;
                double aux1;
                for(int nr=0; nr<Nr; nr++) {
                    aux1 = Ypredr[nr] - getdouble_2D(ancContribr[L-2],nr,n,Nr,N);
                    val += aux1*aux1;
                }
                logWY[n] -= val/sy2;
            }
        }
    }
    
    // Free memory
    mxDestroyArray(Ypred_mx);
    for(int ll=0; ll<L-1; ll++) {
        mxDestroyArray(ancContrib_mx[ll]);
    }
    free(ancContrib_mx);
    free(ancContribr);
    free(ancContribi);
    free(Act_anc);
    free(Act_ancB);
    
    // Return result
    return logWY;
}


double *computeFactor2(int Nt, int N, int t, int T, int Q, double logQ_1, int *xc, int16_T *Xt, double *logam, double *logbm, double *log1_am, double *log1_bm, int lHead, int onOffModel, int offState) {
    // Allocate memory for result
    double *logWZ = (double*)calloc(N,sizeof(double));
    int prev;
    int curr;
    
    // For each particle, sum the log-transition probability of each transmitter
    for(int n=0; n<N; n++) {
        for(int nt=0; nt<Nt; nt++) {
            prev = getint_Xt(Xt,nt,n,t-1,Nt,N,T);
            curr = getint_2D(xc,nt,t,Nt,T);
            // If it was inactive...
            if(prev==0) {
                // ... and it is still inactive
                if(curr==0) {
                    logWZ[n] += logam[nt];
                // ... and it is active now
                } else {
                    if(lHead==0 && curr>0) {
                        logWZ[n] += log1_am[nt]-logQ_1;
                    } else if(lHead>0 && curr==-1) {
                        logWZ[n] += log1_am[nt];
                    } else {
                        logWZ[n] -= myINFINITY;
                    }
                }
            // If it was active and positive...
            } else if(prev>0) {
                // ... and it is inactive now
                if(curr==0 && !onOffModel) {
                    logWZ[n] += logbm[nt];
                // ... and it is in the offState now
                } else if(curr==offState && onOffModel) {
                    logWZ[n] += logbm[nt];
                // ... and it is still active
                } else if(curr>0) {
                    logWZ[n] += log1_bm[nt]-logQ_1;
                // ... and there has been a forbidden transition
                } else {
                    logWZ[n] -= myINFINITY;
                }
            // If it was active and negative...
            } else if(prev>offState) {
                // ... and it is inactive now
                if(curr==0 && !onOffModel) {
                    logWZ[n] += logbm[nt];
                // ... and it is in the offState now
                } else if(curr==offState && onOffModel) {
                    logWZ[n] += logbm[nt];
                // ... and it is still active
                } else if(curr!=offState && curr!=0) {
                    if(-prev>=lHead && curr>0) {
                        logWZ[n] += log1_bm[nt]-logQ_1;
                    } else if(-prev<lHead && curr==prev-1) {
                        logWZ[n] += log1_bm[nt];
                    } else {
                        logWZ[n] -= myINFINITY;
                    }
                // ... and there has been a forbidden transition
                } else {
                    logWZ[n] -= myINFINITY;
                }
            // If it was in the offState...
            } else {
                // ... and it is still in the offState
                if(curr==offState) {
                    logWZ[n] += 0;
                // ... and there has been a forbidden transition
                } else {
                    logWZ[n] -= myINFINITY;
                }
            }
        }
    }
    
    // Return result
    return logWZ;
}


/************************ Misc. Utility functions ************************/

inline double my_exp(double val) {
    return (val<-700.0?0.0:gsl_sf_exp(val));
}

inline double my_log(double val) {
    return (val<=0?-myINFINITY:gsl_sf_log(val));
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

