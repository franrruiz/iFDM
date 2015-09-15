#include "acc_loop.h"

#define input_Nt prhs[0]
#define input_Nr prhs[1]
#define input_N prhs[2]
#define input_T prhs[3]
#define input_L prhs[4]
#define input_t prhs[5]
#define input_Q prhs[6]
#define input_C prhs[7]
#define input_Y prhs[8]
#define input_Xt prhs[9]
#define input_xc prhs[10]
#define input_H prhs[11]
#define input_a_ind prhs[12]
#define input_Act_anc prhs[13]
#define input_sy2 prhs[14]
#define output_logWY plhs[0]

/* INPUTS:
 *
 *  0: Nt,      #transmitters
 *  1: Nr,      #receivers
 *  2: N,       #particles
 *  3: T,       #time steps
 *  4: L,       memory
 *  5: t,       current time instant
 *  6: Q,       constellation order (inc. 0)
 *  7: C,       constellation (1 x Q)
 *  8: Y,       observations (Nt x T)
 *  9: Xt,      particles (Nt x N x T)
 * 10: xc,      particle we condition on (Nt x T)
 * 11: H,       channel (Nr x Nt x L)
 * 12: a_ind,   indexes of the ancestors (N x T)
 * 13: Act_anc, active ancestors in (t-L+1,...,t-2) (N x L-1)
 * 14: sy2,     noise variance
 *
 */

/* OUTPUTS:
 *
 *  0: logWY,    log-product in Eq. (5) (N x 1)
 *
 */


void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] ) {
    
    /* Read inputs */
    int Nt = mxGetScalar(input_Nt);
    int Nr = mxGetScalar(input_Nr);
    int N = mxGetScalar(input_N);
    int T = mxGetScalar(input_T);
    int L = mxGetScalar(input_L);
    int t = mxGetScalar(input_t);
    int Q = mxGetScalar(input_Q);
    double sy2 = mxGetScalar(input_sy2);
    
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
    
    double *Xt = mxGetPr(input_Xt);
    double *xc = mxGetPr(input_xc);
    
    double *a_ind = mxGetPr(input_a_ind);
    double *Act_anc = mxGetPr(input_Act_anc);
    
    /* Allocate memory for output */
    output_logWY = mxCreateDoubleMatrix(N,1,mxREAL);
    double *logWY = mxGetPr(output_logWY);
    
    /* Allocate memory for auxiliary variables */
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
    
    /* LOOP */
    //for tau=t:min(t+L-2,T)
    int maxtau = (t+L-2>T?T:t+L-2);
    if(flagComplex) {
        for(int tau=t; tau<=maxtau; tau++) {
            //Ypred = Y(:,tau);   % Ypred is initialized to the tau-th observation
            for(int nr=0; nr<Nr; nr++) {
                Ypredr[nr] = get_Y(Yr,nr,tau-1,Nr,T);
                Ypredi[nr] = get_Y(Yi,nr,tau-1,Nr,T);
            }
            
            // Remove the contribution of the fixed particle to instant tau
            //  for r=0:tau-t
            //      Ypred = Ypred-H(:,:,r+1)*C(xc(:,tau-r)+1).';
            //  end
            for(int r=0; r<=tau-t; r++) {
                for(int nr=0; nr<Nr; nr++) {
                    for(int nt=0; nt<Nt; nt++) {
                        Ypredr[nr] -= get_H(Hr,nr,nt,r,Nr,Nt,L)*Cr[getint_2D(xc,nt,tau-r-1,Nt,T)] \
                                      -get_H(Hi,nr,nt,r,Nr,Nt,L)*Ci[getint_2D(xc,nt,tau-r-1,Nt,T)];
                        Ypredi[nr] -= get_H(Hr,nr,nt,r,Nr,Nt,L)*Ci[getint_2D(xc,nt,tau-r-1,Nt,T)] \
                                      +get_H(Hi,nr,nt,r,Nr,Nt,L)*Cr[getint_2D(xc,nt,tau-r-1,Nt,T)];
                    }
                }
            }
            
            // Compute the contribution of ancestors to instant tau
            //  minr = tau-t+1;
            //  maxr = min(L-1,tau-1);
            int minr = tau-t+1;
            int maxr = (L>tau?tau-1:L-1);
            // for r=maxr:-1:minr
            //     ii = 1;
            //     while((ii<=N)&&(Act_anc(ii,tau-t+L-r)>0))
            //         n = Act_anc(ii,tau-t+L-r);
            //         ancContrib(:,n,tau-t+L-r) = H(:,:,r+1)*C(1+Xt(:,n,tau-r)).';
            //         if(r<=maxr-1)
            //             nAnt = a_ind(n,tau-r);
            //             ancContrib(:,n,tau-t+L-r) = ancContrib(:,n,tau-t+L-r)+ancContrib(:,nAnt,tau-t+L-r-1);
            //         end
            //         ii = ii+1;
            //     end
            // end
            for(int r=maxr; r>=minr; r--) {
                int ii = 0;
                int indaux = tau-t+L-r-1;
                while((ii<N)&&(getint_2D(Act_anc,ii,indaux,N,L-1)>0)){
                    int n = getint_2D(Act_anc,ii,indaux,N,L-1)-1;
                    for(int nr=0; nr<Nr; nr++) {
                        double valr = 0;
                        double vali = 0;
                        for(int nt=0; nt<Nt; nt++) {
                            valr += get_H(Hr,nr,nt,r,Nr,Nt,L)*Cr[getint_3D(Xt,nt,n,tau-r-1,Nt,N,T)] \
                                    -get_H(Hi,nr,nt,r,Nr,Nt,L)*Ci[getint_3D(Xt,nt,n,tau-r-1,Nt,N,T)];
                            vali += get_H(Hr,nr,nt,r,Nr,Nt,L)*Ci[getint_3D(Xt,nt,n,tau-r-1,Nt,N,T)] \
                                    +get_H(Hi,nr,nt,r,Nr,Nt,L)*Cr[getint_3D(Xt,nt,n,tau-r-1,Nt,N,T)];
                        }
                        setdouble_2D(valr,ancContribr[indaux],nr,n,Nr,N);
                        setdouble_2D(vali,ancContribi[indaux],nr,n,Nr,N);
                    }
                    if(r<=maxr-1) {
                        int nAnt = getint_2D(a_ind,n,tau-r-1,N,T)-1;
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
            //  for n=1:N
            //      logWY2(n) = logWY2(n)-sum(abs(Ypred-ancContrib(:,n,L-1)).^2)/sy2;
            //  end
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
            //Ypred = Y(:,tau);   % Ypred is initialized to the tau-th observation
            for(int nr=0; nr<Nr; nr++) {
                Ypredr[nr] = get_Y(Yr,nr,tau-1,Nr,T);
            }
            
            // Remove the contribution of the fixed particle to instant tau
            //  for r=0:tau-t
            //      Ypred = Ypred-H(:,:,r+1)*C(xc(:,tau-r)+1).';
            //  end
            for(int r=0; r<=tau-t; r++) {
                for(int nr=0; nr<Nr; nr++) {
                    for(int nt=0; nt<Nt; nt++) {
                        Ypredr[nr] -= get_H(Hr,nr,nt,r,Nr,Nt,L)*Cr[getint_2D(xc,nt,tau-r-1,Nt,T)];
                    }
                }
            }
            
            // Compute the contribution of ancestors to instant tau
            //  minr = tau-t+1;
            //  maxr = min(L-1,tau-1);
            int minr = tau-t+1;
            int maxr = (L>tau?tau-1:L-1);
            // for r=maxr:-1:minr
            //     ii = 1;
            //     while((ii<=N)&&(Act_anc(ii,tau-t+L-r)>0))
            //         n = Act_anc(ii,tau-t+L-r);
            //         ancContrib(:,n,tau-t+L-r) = H(:,:,r+1)*C(1+Xt(:,n,tau-r)).';
            //         if(r<=maxr-1)
            //             nAnt = a_ind(n,tau-r);
            //             ancContrib(:,n,tau-t+L-r) = ancContrib(:,n,tau-t+L-r)+ancContrib(:,nAnt,tau-t+L-r-1);
            //         end
            //         ii = ii+1;
            //     end
            // end
            for(int r=maxr; r>=minr; r--) {
                int ii = 0;
                int indaux = tau-t+L-r-1;
                while((ii<N)&&(getint_2D(Act_anc,ii,indaux,N,L-1)>0)){
                    int n = getint_2D(Act_anc,ii,indaux,N,L-1)-1;
                    for(int nr=0; nr<Nr; nr++) {
                        double valr = 0;
                        for(int nt=0; nt<Nt; nt++) {
                            valr += get_H(Hr,nr,nt,r,Nr,Nt,L)*Cr[getint_3D(Xt,nt,n,tau-r-1,Nt,N,T)];
                        }
                        setdouble_2D(valr,ancContribr[indaux],nr,n,Nr,N);
                    }
                    if(r<=maxr-1) {
                        int nAnt = getint_2D(a_ind,n,tau-r-1,N,T)-1;
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
            //  for n=1:N
            //      logWY2(n) = logWY2(n)-sum(abs(Ypred-ancContrib(:,n,L-1)).^2)/sy2;
            //  end
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
    
    /* Free memory */
    mxDestroyArray(Ypred_mx);
    for(int ll=0; ll<L-1; ll++) {
        mxDestroyArray(ancContrib_mx[ll]);
    }
    free(ancContrib_mx);
    free(ancContribr);
    free(ancContribi);
}

inline double get_Y(double *Y, int n, int t, int Nr, int T) {
    return Y[Nr*t+n];
}

inline double get_H(double *H, int nr, int nt, int ll, int Nr, int Nt, int L) {
    return H[Nr*Nt*ll+Nr*nt+nr];
}

inline int getint_2D(double *x, int n1, int n2, int N1, int N2) {
    return (int)(x[N1*n2+n1]);
}

inline int getint_3D(double *x, int n1, int n2, int n3, int N1, int N2, int N3) {
    return (int)(x[N1*N2*n3+N1*n2+n1]);
}

inline void setdouble_2D(double val, double *x, int n1, int n2, int N1, int N2) {
    x[N1*n2+n1] = val;
}

inline double getdouble_2D(double *x, int n1, int n2, int N1, int N2) {
    return x[N1*n2+n1];
}



