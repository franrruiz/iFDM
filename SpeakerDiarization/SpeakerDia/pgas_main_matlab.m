function [Sest SeqEst] = pgas_main_matlab(data,samples,hyper,param)

%% Extract parameters from structs
Nt = size(samples.W,1);
M = param.pgas.Niter;
N_PF = param.pgas.N_PF;
N_PG = param.pgas.N_PG;
returnN = param.pgas.returnNsamples;

Q = param.Q;
if(param.flag0)
    constellation = [0 1:Q];
end

%% Call the pgas function
X_PG = pgas(data.obs,samples.Z,Nt,constellation,Q,param.L,samples.W,samples.s2y,samples.am,samples.bm,N_PF,N_PG,M);

%% Return last obtained samples in a [Nt x T x M] matrix
SeqEst = permute(X_PG(:,M-returnN+1:M,:),[1 3 2]);
Sest = SeqEst;

