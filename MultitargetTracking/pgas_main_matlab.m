function [Sest SeqEst] = pgas_main_matlab(data,samples,hyper,param)

%% Extract parameters from structs
Nt = size(samples.Z,1);
M = param.pgas.Niter;
N_PF = param.pgas.N_PF;
N_PG = param.pgas.N_PG;
returnN = param.pgas.returnNsamples;

%% Call the pgas function
X_PG = pgas(data.obs,data.sensors, data.W, data.Gx, data.Gu, param.Ts, param.d0, param.pathL, data.Ptx,data.s2u,samples.Z,Nt,param.L,samples.s2y,samples.am,samples.bm,N_PF,N_PG,M,data.s2vIni);

%% Return last obtained samples in a [Nt x T x M] matrix
SeqEst = permute(X_PG(:,M-returnN+1:M,:,:,:),[1 3 4 2]);
Sest = SeqEst;

