function [Sest SeqEst] = pgas_main_C(data,samples,hyper,param)

%% Extract parameters from structs
Nt = size(samples.P,2);
M = param.pgas.Niter;
N_PF = param.pgas.N_PF;
N_PG = param.pgas.N_PG;
returnN = param.pgas.returnNsamples;

%% Build auxiliary variables
auxPtrans = zeros(param.Q+1,param.Q+1,Nt);
auxPtrans(1,1,:) = samples.am;
auxPtrans(1,2:param.Q+1,:) = repmat((1-permute(samples.am,[2 3 1])),[1,param.Q,1]).*samples.ptrans(1,:,:);
auxPtrans(2:param.Q+1,1,:) = repmat((permute(samples.bm,[2 3 1])),[param.Q,1,1]);
auxPtrans(2:param.Q+1,2:param.Q+1,:) = repmat((1-permute(samples.bm,[2 3 1])),[param.Q,param.Q,1]).*samples.ptrans(2:param.Q+1,:,:);
auxP = zeros(param.Q+1,Nt);
auxP(2:param.Q+1,:) = samples.P;

%% Call the pgas function
X_PG = pgas_C_parallel(Nt,size(data.obs,1),N_PF,N_PG,M,size(data.obs,2),param.Q,data.obs,1,samples.Z,samples.s2y,auxPtrans,param.pgas.particles,auxP);

%% Return last obtained samples in a [Nt x T x M] matrix
SeqEst = permute(X_PG(:,M-returnN+1:M,:),[1 3 2]);
Sest = SeqEst;

