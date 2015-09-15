function [Sest SeqEst nest] = sample_Z_FFBS(data,samples,hyper,param)

auxConstellation = param.constellation;
if(param.flag0)
    auxConstellation = [0 param.constellation];
end
Q = length(auxConstellation)^(param.L);

Sest = samples.Z;

Mest = size(Sest,1);
A = zeros(Q,Q,Mest);
A(1,1,:) = samples.am;
for m=1:Mest
A(1,2:end,m) = (1-samples.am(m))/(Q-1);
A(2:end,1,m) = samples.bm(m);
A(2:end,2:end,m) = (1-samples.bm(m))/(Q-1);
end
Sest=samples.Z;
Phi(1,:,:)=samples.P;
%% Forward-Filtering Backward-Sampling
for it=1:param.ffbs.Niter    
    Sest=FB(data.obs,samples.s2y,Sest, Phi, A);
end
SeqEst  =Sest;       
%% Update the number of transitions, nest
nest = zeros(2,2,size(Sest,1));
for m=1:Mest
    % From 0 to 0
    nest(1,1,m) = sum([0 SeqEst(m,1:end-1)]==0 & SeqEst(m,:)==0);
    % From 0 to active
    nest(1,2,m) = sum([0 SeqEst(m,1:end-1)]==0 & SeqEst(m,:)~=0);
    % From active to 0
    nest(2,1,m) = sum([0 SeqEst(m,1:end-1)]~=0 & SeqEst(m,:)==0);
    % From active to active
    nest(2,2,m) = sum([0 SeqEst(m,1:end-1)]~=0 & SeqEst(m,:)~=0);
end
