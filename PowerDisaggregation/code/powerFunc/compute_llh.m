function llh = compute_llh(data,samples,hyper,param)

% Obtain parameters from the structs
[Nt T] = size(samples.Z);
Q = param.Q;
Ypred   =   zeros(1,T);
for qq = 1 : Q
     Ypred   =   Ypred   +   samples.P(qq,:)*(samples.Z==qq);
end        
% Evaluate multivariate complex normal distribution
llh = -T*log(pi)-T*log(samples.s2y)-sum(sum(abs(data.obs-Ypred).^2))/samples.s2y;
