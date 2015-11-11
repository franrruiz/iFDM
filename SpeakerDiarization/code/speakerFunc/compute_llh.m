function llh = compute_llh(data,samples,hyper,param)

% Obtain parameters from the structs
D = param.D;
[Nt T] = size(samples.Z);

Ypred   =   zeros(1,T);
% Evaluate multivariate normal distribution
llh = -T*D*log(pi)-T*D*log(samples.s2y)-sum(sum(abs(data.obs-samples.W'*samples.Z).^2))/samples.s2y;
