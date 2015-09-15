function llh = mLik_MIMO(data,Sest,s2y,hyper,param)

% Obtain parameters from the structs
Nr = param.Nr;
[Nt T] = size(Sest);
L = param.L;

% Build matrix S containing the symbols and their shifted replicas
S = zeros(T,Nt*L);
for ll=1:L
    S(ll:T,ll:L:L*Nt) = Sest(:,1:T-ll+1).';
end

% Build the posterior covariance matrix of H
v_1 = (1/hyper.s2h)*exp(hyper.lambda*(0:L-1));  % These are the inverse prior variances
v_1 = repmat(v_1,1,Nt);

GamH = (diag(v_1)+(1/s2y)*(S'*S))\eye(Nt*L);  % Posterior covariance

% Build the posterior precision matrix of Y
PrecY = (eye(T)/s2y-S*GamH'*S'/(s2y^2));

% Evaluate the complex normal distribution
llh = real(-0.5*sum(sum((conj(data.obs)*PrecY).*data.obs,2)) ...
           -0.5*sum(sum((data.obs*conj(PrecY)).*conj(data.obs),2)) ...
           -T*log(pi)+real(logdet(PrecY)));
