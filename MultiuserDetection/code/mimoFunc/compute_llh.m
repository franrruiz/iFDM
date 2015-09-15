function llh = compute_llh(data,samples,hyper,param)

% Obtain parameters from the structs
Nr = param.Nr;
[Nt T] = size(samples.Z);
L = param.L;

% Build matrix S containing the symbols and their shifted replicas,
% and matrix H containing the channel coefficients
S = zeros(T,Nt*L);
H = zeros(Nt*L,Nr);
for ll=1:L
    S(ll:T,ll:L:L*Nt) = samples.Z(:,1:T-ll+1).';
    H(ll:L:L*Nt,:) = samples.H(:,:,ll).';
end

% Evaluate multivariate complex normal distribution
llh = -T*Nr*log(pi)-T*Nr*log(samples.s2y)-sum(sum(abs(data.obs.'-S*H).^2))/samples.s2y;
