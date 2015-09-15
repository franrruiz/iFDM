function s2y = sample_post_s2y(data,samples,hyper,param)

if(~param.infer.sampleNoiseVar)
    s2y = samples.s2y;
    return;
end

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

% Posterior parameters
nuP = hyper.nu+T*Nr;
tauP = hyper.tau+sum(sum(abs(data.obs.'-S*H).^2));
% [Note: There is a typo in the ModelNotes.pdf file (in that file, there is
%  an extra 0.5 factor for tauP). Here it has been fixed]

% Sample the noise variance from an Inverse-Gamma distribution
s2y = 1/gamrnd(nuP,1/tauP);

