function bX = sample_post_bX(data,samples,hyper,param)

if(~param.infer.sampleVarX)
    bX = samples.bX;
    return;
end

% Obtain parameters from the structs
Nr = param.D;
[Nt T] = size(samples.Z);

% Posterior parameters
nuP = hyper.nubX+sum(sum(samples.Z~=0));
tauP = hyper.taubX+sum(sum((samples.Z~=0).*abs(samples.Z-[zeros(Nt,1) samples.Z(:,1:end-1)])));

% Sample the variance from an Inverse-Gamma distribution
bX = 1./gamrnd(nuP,1./tauP);
