function bX = sample_post_s2X(data,samples,hyper,param)

if(~param.infer.sampleVarX)
    bX = samples.bX;
    return;
end

% Obtain parameters from the structs
Nr = param.D;
[Nt T] = size(samples.Z);

% Posterior parameters
nuP = hyper.nubX+0.5*sum(sum(samples.Z~=0));
tauP = hyper.taubX+0.5*sum(sum((samples.Z).^2));

% Sample the variance from an Inverse-Gamma distribution
bX = 1./gamrnd(nuP,1./tauP);
