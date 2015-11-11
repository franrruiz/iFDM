function s2W = sample_post_s2W(data,samples,hyper,param)

if(~param.infer.sampleWVar)
    s2W = samples.s2W;
    return;
end

% Obtain parameters from the structs
Nr = param.D;
[Nt T] = size(samples.Z);

% Posterior parameters
nuP = hyper.nuW+Nt*Nr/2;
tauP = hyper.tauW+0.5*sum(sum(abs(samples.W).^2));

% Sample the variance from an Inverse-Gamma distribution
s2W = 1./gamrnd(nuP,1./tauP);
