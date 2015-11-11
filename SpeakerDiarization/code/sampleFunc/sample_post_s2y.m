function s2y = sample_post_s2y(data,samples,hyper,param)

if(~param.infer.sampleNoiseVar)
    s2y = samples.s2y;
    return;
end

% Obtain parameters from the structs
Nr = param.D;
[Nt T] = size(samples.Z);

% Posterior parameters
nuP = hyper.nu+T*Nr/2;
tauP = hyper.tau +1/2*sum(sum(abs(data.obs-samples.W'*samples.Z).^2));

% Sample the noise variance from an Inverse-Gamma distribution
s2y = 1/gamrnd(nuP,1/tauP);

