function s2H = sample_post_s2H(data,samples,hyper,param)

if(~param.infer.sampleVarH)
    s2H = samples.s2H;
    return;
end

% Obtain parameters from the structs
Nr = param.Nr;
Nt = size(samples.Z,1);

% Prior hyperparameters:
tau0 = 2+1/(hyper.kappa^2);
nu0 = (tau0-1)*hyper.s2h*exp(-hyper.lambda*(0:param.L-1));

% Posterior parameters
nuP = nu0+squeeze(sum(sum(abs(samples.H).^2,1),2)).';
tauP = tau0+Nt*Nr;

% Sample the variance from an Inverse-Gamma distribution
s2H = 1./gamrnd(tauP,1./nuP);
