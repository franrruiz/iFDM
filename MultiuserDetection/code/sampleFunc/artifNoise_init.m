function [obsWithoutNoise obs artifNoise s2y] = artifNoise_init(data,samples,hyper,param)

if(~param.infer.addArtificialNoise)
    obs = data.obs;
    obsWithoutNoise = NaN;
    artifNoise = NaN;
    s2y = samples.s2y;
    return;
end

% Save the original observations
obsWithoutNoise = data.obs;

% Create artificial noise
artifNoise = randn(size(data.obs))+1i*randn(size(data.obs));
if(sum(abs(real(param.constellation)))>1e3*sum(abs(imag(param.constellation))))
    artifNoise = real(artifNoise);
else
    artifNoise = artifNoise/sqrt(2);
end

% If the current SNR is too bad, do nothing
curSNR = -10*log10(samples.s2y);
if(curSNR<param.artifNoise.iniSNR)
    s2y = samples.s2y;
    obs = obsWithoutNoise;
    return;
end

% Else, add artificial noise
targetVar = 10^(-param.artifNoise.iniSNR/10);
curVar = samples.s2y;
extraVar = targetVar-curVar;

obs = obsWithoutNoise+sqrt(extraVar)*artifNoise;
s2y = targetVar;
