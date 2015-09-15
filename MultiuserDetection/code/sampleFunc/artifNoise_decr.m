function [obs s2y] = artifNoise_decr(data,samples,hyper,param)

if(~param.infer.addArtificialNoise)
    obs = data.obs;
    s2y = samples.s2y;
    return;
end

% If the current SNR is the highest, do nothing
curSNR = -10*log10(samples.s2y);
if(curSNR>=param.artifNoise.finalSNR)
    s2y = samples.s2y;
    obs = data.obsWithoutNoise;
    return;
end

% Else, decrease the variance of the artificial noise
if(curSNR+param.artifNoise.stepDB>param.artifNoise.finalSNR)
    % If currentSNR+stepSNR>finalSNR, then remove all artificial noise
    targetVar = 10^(-param.artifNoise.finalSNR/10);
    extraVar = 0;
else
    % If currentSNR+stepSNR<finalSNR, then just decrease the variance of the artificial noise
    targetVar = 10^(-(curSNR+param.artifNoise.stepDB)/10);
    extraVar = targetVar-10^(-param.artifNoise.finalSNR/10);
end

obs = data.obsWithoutNoise+sqrt(extraVar)*data.artifNoise;
s2y = targetVar;
