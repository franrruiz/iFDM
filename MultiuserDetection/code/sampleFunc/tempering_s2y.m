function s2y = tempering_s2y(data,samples,hyper,param)

temperNoiseCur = samples.s2y;
temperIdxCur = find(param.temper.s2yValues==temperNoiseCur);
if((temperIdxCur<length(param.temper.s2yValues)) && (temperIdxCur>1))
    if(rand()<param.temper.pNext)
        % Decrease temperature
        temperNoiseNext = param.temper.s2yValues(temperIdxCur+1);
        auxlogp2 = log(param.temper.pNext);
        auxlogp1 = log(1-param.temper.pNext);
    else
        % Increase temperature
        temperNoiseNext = param.temper.s2yValues(temperIdxCur-1);
        auxlogp2 = log(1-param.temper.pNext);
        auxlogp1 = log(param.temper.pNext);
    end
elseif(temperIdxCur==length(param.temper.s2yValues))
    temperNoiseNext = param.temper.s2yValues(temperIdxCur-1);
    auxlogp2 = 0;
    auxlogp1 = log(param.temper.pNext);
elseif(temperIdxCur==1)
    temperNoiseNext = param.temper.s2yValues(temperIdxCur+1);
    auxlogp2 = 0;
    auxlogp1 = log(1-param.temper.pNext);
end
samples.s2y = temperNoiseCur;
auxLLH2 = compute_llh(data,samples,hyper,param);
samples.s2y = temperNoiseNext;
auxLLH1 = compute_llh(data,samples,hyper,param);
if(rand()<exp(auxlogp1-auxlogp2+auxLLH1-auxLLH2))
    % Accept
    s2y = temperNoiseNext;
else
    % Reject
    s2y = temperNoiseCur;
end