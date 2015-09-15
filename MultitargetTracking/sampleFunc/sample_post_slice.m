function s = sample_post_slice(data,samples,hyper,param)

c = 1-samples.am;
cmin = min(c);
s = betarnd(param.bnp.betaSlice1,param.bnp.betaSlice2)*cmin;
