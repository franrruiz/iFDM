function [am bm] = sample_post_transitionProb(data,samples,hyper,param)

M = size(samples.Z,1);
am = zeros(M,1);
bm = zeros(M,1);
for m=1:M
    % Sample am
    am(m) = 1-betarnd(samples.nest(1,2,m),samples.nest(1,1,m)+1);
    % Sample bm
    bm(m) = betarnd(hyper.gamma1+samples.nest(2,1,m),hyper.gamma2+samples.nest(2,2,m));
end
