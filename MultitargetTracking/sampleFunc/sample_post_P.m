function [P MuP] = sample_post_H(data,samples,hyper,param)

% Obtain parameters from the structs
[Nt T] = size(samples.Z);
Q=param.Q;

S=[];
for qq=1:Q
    S=[S ; samples.Z==qq];
end

PrecP=((1/samples.s2y)* (S*S')+ 1/hyper.s2P*diag(ones(Nt*Q,1)));
MuP= PrecP\(1/samples.s2y*(S*data.obs') +hyper.muP/hyper.s2P*ones(Nt*Q,1));

Paux = mvnrnd2(MuP,inv(PrecP),1);


P=reshape(Paux,Nt,Q)';



