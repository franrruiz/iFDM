function [W] = sample_post_W(data,samples,hyper,param)

% Obtain parameters from the structs
[Nt T] = size(samples.Z);
D=param.D;
W=zeros(Nt,D);
ZZ=samples.Z*samples.Z';
for d=1:D
    PrecW=((1/samples.s2y)* ZZ+ 1/samples.s2W*diag(ones(Nt,1)));
    MuW= PrecW\(1/samples.s2y*(samples.Z*data.obs(d,:)'));
    W(:,d) = mvnrnd2(MuW,inv(PrecW),1);
end




