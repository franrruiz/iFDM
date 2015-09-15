function [ptrans] = sample_post_ptrans(data,samples,hyper,param)

% Obtain parameters from the structs
[Nt T] = size(samples.Z);
Q=param.Q;
ptrans=zeros(Q+1,Q,Nt);

S=samples.Z();
Sd=[ones(Nt,1) samples.Z(:,1:T-1)];
nqq=zeros(Q+1,Q,Nt);
for qq=0:Q
    for qq2=1:Q
        nqq(1+qq,qq2,:)=sum(Sd==qq & S==qq2,2);
    end
    
    for mm=1:Nt
        ptrans(qq+1,:,mm) = dirichletrnd(hyper.gamma*ones(1,Q)+nqq(qq+1,:,mm), 1);
    end

end


