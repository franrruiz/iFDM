function llh = compute_llh(data,samples,hyper,param)

% Obtain parameters from the structs
[Nt aux T] = size(samples.Z);
Nr= size(data.obs,1);

Ptot=zeros(Nr,T);
for jj=1:Nr
    d= 1./((repmat(data.sensors(jj,1),Nt,T)-reshape(samples.Z(:,1,:),[Nt,T])).^2 +(repmat(data.sensors(jj,2),Nt,T)-reshape(samples.Z(:,2,:),[Nt,T])).^2).^(param.pathL/2);
    d(squeeze(samples.Z(:,1,:))==0)=0;
    Ptot(jj,:)=   10^(data.Ptx/10)*param.d0^param.pathL*sum(d,1);
end

llh = -T*log(pi)-T*log(samples.s2y)-sum(sum(abs(data.obs-Ptot).^2))/samples.s2y;