function s2y = sample_post_s2y(data,samples,hyper,param)

if(~param.infer.sampleNoiseVar)
    s2y = samples.s2y;
    return;
end

% Obtain parameters from the structs
Nr = param.D;
[Nt aux T] = size(samples.Z);

% Build matrix S containing the symbols and their shifted replicas,
% and matrix H containing the channel coefficients
Ptot=zeros(Nr,T);
for jj=1:Nr
    d= 1./((repmat(data.sensors(jj,1),Nt,T)-reshape(samples.Z(:,1,:),[Nt,T])).^2 +(repmat(data.sensors(jj,2),Nt,T)-reshape(samples.Z(:,2,:),[Nt,T])).^2).^(param.pathL/2);
    d(squeeze(samples.Z(:,1,:))==0)=0;
    Ptot(jj,:)=   10^(data.Ptx/10)*param.d0^param.pathL*sum(d,1);
end

% Posterior parameters
nuP = hyper.nu+T*Nr;
tauP = hyper.tau+sum(sum(abs(data.obs-Ptot).^2));
% [Note: There is a typo in the ModelNotes.pdf file (in that file, there is
%  an extra 0.5 factor for tauP). Here it has been fixed]

% Sample the noise variance from an Inverse-Gamma distribution
s2y = 1/gamrnd(nuP,1/tauP);

