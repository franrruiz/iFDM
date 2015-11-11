function [Sest SeqEst nest] = sample_Z_FFBS_gauss(data,samples,hyper,param)

%% Extract parameters from structs
Q = 2;
s2y = samples.s2y;
s2x = samples.bX;
T = param.T;

Sest = samples.Z;           % Size = [M x T]
Mest = size(Sest,1);
SeqEst = zeros(Mest,T);     % Size = [M x T]
W = samples.W;              % Size = [M x D]

%% Forward-Filtering Backward-Sampling
fw = zeros(Q,T);
bw = zeros(Q,T);
for it=1:param.ffbs.Niter
    for m=1:Mest
        % Remove the contribution of all transmitters but the m-th one
        suma = zeros(param.D,T);
        Zq = Sest([1:m-1 m+1:end],:);
        suma = suma+W([1:m-1 m+1:end],:)'*Zq;

        Xhat = data.obs-suma;    % Size = [D x T]

        % Build transition matrix A
        A = zeros(Q,Q);
        A(1,1) = samples.am(m);
        A(1,2) = 1-samples.am(m);
        A(2,1) = samples.bm(m);
        A(2,2) = 1-samples.bm(m);
        
        % Forward filtering
        % (i) t=1
        t = 1;
        s2P = 1/((1/s2x)+(sum(W(m,:).^2)/s2y));
        muP = (s2P/s2y)*W(m,:)*Xhat(:,t);
        fw(1,t) = -1/(2*s2y)*sum(Xhat(:,t).^2)+log(A(1,1));
        fw(2,t) = 0.5*log(s2P/s2x)-1/(2*s2y)*sum(Xhat(:,t).^2)+muP^2/(2*s2P)+log(A(1,2));
        fw(:,t) = exp(fw(:,t)-max(fw(:,t)));
        fw(:,t) = fw(:,t)/sum(fw(:,t));
        % (ii) t>1
        for t=2:T
            % s2P doesn't change
            muP = (s2P/s2y)*W(m,:)*Xhat(:,t);
            fw(1,t) = -1/(2*s2y)*sum(Xhat(:,t).^2)...
                      +log(fw(1,t-1)*A(1,1)+fw(2,t-1)*A(2,1));
            fw(2,t) = 0.5*log(s2P/s2x)-1/(2*s2y)*sum(Xhat(:,t).^2)+muP^2/(2*s2P)...
                      +log(fw(1,t-1)*A(1,2)+fw(2,t-1)*A(2,2));
            fw(:,t) = exp(fw(:,t)-max(fw(:,t)));
            fw(:,t) = fw(:,t)/sum(fw(:,t));
        end

        % Backward sampling
        t = T;
        muP = (s2P/s2y)*W(m,:)*Xhat(:,t);
        aux = mnrnd(1,fw(:,t))*(0:Q-1).';
        SeqEst(m,t) = aux;
        Sest(m,t) = aux*(muP+sqrt(s2P)*randn(1));
        for t=T-1:-1:1
            muP = (s2P/s2y)*W(m,:)*Xhat(:,t);
            bw(:,t) = fw(:,t).*A(:,SeqEst(m,t+1)+1);
            bw(:,t) = bw(:,t)/sum(bw(:,t));
            aux = mnrnd(1,bw(:,t))*(0:Q-1).';
            SeqEst(m,t) = aux;
            Sest(m,t) = aux*(muP+sqrt(s2P)*randn(1));
        end
    end
end

%% Update the number of transitions, nest
nest = zeros(2,2,size(Sest,1));
for m=1:Mest
    % From 0 to 0
    nest(1,1,m) = sum([0 SeqEst(m,1:end-1)]==0 & SeqEst(m,:)==0);
    % From 0 to active
    nest(1,2,m) = sum([0 SeqEst(m,1:end-1)]==0 & SeqEst(m,:)~=0);
    % From active to 0
    nest(2,1,m) = sum([0 SeqEst(m,1:end-1)]~=0 & SeqEst(m,:)==0);
    % From active to active
    nest(2,2,m) = sum([0 SeqEst(m,1:end-1)]~=0 & SeqEst(m,:)~=0);
end
