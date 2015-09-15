function [Sest SeqEst nest] = sample_Z_FFBS(data,samples,hyper,param)

auxConstellation = param.constellation;
if(param.flag0)
    auxConstellation = [0 param.constellation];
end
Q = length(auxConstellation)^(param.L);

Sest = samples.Z;
SeqEst = samples.seq;

Mest = size(Sest,1);

%% Forward-Filtering Backward-Sampling
fw = zeros(Q,param.T);
bw = zeros(Q,param.T);
for it=1:param.ffbs.Niter
    for m=1:Mest
        % Extended state variable for m-th chain:
        extZ = zeros(1,param.T);

        % Remove the contribution of all transmitters but the m-th one
        suma = zeros(param.Nr,param.T);
        for ll=0:param.L-1
            Zq = Sest([1:m-1 m+1:end],:);
            Zq = cat(2,zeros(size(Zq,1),ll),Zq(:,1:(param.T-ll)));
            suma = suma+samples.H(:,[1:m-1 m+1:end],ll+1)*Zq;
        end
        Xhat = data.obs-suma;

        % Build matrix Phi and transition matrix A
        Phi = zeros(param.Nr,Q);
        A = sparse(Q,Q);
        for q=1:Q
            % Current state
            currSt = de2bi(q-1,param.L,length(auxConstellation));
            Phi(:,q) = squeeze(samples.H(:,m,:))*(auxConstellation(1+currSt).');
            % Allowed next states
            for r=0:length(auxConstellation)-1
                nextSt = [r currSt(1:end-1)];
                qq = 1+bi2de(nextSt,length(auxConstellation));
                if(param.flag0)
                    if(r==0 && currSt(1)==0)
                        A(q,qq) = samples.am(m);
                    elseif(r==0 && currSt(1)~=0)
                        A(q,qq) = samples.bm(m);
                    elseif(r~=0 && currSt(1)==0)
                        A(q,qq) = (1-samples.am(m))/(length(auxConstellation)-1);
                    elseif(r~=0 && currSt(1)~=0)
                        A(q,qq) = (1-samples.bm(m))/(length(auxConstellation)-1);
                    end
                else
                    A(q,qq) = 1/length(auxConstellation);
                end
            end
        end

        % Forward
        t = 1;
        fw(:,t) = -(1/samples.s2y)*sum(abs(repmat(Xhat(:,t),1,Q)-Phi).^2,1).'+log(A(1,:).');
        fw(:,t) = exp(fw(:,t)-max(fw(:,t)));
        fw(:,t) = fw(:,t)/sum(fw(:,t));
        for t=2:param.T
            fw(:,t) = -(1/samples.s2y)*sum(abs(repmat(Xhat(:,t),1,Q)-Phi).^2,1).'+log(A.'*fw(:,t-1));
            fw(:,t) = exp(fw(:,t)-max(fw(:,t)));
            fw(:,t) = fw(:,t)/sum(fw(:,t));
        end

        % Backward
        t = param.T;
        extZ(t) = mnrnd(1,fw(:,t))*(0:Q-1).';
        aux = 1+de2bi(extZ(t),param.L,length(auxConstellation));
        SeqEst(m,t) = aux(1)-param.flag0;
        Sest(m,t) = auxConstellation(SeqEst(m,t)+param.flag0);
        for t=param.T-1:-1:1
            bw(:,t) = fw(:,t).*A(:,extZ(t+1)+1); 
            bw(:,t) = bw(:,t)/sum(bw(:,t));
            extZ(t) = mnrnd(1,bw(:,t))*(0:Q-1).';
            aux = 1+de2bi(extZ(t),param.L,length(auxConstellation));
            SeqEst(m,t) = aux(1)-param.flag0;
            Sest(m,t) = auxConstellation(SeqEst(m,t)+param.flag0);
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
