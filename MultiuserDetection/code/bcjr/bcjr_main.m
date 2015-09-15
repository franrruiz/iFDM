function [Sest qt_red Simb_red] = bcjr_main(data,samples,hyper,param)

%% Extract parameters from param
p1 = param.bcjr.p1;
p2 = param.bcjr.p2;
Nt = size(samples.H,2);

%% Extend the constellation with 0 if desired
if(param.flag0)
    Q = length(param.constellation)+1;
    constellation = [0 param.constellation];
else
    Q = length(param.constellation);
    constellation = param.constellation;
end

%% Build constellationT for a single transmitter. Size = [Q^param.L,param.L]
constellationT = zeros(Q^param.L,param.L);
for l=1:param.L
    aux = [];
    for s=1:length(constellation)
        aux = [aux; repmat(constellation(s),Q^(param.L-l),1)];
    end
    constellationT(:,l) = repmat(aux,Q^(l-1),1);
end

%% Build transition probabilities (for one transmitter)
ptrans = zeros(Q^param.L,Q^param.L);
for i=1:Q^param.L
    currentSt = constellationT(i,:);
    for q=1:Q
        newSt = [constellation(q) currentSt(1:end-1)];
        [valnul idx] = min(sum(abs(constellationT-repmat(newSt,Q^param.L,1)),2));
        if(constellation(q)==0 && currentSt(1)==0)
            ptrans(i,idx) = p1;
        elseif(constellation(q)~=0 && currentSt(1)==0)
            ptrans(i,idx) = (1-p1)/(Q-1);
        elseif(constellation(q)==0 && currentSt(1)~=0)
            ptrans(i,idx) = p2;
        elseif(constellation(q)~=0 && currentSt(1)~=0)
            if(param.flag0)
                ptrans(i,idx) = (1-p2)/(Q-1);
            else
                ptrans(i,idx) = (1-p2)/Q;
            end
        end
    end
end

%% Build transition probabilities for all transmitters
ptransAux = ptrans;
for i=1:Nt-1
    ptrans = kron(ptransAux,ptrans);
end
Qtot = size(ptrans,1);
ptrans=ptrans./(repmat(sum(ptrans,2),1,Qtot));

%% Build matrix A with symbols for each state. Size=[#states(total),  Nt*L]
A = zeros(Qtot,Nt*param.L);
if Nt==2
    for i=1:Q^param.L
        for i2=1:Q^param.L
            A((i-1)*Q^param.L+i2,:)=[constellationT(i,:) constellationT(i2,:)];
        end
    end
elseif Nt==3
    for i=1:Q^param.L
        for i2=1:Q^param.L
            for i3=1:Q^param.L
                A((i-1)*Q^(2*param.L)+(i2-1)*Q^param.L+i3,:)=[constellationT(i,:) constellationT(i2,:) constellationT(i3,:)];
            end
        end
    end
elseif Nt==4
    for i=1:Q^param.L
        for i2=1:Q^param.L
            for i3=1:Q^param.L
                for i4=1:Q^param.L
                    A((i-1)*Q^(3*param.L)+(i2-1)*Q^(2*param.L)+(i3-1)*Q^param.L+i4,:)=[constellationT(i,:) constellationT(i2,:) constellationT(i3,:) constellationT(i4,:)];
                end
            end
        end
    end
else
    error('Not implemented: too many transmitters');
end

%% Build the auxiliary channel. Size = [Nr,  Nt*L]
Haux = zeros(param.Nr,param.L*Nt);
for tx=1:Nt
    for m=1:param.L
        Haux(:,(tx-1)*param.L+m) = samples.H(:,tx,m);
    end
end

%% Run BCJR
qt = bcjr(data.obs,Haux,samples.s2y,ptrans,A);

%% Probabilites qt_red (reduced to only consider time instant t), i.e., marginalize t-1, t-2, ..., t-L+1
qt_red = zeros(Q^Nt,param.T);
Simb_red = zeros(Q^Nt,Nt);
for i=1:Q^Nt
    aux = de2bi(i-1,Nt,Q)+1;
    Simb_red(i,:) = constellation(aux);
    indices = find(sum(repmat(constellation(aux),Qtot,1)==A(:,1:param.L:Nt*param.L),2)==Nt);
    for t=1:param.T
        qt_red(i,t) = sum(qt(indices,t));
    end
end

%% Estimated symbols. Size = [Nt, T]
[val idx] = max(qt_red,[],1);
Sest = zeros(Nt,param.T);
for t=1:param.T
    Sest(:,t) = Simb_red(idx(t),:).';
end
