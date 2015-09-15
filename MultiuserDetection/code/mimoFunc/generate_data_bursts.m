function data = generate_data_bursts(param)

% Generate the channel:
data.channel = sqrt(param.gen.varH/2)*randn(param.Nr,param.gen.Nt,max(param.gen.L_true))+1i*sqrt(param.gen.varH/2)*randn(param.Nr,param.gen.Nt,max(param.gen.L_true));
for i=1:param.gen.Nt
    data.channel(:,i,param.gen.L_true(i)+1:max(param.gen.L_true)) = 0;
end
auxSparse = (rand([param.Nr,param.gen.Nt,max(param.gen.L_true)-1])>=param.gen.sparsityH);
data.channel(:,:,2:end) = data.channel(:,:,2:end).*auxSparse;

% Generate the sequence of symbols:
C = length(param.constellation);
idxCpositiveOrthant = find(real(param.constellation)>=0 & imag(param.constellation)>=0);
data.seq = randint(param.gen.Nt,param.T,[1 C]);
data.symbols = zeros(param.gen.Nt,param.T);
data.symbols(:) = param.constellation(data.seq(:));

% Fill symbols with 0's to simulate bursts
Tini_vec = zeros(1,param.gen.Nt);
for i=1:param.gen.Nt
    Tini = 1+randint(1,1,[1 round(param.T/2)]);
    Tini_vec(i) = Tini;
    if(isinf(param.gen.burstLengthStdFactor))
        Tend = min(Tini+param.gen.burstLength(i)-1,param.T);
    elseif(param.gen.burstLength(i)<=param.gen.burstLengthStdFactor^2)
        Tend = min(Tini+1+poissrnd(param.gen.burstLength(i)),param.T);
    else
        p = param.gen.burstLengthStdFactor^2/param.gen.burstLength(i);
        r = param.gen.burstLengthStdFactor^2/(1-p);
        Tend = min(Tini+1+nbinrnd(r,p),param.T);
    end
    idxBurst = Tini:Tend;
    idx0 = setdiff(1:param.T,idxBurst);
    data.seq(i,idx0) = 0;
    if(Tini+length(param.header)-1>param.T)
        data.seq(i,Tini:param.T) = param.header(1:param.T-Tini+1);
        data.symbols(i,Tini:param.T) = param.constellation(param.header(1:param.T-Tini+1));
    else
        data.seq(i,Tini:Tini+length(param.header)-1) = param.header;
        data.symbols(i,Tini:Tini+length(param.header)-1) = param.constellation(param.header);
    end
    
    % Ensure that first transmitted symbol is in the positive orthant
    if(isempty(param.header))
        aux = randint(1,1,[1 length(idxCpositiveOrthant)]);
        data.symbols(i,Tini) = param.constellation(idxCpositiveOrthant(aux));
        data.seq(i,Tini) = idxCpositiveOrthant(aux);
    end
end
data.symbols(data.seq==0) = 0;

% Generate noise:
noise = sqrt(param.gen.s2n/2)*randn(param.Nr,param.T)+1i*sqrt(param.gen.s2n/2)*randn(param.Nr,param.T);

% Generate the observations:
data.obs = zeros(param.Nr,param.T);
for t=1:param.T
    for i=0:max(param.gen.L_true)-1
        if(t-i<=0)
            data.obs(:,t) = data.obs(:,t) + data.channel(:,:,i+1)*param.gen.symbol0*ones(param.gen.Nt,1);
        else
            data.obs(:,t) = data.obs(:,t) + data.channel(:,:,i+1)*data.symbols(:,t-i);
        end
    end
end

data.obs = data.obs+noise;
% If the constellation is real, remove imaginary part of H and Z
if(sum(abs(imag(param.constellation)))<1e-5*sum(abs(real(param.constellation))))
    data.channel = real(data.channel);
    data.obs = real(data.obs);
end

% Ensure that transmitters become active sequentially
[valnul idxOrder] = sort(Tini_vec);
data.seq = data.seq(idxOrder,:);
data.symbols = data.symbols(idxOrder,:);
data.channel = data.channel(:,idxOrder,:);

