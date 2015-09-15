function data = generate_data_full(param)

% Generate the channel:
data.channel = sqrt(param.gen.varH/2)*randn(param.Nr,param.gen.Nt,max(param.gen.L_true))+1i*sqrt(param.gen.varH/2)*randn(param.Nr,param.gen.Nt,max(param.gen.L_true));
for i=1:param.gen.Nt
    data.channel(:,i,param.gen.L_true(i)+1:max(param.gen.L_true)) = 0;
end
auxSparse = (rand([param.Nr,param.gen.Nt,max(param.gen.L_true)-1])>=param.gen.sparsityH);
data.channel(:,:,2:end) = data.channel(:,:,2:end).*auxSparse;

% Generate the sequence of symbols:
C = length(param.constellation);
data.seq = randint(param.gen.Nt,param.T,[1 C]);
data.symbols = zeros(param.gen.Nt,param.T);
data.symbols(:) = param.constellation(data.seq(:));

% Generate noise:
noise = sqrt(param.gen.s2n/2)*randn(param.Nr,param.T+max(param.gen.L_true))+1i*sqrt(param.gen.s2n/2)*randn(param.Nr,param.T+max(param.gen.L_true));

% Generate the observations:
data.obs = zeros(param.Nr,param.T+max(param.gen.L_true));
for t=1:param.T
    for i=0:max(param.gen.L_true)-1
        if(t-i<=0)
            data.obs(:,t) = data.obs(:,t) + data.channel(:,:,i+1)*param.gen.symbol0*ones(param.gen.Nt,1);
        else
            data.obs(:,t) = data.obs(:,t) + data.channel(:,:,i+1)*data.symbols(:,t-i);
        end
    end
end
% Send L known symbols to recover the initial state:
for t=param.T+1:param.T+max(param.gen.L_true)
    for i=0:param.gen.L_true-1
        if(t-i<=param.T)
            data.obs(:,t) = data.obs(:,t) + data.channel(:,:,i+1)*data.symbols(:,t-i);
        else
            data.obs(:,t) = data.obs(:,t) + data.channel(:,:,i+1)*param.gen.symbol0*ones(param.gen.Nt,1);
        end
    end
end

data.obs = data.obs+noise;
% If the constellation is real, remove imaginary part of H and Y
if(sum(abs(imag(param.constellation)))<1e-5*sum(abs(real(param.constellation))))
    data.channel = real(data.channel);
    data.obs = real(data.obs);
end
