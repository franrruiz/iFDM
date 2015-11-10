%% Add code folders to current path
addpath(genpath('./code'));


%% Configuration parameters
param.Nd = 3;                         % Number of targets
param.D = 25;                         % Number of sensors
param.T  = 300;                       % Length of the sequence
param.d0 = 1;                         % Reference distance
param.Ts  = 0.04;                     % Sampling time
param.pathL= 2;                       % Path loss exponent
param.L = 1;
param.flag0 = 1;                      % Consider symbol 0 as part of the constellation (if false, transmitters are always active)
param.Niter = 1000;                   % Number of iterations of the sampler
param.saveCycle = 200;
param.storeIters = 2000;
param.constellation=1;

%% Generate data
data.s2y=2;
data.s2u=1;
data.s2vIni=0.01;
data.Ptx= 50; % Transmitted power in dB
data.W=800;
data.Gx=[1 0 param.Ts 0; 0 1 0 param.Ts; 0 0 1 0; 0 0 0 1];
data.Gu= [param.Ts^2/2 0; 0 param.Ts^2/2; param.Ts 0; 0 param.Ts];
pos = 0:round(data.W/(sqrt(param.D)-1)):data.W;
[a b] = meshgrid(pos);
data.sensors = [a(:) b(:)];

data.state=zeros(param.D,4,param.T);
for nt=1:param.Nd
    Tini = 1;%randi([1 round(T/2)],1,1);
    Tend = param.T;%min(Tini+round(T/2)-1,T);
    data.state(nt,:,Tini)= [data.W*rand(2,1); sqrt(data.s2vIni)*randn(2,1)];
    flagTxDone = 0;
    while(~flagTxDone)
        for t=Tini+1:Tend
            data.state(nt,:,t) = data.Gx*squeeze(data.state(nt,:,t-1))'+sqrt(data.s2u)*data.Gu*randn(2,1);
        end
        if(~(sum(data.state(nt,1,:)<-data.W/5)>0 || sum(data.state(nt,1,:)>data.W+data.W/5)>0 || sum(data.state(nt,2,:)<-data.W/5)>0 || sum(data.state(nt,2,:)>data.W+data.W/5)>0))
            flagTxDone = 1;
        end
    end
end

Ptx=zeros(param.D,param.T);
for t=1:param.T
    d=zeros(param.D,1);
    for nt=1:param.Nd
        if(data.state(nt,1,t)~=0)
        	d= d+1./(((data.sensors(:,1)-data.state (nt,1,t)).^2 +(data.sensors(:,2)-data.state (nt,2,t)).^2).^(param.pathL/2));
        end
    end
    Ptx(d>0,t)= 10^(data.Ptx/10)*param.d0^param.pathL*d(d>0);
end

data.obs = Ptx+sqrt(data.s2y)*randn(param.D,param.T);



%% Configuration parameters for BCJR, PGAS, EP, FFBS and collapsed Gibbs
param.bcjr.p1 = 0.95;
param.bcjr.p2 = 0.05;
param.pgas.N_PF = 3000;
param.pgas.N_PG = 10000;
param.pgas.Niter = 1;
param.pgas.returnNsamples = 1;
param.pgas.maxM = 40;
%param.pgas.particles = zeros(param.pgas.maxM,max(param.pgas.N_PF,param.pgas.N_PG),param.T,'int16');
param.ep.eps = 5e-7;
param.ep.beta = 0.2;
param.ep.Niter = 15;
param.colGibbs.Niter = 1;
param.ffbs.Niter = 1;

%% Configuration parameters for BNP and inference method
param.infer.symbolMethod = 'pgas';
param.infer.sampleNoiseVar = 1;
param.bnp.betaSlice1 = 0.5;
param.bnp.betaSlice2 = 5;
param.bnp.maxMnew = 15;
param.bnp.Mini = 1;


%% Hyperparameters
hyper.alpha = 1;    % Concentration parameter for Z ~ IBP(alpha)
hyper.gamma1 = 0.1; % Parameter for bm ~ Beta(gamma1,gamma2)
hyper.gamma2 = 2;   % Parameter for bm ~ Beta(gamma1,gamma2)
hyper.tau = 1;      % Parameter for s2y ~ IG(tau,nu)
hyper.nu = 1;       % Parameter for s2y ~ IG(tau,nu)

%% Initialization
if param.infer.sampleNoiseVar
    init.s2y = 100*rand;    % INITIALIZE s2y at random
else
    init.s2y = data.s2y;    % INITIALIZE s2y TO THE GROUND TRUTH
end
init.am = 0.95*ones(param.bnp.Mini,1);
init.bm = 0.05*ones(param.bnp.Mini,1);
init.Z = zeros(param.bnp.Mini,4,param.T);
init.nest = zeros(2,2,param.bnp.Mini);
init.nest(1,1,:) = param.T;
init.slice = 0;
init.epAcc = 0;
samples = init;
samplesAll = cell(1,param.storeIters);

%% Inference
for it=1:param.Niter
    %% Algorithm
    
    % Step 1)
    % -Sample the slice variable
    samples.slice = sample_post_slice(data,samples,hyper,param);
    % -Sample new sticks (and the corresponding new parameters)
    samples = sample_newsticks(data,samples,hyper,param);
    
    % Step 2)
    % -Sample the symbols Z
    [samples.Z samples.seq samples.nest out] = sample_post_Z(data,samples,hyper,param);
    % -Compute some statistics of interest
    if(strcmp(param.infer.symbolMethod,'pgas'))
        
    elseif(strcmp(param.infer.symbolMethod,'ep'))
        samples.epAcc = samples.epAcc+out;
    end
    
    % Step 3)
    % -Remove unused chains
    samples = sample_remove_unused(data,samples,hyper,param);
    
    % Step 4)
    % -Sample the transition probabilities (semi-ordered construction)
    [samples.am samples.bm]= sample_post_transitionProb(data,samples,hyper,param);
    
    % Step 5)
    % -Sample the noise variance
    samples.s2y = sample_post_s2y(data,samples,hyper,param);

    %% Evaluation
%     % Trace of the estimated number of transmitters
     M_EST(it) = sum(sum(samples.Z(:,1,:)~=0,3)>0);
%     % Trace of the log-likelihood
     LLH(it) = compute_llh(data,samples,hyper,param);
    
end

