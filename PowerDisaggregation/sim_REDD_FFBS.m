function sim_REDD_FFBS(H,noiseVar,Tini,Tend,Nd,Q,Niter,LastIt)

addpath(genpath('sampleFunc/'));
addpath(genpath('auxFunc/'));

randn('seed',round(sum(1e5*clock)));
rand('seed',round(sum(1e5*clock)));

flagRecovered=0;
if LastIt>0
    flagRecovered=1;
end

%% Configuration parameters
param.Nd = Nd;                        % Number of devices
param.D = 1;                        %Dimensionality of the observations
param.T  = Tend-Tini+1;                         % Length of the sequence
param.Q = Q;
param.L = 1;
param.flag0 = 1;    % Consider symbol 0 as part of the constellation (if false, transmitters are always active)
param.Niter = Niter;  % Number of iterations of the sampler
param.saveCycle = 200;
param.storeIters = 2000;
%param.onOffModel = onOffModel;
param.constellation=1:Q;

idxDevOrder = [1,2,3,4,7,14,16,10,5,6,8,9,11,13,15,18,19,12,17];

%% Load data
BASEDIR1=['REDD/resultsFFBS/House' num2str(H) '_M' num2str(param.Nd) '_T' num2str(Tini) '_' num2str(Tend)];
if(~isdir(BASEDIR1))
    mkdir(BASEDIR1);
end
load('REDD/data/idxDevNew.mat');
load(['REDD/data/datosH' num2str(H) '_2d_t30.mat']);
dev = [];
X = zeros(1,size(devices,2));
for i=1:Nd
    X = X+sum(devices(idxDevNew{H,idxDevOrder(i)},:),1);
    dev = [dev;sum(devices(idxDevNew{H,idxDevOrder(i)},:),1)];
end
data.obs = X(:,Tini:Tend)/100;
data.devices = dev(:,Tini:Tend)/100;

%% Configuration parameters for BCJR, PGAS, EP, FFBS and collapsed Gibbs
param.bcjr.p1 = 0.95;
param.bcjr.p2 = 0.05;
param.pgas.N_PF = 3000;
param.pgas.N_PG = 3000;
param.pgas.Niter = 1;
param.pgas.returnNsamples = 1;
param.pgas.maxM = 40;
param.pgas.particles = zeros(param.pgas.maxM,max(param.pgas.N_PF,param.pgas.N_PG),param.T,'int16');
param.ep.eps = 5e-7;
param.ep.beta = 0.2;
param.ep.Niter = 15;
param.colGibbs.Niter = 1;
param.ffbs.Niter = 1;

%% Configuration parameters for BNP and inference method
param.infer.symbolMethod = 'ffbs';
param.infer.sampleNoiseVar = 0;
param.infer.sampleP = 1;
param.infer.sampleVarP = 0;
param.bnp.betaSlice1 = 0.5;
param.bnp.betaSlice2 = 5;
param.bnp.maxMnew = 15;
param.bnp.Mini = 1;


%% Hyperparameters
hyper.s2P = 10;      % Prior Pqm, power of state q in chain m is gaussian distributed
hyper.muP = 15;
hyper.gamma = 1;    % prior over the transition probabilities from x_t-1 to x_t forllows a dirichlet with Q components and parameter gamma
%hyper.kappa = 1;    % Std[s2H(r)]=kappa*E[s2H(r)]
hyper.alpha = 1;    % Concentration parameter for Z ~ IBP(alpha)
hyper.gamma1 = 0.1; % Parameter for bm ~ Beta(gamma1,gamma2)
hyper.gamma2 = 2;   % Parameter for bm ~ Beta(gamma1,gamma2)
hyper.tau = 1;      % Parameter for s2y ~ IG(tau,nu)
hyper.nu = 1;       % Parameter for s2y ~ IG(tau,nu)

%% Initialization
if(~flagRecovered)
    init.P = hyper.muP+sqrt(hyper.s2P)*randn(param.Q,param.bnp.Mini);
    for mm=1:param.bnp.Mini
        init.ptrans(:,:,mm) = dirichletrnd(hyper.gamma*ones(1,param.Q), param.Q+1);
    end
    init.s2y = noiseVar;      % INITIALIZE s2y TO THE GROUND TRUTH
    init.am = 0.95*ones(param.bnp.Mini,1);
    init.bm = 0.05*ones(param.bnp.Mini,1);
    init.Z = zeros(param.bnp.Mini,param.T);
    init.nest = zeros(2,2,param.bnp.Mini);
    init.nest(1,1,:) = param.T;
    init.slice = 0;
    init.epAcc = 0;
    samples = init;
    samplesAll = cell(1,param.storeIters);
    LLH = zeros(1,param.Niter+1);
    M_EST = zeros(1,param.Niter+1);
else
    load([BASEDIR1 '/it' num2str(LastIt) '.mat'],'data','init','samples','samplesAll');
    
end

%% Inference
for it=LastIt+1:param.Niter
    %tic
    %% Algorithm
    
    % Step 1)
    % -Sample the slice variable
    samples.slice = sample_post_slice(data,samples,hyper,param);
    % -Sample new sticks (and the corresponding new parameters)
    samples = sample_newsticks(data,samples,hyper,param);
    
    % For PGAS, check that the number of current chains does not exceed maxM
    if(strcmp(param.infer.symbolMethod,'pgas'))
        if(size(samples.Z,1)>param.pgas.maxM)
            param.pgas.maxM = size(samples.Z,1);
            param.pgas.particles = zeros(param.pgas.maxM,max(param.pgas.N_PF,param.pgas.N_PG),param.T,'int16');
        end
    end
    
    % Step 2)
    % -Sample the symbols Z
    %tic;
    [samples.Z samples.seq samples.nest out] = sample_post_Z(data,samples,hyper,param);
    %toc
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
    % -Sample the mean power associated to each device
    samples.P = sample_post_P(data,samples,hyper,param);
    % -Sample tptrans
    samples.ptrans = sample_post_ptrans(data,samples,hyper,param);
    
    %% Store current sample
    if(it>param.Niter-param.storeIters)
        samplesAll{it-param.Niter+param.storeIters} = samples;
    end
    
%     %% Evaluation
%     % Trace of the estimated number of transmitters
     M_EST(it) = sum(sum(samples.Z~=0,2)>0);
%     % Trace of the log-likelihood
     LLH(it) = compute_llh(data,samples,hyper,param);
    
    %% Save temporary result file
    if(mod(it,param.saveCycle)==0)
        save([BASEDIR1 '/it' num2str(it) '.mat'],'data','init','samples','samplesAll', 'LLH', 'M_EST');
        % If successfully saved, delete previous temporary file
%         if(exist([saveFile '/it' num2str(it-param.saveCycle) '.mat'],'file'))
%             delete([saveFile '/it' num2str(it-param.saveCycle) '.mat']);
%         end
    end
    %toc
end
save([BASEDIR1 '/Final.mat'],'data','init','samples','samplesAll', 'LLH', 'M_EST');

