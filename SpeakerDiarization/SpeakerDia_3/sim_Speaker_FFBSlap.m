function sim_Speaker_FFBSlap(Tsubsample,Nd,D,Niter,LastIt)

addpath(genpath('sampleFunc/'));
addpath(genpath('auxFunc/'));

randn('seed',round(sum(1e5*clock)));
rand('seed',round(sum(1e5*clock)));

flagRecovered=0;
if LastIt>0
    flagRecovered=1;
end

%% Configuration parameters
param.D = D;
param.Nd = Nd;                        % Number of devices
param.L = 1;
param.flag0 = 1;    % Consider symbol 0 as part of the constellation (if false, transmitters are always active)
param.Niter = Niter;  % Number of iterations of the sampler
param.saveCycle = 200;
param.storeIters = 2000;
param.Q=1;
%param.constellation=1:Q;

load('PCCdata16kHz_isolated/data/dataF4.mat','speakers','obs','W','s2y');
data.speakers = (speakers(1:Tsubsample:end,1:param.Nd));
[param.T aux1 aux2]=size(data.speakers);
data.W=W;
data.s2y=s2y;
data.obs = obs(:,1:Tsubsample:end);

BASEDIR1=['PCCdata16kHz_isolated/resultsFFBS_lap/S' num2str(param.Nd) '_T' num2str(param.T) '_Tsub' num2str(Tsubsample)];
if(~isdir(BASEDIR1))
    mkdir(BASEDIR1);
end
%% Configuration parameters for BCJR, PGAS, EP, FFBS and collapsed Gibbs
param.bcjr.p1 = 0.95;
param.bcjr.p2 = 0.05;
param.pgas.N_PF = 3000;
param.pgas.N_PG = 100;
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
param.infer.symbolMethod = 'ffbs_laplace';
param.infer.sampleNoiseVar = 1;
param.infer.sampleWVar = 0;
param.infer.sampleVarX = 0;
param.bnp.betaSlice1 = 0.5;
param.bnp.betaSlice2 = 5;
param.bnp.maxMnew = 15;
param.bnp.Mini = 1;


%% Hyperparameters
hyper.s2W = 1;      % Prior Pqm, power of state q in chain m is gaussian distributed
hyper.bX = 2;      % Prior X, which is Gaussian(0,bX)
hyper.alpha = 1;    % Concentration parameter for Z ~ IBP(alpha)
hyper.gamma1 = 0.1; % Parameter for bm ~ Beta(gamma1,gamma2)
hyper.gamma2 = 2;   % Parameter for bm ~ Beta(gamma1,gamma2)
hyper.tau = 2;      % Parameter for s2y ~ IG(tau,nu)
hyper.nu = 1;       % Parameter for s2y ~ IG(tau,nu)
hyper.tauW = 2;      % Parameter for s2W ~ IG(tauW,nuW)
hyper.nuW = 1;       % Parameter for s2W ~ IG(tauW,nuW)
hyper.taubX = 2;      % Parameter for bX ~ IG(taubX,nubX)
hyper.nubX = 1;       % Parameter for bX ~ IG(taubX,nubX)


%% Initialization
if(~flagRecovered)
    if param.infer.sampleWVar
        init.s2W=2*rand;
    else
        init.s2W=hyper.s2W;
    end
    if param.infer.sampleVarX
        init.bX=5*rand;
    else
        init.bX=hyper.bX;
    end
    
    init.W = sqrt(init.s2W)*rand(param.bnp.Mini, param.D);   
    if param.infer.sampleNoiseVar
        init.s2y = 20*rand;      % INITIALIZE s2y TO THE GROUND TRUTH
    else
        init.s2y = data.s2y;
    end
    init.am = 0.95*ones(param.bnp.Mini,1);
    init.bm = 0.05*ones(param.bnp.Mini,1);
    init.Z = zeros(param.bnp.Mini,param.T);
    init.seq = zeros(param.bnp.Mini,param.T);
    init.nest = zeros(2,2,param.bnp.Mini);
    init.nest(1,1,:) = param.T;
    init.slice = 0;
    init.epAcc = 0;
    samples = init;
    samplesAll = cell(1,param.storeIters);
else
    load([BASEDIR1 '/it' num2str(LastIt) '.mat'],'data','init','samples','samplesAll','LLH', 'M_EST');
    
end

%% Inference
for it=LastIt+1:param.Niter
    %% Algorithm
    
    % Step 1)
    % -Sample the slice variable
    samples.slice = sample_post_slice(data,samples,hyper,param);
    % -Sample new sticks (and the corresponding new parameters)
    samples = sample_newsticks(data,samples,hyper,param);
    
    % Step 2)
    % -Sample the symbols Z
    [samples.Z samples.seq samples.nest out] = sample_post_Z_NEW(data,samples,hyper,param);
    
    % Step 3)
    % -Remove unused chains
    samples = sample_remove_unused(data,samples,hyper,param);
    
    % Step 4)
    % -Sample the transition probabilities (semi-ordered construction)
    [samples.am samples.bm]= sample_post_transitionProb(data,samples,hyper,param);
    
    % Step 5)
    % -Sample the mean power associated to each device
    samples.W = sample_post_W(data,samples,hyper,param);
   
    % -Sample the noise variance
    samples.s2y = sample_post_s2y(data,samples,hyper,param);
    % -Sample the variance of the channel coefficients
    samples.s2W = sample_post_s2W(data,samples,hyper,param);
    % -Sample the variance of the channel coefficients
    samples.bX = sample_post_s2X(data,samples,hyper,param);
    
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
end
save([BASEDIR1 '/Final.mat'],'data','init','samples','samplesAll', 'LLH', 'M_EST');


