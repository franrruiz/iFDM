%% Add code folders to current path
addpath(genpath('./code'));
addpath(genpath('./REDD'));
%addpath(genpath('./AMPds'));

%% Configuration parameters
H = 1;                              % House from REDD data
Tini = 1;                           % Intinial time of the sequence
Tend = 2880;                        % Final time of the sequence
param.Nd = 4;                       % Number of devices
param.D = 1;                        % Dimensionality of the observations
param.T  = Tend-Tini+1;             % Length of the sequence
param.Q = 4;                        % Num,ber of states
param.flag0 = 1;                    % Consider symbol 0 as part of the constellation (if false, transmitters are always active)
param.Niter = 10000;                % Number of iterations of the sampler
%param.constellation=1:param.Q;

idxDevOrder = [1,2,3,4,7,14,16,10,5,6,8,9,11,13,15,18,19,12,17];

%% Load REDD data 
load('REDD/data/idxDevNew.mat');
load(['REDD/data/dataH' num2str(H) '_2d_t30.mat']);
dev = [];
X = zeros(1,size(devices,2));
for i=1:param.Nd
    X = X+sum(devices(idxDevNew{H,idxDevOrder(i)},:),1);
    dev = [dev;sum(devices(idxDevNew{H,idxDevOrder(i)},:),1)];
end
data.obs = X(:,Tini:Tend)/100;
data.devices = dev(:,Tini:Tend)/100;

% %% Load AMPds data
% %% NOTE: AMPds is *not* available because the usage of this dataset requires
% %%       approval by the author. See http://ampds.org
% idxDevOrder = [3     4     7    10    13    15    17    19];
% Nd = length(idxDevOrder);
% load('AMPds/data/AMPds_data.mat','devices');
% devices = devices(idxDevOrder(1:Nd),Tini:Tend)/100;
% data.obs = sum(devices,1);
% data.devices = devices;


%% Configuration parameters for BCJR, PGAS, EP, FFBS and collapsed Gibbs
param.bcjr.p1 = 0.95;
param.bcjr.p2 = 0.05;
param.pgas.N_PF = 3000;
param.pgas.N_PG = 3000;
param.pgas.Niter = 1;
param.pgas.returnNsamples = 1;
param.pgas.maxM = 40;
param.pgas.particles = zeros(param.pgas.maxM,max(param.pgas.N_PF,param.pgas.N_PG),param.T,'int16');
param.ffbs.Niter = 1;

%% Configuration parameters for BNP and inference method
param.infer.symbolMethod = 'pgas'; %Set param.infer.symbolMethod = 'ffbs'; for running the forward filtering-backward sampling algorithm
param.infer.sampleNoiseVar = 0;
param.infer.sampleP = 1;
param.infer.sampleVarP = 0;
param.bnp.betaSlice1 = 0.5;
param.bnp.betaSlice2 = 5;
param.bnp.maxMnew = 15;
param.bnp.Mini = 1;


%% Hyperparameters
hyper.s2P = 10;     % Prior Pqm, power of state q in chain m is gaussian distributed
hyper.muP = 15;     % Prior Pqm, power of state q in chain m is gaussian distributed
hyper.gamma = 1;    % prior over the transition probabilities from x_t-1 to x_t forllows a dirichlet with Q components and parameter gamma
hyper.alpha = 1;    % Concentration parameter for Z ~ IBP(alpha)
hyper.gamma1 = 0.1; % Parameter for bm ~ Beta(gamma1,gamma2)
hyper.gamma2 = 2;   % Parameter for bm ~ Beta(gamma1,gamma2)
hyper.tau = 1;      % Parameter for s2y ~ IG(tau,nu)
hyper.nu = 1;       % Parameter for s2y ~ IG(tau,nu)

%% Initialization
init.P = hyper.muP+sqrt(hyper.s2P)*randn(param.Q,param.bnp.Mini);
for mm=1:param.bnp.Mini
    init.ptrans(:,:,mm) = dirichletrnd(hyper.gamma*ones(1,param.Q), param.Q+1);
end
init.s2y = 0.5;      % INITIALIZE s2y TO THE GROUND TRUTH
init.am = 0.95*ones(param.bnp.Mini,1);
init.bm = 0.05*ones(param.bnp.Mini,1);
init.Z = zeros(param.bnp.Mini,param.T);
init.nest = zeros(2,2,param.bnp.Mini);
init.nest(1,1,:) = param.T;
init.slice = 0;
init.epAcc = 0;
samples = init;
LLH = zeros(1,param.Niter+1);
M_EST = zeros(1,param.Niter+1);    


%% Inference
for it=1:param.Niter
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
    [samples.Z samples.seq samples.nest out] = sample_post_Z(data,samples,hyper,param);
    
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
    
     %% Evaluation
     % Trace of the estimated number of transmitters
     M_EST(it) = sum(sum(samples.Z~=0,2)>0);
     % Trace of the log-likelihood
     LLH(it) = compute_llh(data,samples,hyper,param);
end

%% Compute Accuracy
[ACC cad_ord]= computeAccuracy(samples.Z,data.devices);

