%% Add code folders to current path
addpath(genpath('./code'));

%% Configuration parameters
param.Nr = 20;                       % Number of receivers
param.T  = 100;                      % Length of the sequence
param.constellation = qammod(0:3,4,[],'gray'); % Normalized 4-QAM constellation
param.constellation = param.constellation/sqrt(mean(abs(param.constellation.^2)));
param.flag0 = 1;      % Consider symbol 0 as part of the constellation (if false, transmitters are always active)
param.L = 1;          % Channel memory to be considered during inference
param.Niter = 10000;  % Number of iterations of the sampler
param.saveCycle = inf;   % Save temporary results every saveCycle iterations
param.storeIters = 2000; % How many final iterations to be kept in memory
param.header = [];     % Use a fixed header at the beginning of each data frame?
param.onOffModel = 0;  % A transmitter that switches off may switch on again

%% Generate synthetic data
SNR = 0;               % SNR = -10*log10(noiseVariance)
Ltrue = 1;             % True channel length
Nt = 2;                % Number of transmitters
M = length(param.constellation);
noiseVar = 10^(-SNR/10);
if(log2(M)==1)
    noiseVar = 2*noiseVar;
end
param.gen.s2n = noiseVar;
param.gen.L_true = Ltrue*ones(1,Nt);
param.gen.varH = 1;     % Variance of the channel coefficients
if(log2(M)==1)
    param.gen.varH = 2*param.gen.varH;
end
param.gen.Nt = Nt;
param.gen.burstLength = round(param.T/2)*ones(1,param.gen.Nt);
param.gen.burstLengthStdFactor = inf;
param.gen.symbol0 = 0;
param.gen.sparsityH = 0;

data = generate_data_bursts(param);

%% Configuration parameters for BCJR, PGAS, EP, FFBS and collapsed Gibbs
Nparticles = 3000;     % Number of particles
param.pgas.N_PF = Nparticles;
param.pgas.N_PG = Nparticles;
param.pgas.Niter = 1;
param.pgas.returnNsamples = 1;
param.pgas.maxM = 40;
% Allocate memory for a large number of FSMs (to speed up the C code)
param.pgas.particles = zeros(param.pgas.maxM,max(param.pgas.N_PF,param.pgas.N_PG),param.T,'int16');
param.pgas.flagParallel = 0;  % If 1, it can only be run over multicore machines
param.pgas.blockNtSize = inf; % All transmitters are sampled simultaneously

%% Configuration parameters for BNP and inference method
param.infer.symbolMethod = 'pgas';  % Use PGAS
param.infer.sampleNoiseVar = 0;     % Sample noise variance?
param.infer.sampleChannel = 1;      % Sample channel coefficients?
param.infer.sampleVarH = 1;         % Sample variance of channel coefficients?
param.infer.simulatedTempering = 0;
param.infer.addArtificialNoise = 1; % Use tempering procedure
param.artifNoise.itCycle = 1;
param.artifNoise.stepDB = 12/6000; % This goes from -12dB to 0dB in 6000 iters
param.artifNoise.iniSNR = -12;
param.artifNoise.finalSNR = SNR;

param.bnp.betaSlice1 = 0.5;  % Config parameters for adding new FSMs
param.bnp.betaSlice2 = 5;
param.bnp.maxMnew = 15;
param.bnp.Mini = 1;          % Initial number of FSMs

%% Hyperparameters
hyper.s2h = 1;      % E[s2H(r)]=s2h*exp(-lambda*(r-1))
if(log2(M)==1)
    hyper.s2h = 2*hyper.s2h;
end
hyper.lambda = 0.5; % E[s2H(r)]=s2h*exp(-lambda*(r-1))
hyper.kappa = 1;    % Std[s2H(r)]=kappa*E[s2H(r)]
hyper.alpha = 1;    % Concentration parameter for IBP(alpha)
hyper.gamma1 = 0.1; % Parameter for bm ~ Beta(gamma1,gamma2)
hyper.gamma2 = 2;   % Parameter for bm ~ Beta(gamma1,gamma2)
hyper.tau = 1;      % Parameter for s2y ~ IG(tau,nu)
hyper.nu = 1;       % Parameter for s2y ~ IG(tau,nu)

%% Initialization of the inference
init.H = zeros(param.Nr,param.bnp.Mini,param.L);
for ll=1:param.L
    init.H(:,:,ll) = sqrt(hyper.s2h*exp(-hyper.lambda*(ll-1)))*(randn(param.Nr,param.bnp.Mini,1)+1i*randn(param.Nr,param.bnp.Mini,1));
end
if(~param.infer.sampleNoiseVar)
    init.s2y = noiseVar;      % Initialize s2y to the ground truth
else
    init.s2y = 20*rand(1);    % Initialize s2y to some large value
end
if(param.infer.simulatedTempering)
    init.s2y = param.temper.s2yValues(1);    % Initialize s2y to the largest temperature
end
init.s2H = hyper.s2h*exp(-hyper.lambda*(0:param.L-1));  % Initialize s2H to its mean value
init.am = 0.95*ones(param.bnp.Mini,1);
init.bm = 0.05*ones(param.bnp.Mini,1);
init.Z = zeros(param.bnp.Mini,param.T);
init.seq = zeros(param.bnp.Mini,param.T);
init.nest = zeros(2,2,param.bnp.Mini);
init.nest(1,1,:) = param.T;
init.slice = 0;
init.epAcc = 0;
samples = init;

%% Inference

% Statistics of interest:
LLH = zeros(1,param.Niter);     % log-likelihood at each iteration
M_EST = zeros(1,param.Niter);   % number of estimated FSMs (transmitters) at each iteration
samplesAll = cell(1,param.storeIters); % cell to store the last storeIters samples

% Run inference
for it=1:param.Niter
    %% Algorithm
    if(mod(it,10)==0)
        disp(['Starting iteration ' num2str(it) ' of ' num2str(param.Niter) '...']);
    end
    
    % Step 0a) Add artificial noise
    if(param.infer.addArtificialNoise)
        if(it==1)
            [data.obsWithoutNoise data.obs data.artifNoise samples.s2y] = artifNoise_init(data,samples,hyper,param);
        elseif(mod(it,param.artifNoise.itCycle)==0)
            [data.obs samples.s2y] = artifNoise_decr(data,samples,hyper,param);
        end
    end
    
    % Step 1)
    % -Sample the slice variable
    samples.slice = sample_post_slice(data,samples,hyper,param);
    % -Sample new sticks (and the corresponding new parameters)
    samples = sample_newsticks(data,samples,hyper,param);

    % For PGAS, check that the number of current chains does not exceed maxM
    if(strcmp(param.infer.symbolMethod,'pgas'))
        if(size(samples.seq,1)>param.pgas.maxM)
            param.pgas.maxM = size(samples.seq,1);
            param.pgas.particles = zeros(param.pgas.maxM,max(param.pgas.N_PF,param.pgas.N_PG),param.T,'int16');
        end
    end

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
    % -Sample the channel H
    samples.H = sample_post_H(data,samples,hyper,param);
    % -Sample the noise variance
    samples.s2y = sample_post_s2y(data,samples,hyper,param);
    % -Sample the variance of the channel coefficients
    samples.s2H = sample_post_s2H(data,samples,hyper,param);
        
    %% Store current sample
    if(it>param.Niter-param.storeIters)
        samplesAll{it-param.Niter+param.storeIters} = samples;
    end
    
    %% Evaluation
    % Trace of the estimated number of transmitters
    M_EST(it) = sum(sum(samples.seq~=0,2)>0);
    % Trace of the log-likelihood
    LLH(it) = compute_llh(data,samples,hyper,param);

end

%% Final evaluation of performance (average last storeIters iterations)
Zaux = zeros(size(samples.Z,1),param.T,1+length(param.constellation));
auxSample.s2H = zeros(1,param.L);
for it=1:param.storeIters
    if(size(samplesAll{it}.seq,1)>size(Zaux,1))
        Zaux = cat(1,Zaux,zeros(size(samplesAll{it}.seq,1)-size(Zaux,1),param.T,1+length(param.constellation)));
    end
    for t=1:param.T
        for m=1:size(samplesAll{it}.seq,1)
            Zaux(m,t,samplesAll{it}.seq(m,t)+1) = 1+Zaux(m,t,samplesAll{it}.seq(m,t)+1);
        end
    end
    auxSample.s2H = auxSample.s2H+(samplesAll{it}.s2H/param.storeIters);
end
[valnul auxIdx] = max(Zaux,[],3);
auxConstellation = [0 param.constellation];
auxSample.seq = auxIdx-1;
auxSample.Z = auxConstellation(auxIdx);
auxSample.s2y = samples.s2y;
[valnul auxSample.H] = sample_post_H(data,auxSample,hyper,param);

% Evaluate performance:
% -ADER: activity detection error rate
% -SER_ALL: symbol error rate
% -MMSE: mean square error of the channel coefficients
% -*_indiv: contains the values as above, computed for each individual transmitter
[ADER SER_ALL valnul MMSE vec_ord rot ADER_indiv SER_ALL_indiv valnul2 MMSE_indiv] = compute_error_rates(data,auxSample,hyper,param);
ADER_indiv
SER_ALL_indiv
MMSE_indiv

% Plotting
figure;
plot(M_EST);
xlabel('Iteration');
ylabel('Inferred transmitters')

figure;
plot(LLH);
xlabel('Iteration');
ylabel('Log-likelihood')
