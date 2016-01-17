function power_disaggregation(input_file, output_file, config_file)
%
% This function applies the iFDM model [1] for power disaggregation.
% 
% Usage:
%     power_disaggregation(input_file, output_file [, config_file])
% 
%  +input_file: Path to the input file. It should contain a matrix in which
%   the first row corresponds to the observations (aggregated signal),
%   while the rest of the rows correspond to the individual consuption of
%   the devices. The number of columns of the matrix corresponds to the
%   number of time steps.
% 
%  +output_file: Path to the output file. Two variables are saved in the
%   output file, detailed below. Warning: These results correspond to the
%   last iteration of the sampler, they have not been averaged over several
%   iterations.
%   -acc: The resulting accuracy.
%   -output: The inferred disaggregated signals, sorted to match the input.
% 
%  +config_file: Path to a JSON configuration file (optional). If specified,
%   it may contain the following variables:
%   -Q: Number of active states of the devices (default=4).
%   -Niter: Number of iterations of the algorithm (default=1000).
%   -verboseCycle: Print progress on stdout every verboseCycle iterations
%    (default=10).
%   -pgas: Struct to configure the PGAS kernel. It may contain the 
%          following fields:
%          -N_PG: Number of particles for the PGAS kernel (default=3000).
%          -Niter: Number of iterations per call to the PGAS function
%                  (default=1).
%          -maxM: Maximum number of parallel chains allowed before resizing
%                 matrices. A high value is recommended for memory
%                 efficiency reasons (default=40).
%   -bnp: Struct to configure the BNP-specific part of the algorithm. It may
%         contain the following fields:
%         -betaSlice1: See below (default=0.5).
%         -betaSlice2: Parameters of the Beta distribution to sample the
%                      slice mIBP variable (default=5).
%         -maxMnew: Maximum number of new chains allowed per iteration
%                   (default=15).
%         -Mini: Initial number of parallel chains (default=1).
%   -hyper: Struct to set the hyperparameters of the model. It may contain
%           the following fields:
%           -s2y: Noise variance of the observations (default=0.5).
%           -s2P: Variance of the power values (default=10).
%           -muP: Mean of the power values (default=15).
%           -gamma: Hyperparameter over the Dirichlet distribution for the
%                   transition probabilities of the active states (default=1).
%           -alpha: Concentration parameter of the mIBP (default=1).
%           -gamma1: See below (default=0.1).
%           -gamma2: Hyperparameters of the Beta distribution over the
%                    transition probabilities from active to inactive
%                    (default=2).
% 
% References:
% [1] I. Valera, F. J. R. Ruiz, L. Svensson, F. Perez-Cruz. "Infinite
%     Factorial Dynamical Model". NIPS 2015.
% 

%% Load data 
input=importdata(input_file);

data.obs = input(1,:);
data.devices = input(2:end,:);
clear input;
%% Add code folders to current path
addpath(genpath('./code'));
addpath(genpath('./jsonlab-1.2'));

%% Configuration parameters
if(nargin>=3 && ~isempty(config_file))
    fname=sprintf(config_file);
    json=savejson([],loadjson(fname));
    param=loadjson(json);
end
param.Nd = size(data.devices,1);     % Number of devices
param.D = 1;                         % Dimensionality of the observations
param.T  = size(data.devices,2);     % Length of the sequence
param.flag0 = 1;

% Global parameters for PGAS
if ~isfield(param,'Q')
    param.Q = 4;
end
if ~isfield(param,'Niter')
    param.Niter = 1000;
end
if ~isfield(param,'verboseCycle')
    param.verboseCycle = 10;
end

% Configuration parameters for PGAS
if ~isfield(param,'pgas')
    param.pgas=[];
end
if ~isfield(param.pgas,'N_PG')
    param.pgas.N_PG = 3000;
end
if ~isfield(param.pgas,'Niter')
    param.pgas.Niter = 1;
end
if ~isfield(param.pgas,'maxM')
    param.pgas.maxM = 40;
end
param.pgas.N_PF = 3000;
param.pgas.returnNsamples = 1;
param.pgas.particles = zeros(param.pgas.maxM,max(param.pgas.N_PF,param.pgas.N_PG),param.T,'int16');

% Configuration parameters for BNP and inference method
param.infer.symbolMethod = 'pgas'; %Set param.infer.symbolMethod = 'ffbs'; for running the forward filtering-backward sampling algorithm
param.infer.sampleNoiseVar = 0;    % We do not sample the noise variance
param.infer.sampleP = 1;
param.infer.sampleVarP = 0;
if ~isfield(param,'bnp')
    param.bnp=[];
end
if ~isfield(param.bnp,'betaSlice1')
    param.bnp.betaSlice1 = 0.5;
end
if ~isfield(param.bnp,'betaSlice2')
    param.bnp.betaSlice2 = 5;
end
if ~isfield(param.bnp,'maxMnew')
    param.pgas.maxMnew = 15;
end
if ~isfield(param.bnp,'Mini')
    param.pgas.Mini = 1;
end

% Hyperparameters
if ~isfield(param,'hyper')
    param.hyper = [];
end
if ~isfield(param.hyper,'s2y')
    param.hyper.s2y = 0.5; % Noise variance
end
if ~isfield(param.hyper,'s2P')
    param.hyper.s2P = 10; % Prior Pqm, power of state q in chain m is gaussian distributed
end
if ~isfield(param.hyper,'muP')
    param.hyper.muP = 15; % Prior Pqm, power of state q in chain m is gaussian distributed
end
if ~isfield(param.hyper,'gamma')
    param.hyper.gamma = 1; % prior over the transition probabilities from x_t-1 to x_t forllows a dirichlet with Q components and parameter gamma
end
if ~isfield(param.hyper,'alpha')
    param.hyper.alpha = 1; % Concentration parameter for Z ~ IBP(alpha)
end
if ~isfield(param.hyper,'gamma1')
    param.hyper.gamma1 = 0.1; % Parameter for bm ~ Beta(gamma1,gamma2)
end
if ~isfield(param.hyper,'gamma2')
    param.hyper.gamma2 = 2; % Parameter for bm ~ Beta(gamma1,gamma2)
end

%% Initialization
init.P = param.hyper.muP+sqrt(param.hyper.s2P)*randn(param.Q,param.bnp.Mini);
for mm=1:param.bnp.Mini
    init.ptrans(:,:,mm) = dirichletrnd(param.hyper.gamma*ones(1,param.Q), param.Q+1);
end
init.s2y = param.hyper.s2y;      % INITIALIZE s2y TO THE GROUND TRUTH
init.am = 0.95*ones(param.bnp.Mini,1);
init.bm = 0.05*ones(param.bnp.Mini,1);
init.Z = zeros(param.bnp.Mini,param.T);
init.nest = zeros(2,2,param.bnp.Mini);
init.nest(1,1,:) = param.T;
init.slice = 0;
samples = init;
LLH = zeros(1,param.Niter+1);
M_EST = zeros(1,param.Niter+1);    


%% Inference
fprintf(1,['Initializing inference...' char(13)]);
for it=1:param.Niter
    %% Print progress
    if(mod(it,param.verboseCycle)==0)
        fprintf(1,['Iteration ' num2str(it) ' out of ' num2str(param.Niter) ' completed' char(13)]);
    end
    
    %% Algorithm
    % Step 1)
    % -Sample the slice variable
    samples.slice = sample_post_slice(data,samples,param.hyper,param);
    % -Sample new sticks (and the corresponding new parameters)
    samples = sample_newsticks(data,samples,param.hyper,param);
    
    % For PGAS, check that the number of current chains does not exceed maxM
    if(strcmp(param.infer.symbolMethod,'pgas'))
        if(size(samples.Z,1)>param.pgas.maxM)
            param.pgas.maxM = size(samples.Z,1);
            param.pgas.particles = zeros(param.pgas.maxM,max(param.pgas.N_PF,param.pgas.N_PG),param.T,'int16');
        end
    end
    
    % Step 2)
    % -Sample the symbols Z
    [samples.Z samples.seq samples.nest out] = sample_post_Z(data,samples,param.hyper,param);
    
    % Step 3)
    % -Remove unused chains
    samples = sample_remove_unused(data,samples,param.hyper,param);
    
    % Step 4)
    % -Sample the transition probabilities (semi-ordered construction)
    [samples.am samples.bm]= sample_post_transitionProb(data,samples,param.hyper,param);
    
    % Step 5)
    % -Sample the mean power associated to each device
    samples.P = sample_post_P(data,samples,param.hyper,param);
    % -Sample tptrans
    samples.ptrans = sample_post_ptrans(data,samples,param.hyper,param);
    
    %% Evaluation
    % Trace of the estimated number of transmitters
    M_EST(it) = sum(sum(samples.Z~=0,2)>0);
    % Trace of the log-likelihood
    LLH(it) = compute_llh(data,samples,param.hyper,param);
end

%% Obtain individual inferred chains in terms of power, not states
samples.disaggregated = samples.Z;
samples.Paux = [zeros(1,size(samples.Z,1)); samples.P];
for m=1:size(samples.Z,1)
    samples.disaggregated(m,:) = samples.Paux(samples.Z(m,:)+1,m);
end
%% Compute Accuracy
[acc output] = computeAccuracy(samples.disaggregated,data.devices);

save(output_file,'acc','output');
