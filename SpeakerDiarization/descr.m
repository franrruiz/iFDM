%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Description of variables and structs %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%% DATA
% 
% data
% -obs: Observations. Size = [1 x T]
% -speakers: Decices power signal. Size = [Nd x T]
% -W: Mixture matrix. Size = [Nd x D]
% -s2y: noise variance

%% SAMPLES
% 
% samples
% -ptrans: Transition probabilities p(x_t!x_t-1). Size = [(Q+1) x Q x Nt]
% -s2y: Noise variance.
% -am: Transition probabilities from 0 to 0. Size = [Nt x 1]
% -bm: Transition probabilities from active to 0. Size = [Nt x 1]
% -Z: State matrix. Size = [Nt x T]
% -W: Mixture matrix. Size = [Nd x D]
% -bX: Lplace distribution parameter
% -nest: Number of jumps from 0/1 to 0/1. Size = [2 x 2 x Nt]
% -slice: Slice variable for the blocked sampling approach.
% -epAcc: Number of times that the MH algorithm accepts the EP proposal (only used if param.infer.symbolMethod='ep').
% 
%% CONFIGURATION PARAMETERS
% 
% param
% -Nd: Number of speakers.
% -D: Dimensionality of the observations, ie., number of microphones.
% -T: Length of the symbol sequence.
% -Niter: Number of iterations of the sampler.
% -ffbs: Struct with the configuration parameters for FFBS algorithm.
%   +Niter: Number of iterations to run.
% -infer: Struct with the options for the inference algorithm.
%   +symbolMethod: Selects the method to infer the symbols. It can be 'ep', 'pgas', 'ffbs' or 'gibbs'.
%   +sampleNoiseVar: Flag to indicate if the noise variace should be sampled.
%   +sampleChannel: Flag to indicate if the channel coefficients (and Nt) should be sampled.
%   +sampleVarH: Flag to indicate if the variance of the channel coefficients should be samples.
% 
%% HYPERPARAMETERS
% 
% hyper
% -bX: variance of X.
% -s2W: Variance for the prior on P.
% -alpha: Concentration parameter of the IBP.
% -gamma1: Hyperparameter for variables, bm~Beta(gamma1,gamma2).
% -gamma2: Hyperparameter for variables, bm~Beta(gamma1,gamma2).
% -nu: Hyperparameter for the noise variance, s2y~InvGamma(nu,tau).
% -tau: Hyperparameter for the noise variance, s2y~InvGamma(nu,tau).
% -nuX: Hyperparameter for the bX variance, bX~InvGamma(nu,tau).
% -tauX: Hyperparameter for the X variance, bX~InvGamma(nu,tau).

% 

