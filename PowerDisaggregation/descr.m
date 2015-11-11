%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Description of variables and structs %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%% DATA
% 
% data
% -obs: Observations. Size = [1 x T]
% -devices: Decices power signal. Size = [Nd x T]
% 

%% SAMPLES
% 
% samples
% -P: Power values. Size = [Q x Nd]
% -ptrans: Transition probabilities p(x_t!x_t-1). Size = [(Q+1) x Q x Nt]
% -s2y: Noise variance.
% -am: Transition probabilities from 0 to 0. Size = [Nt x 1]
% -bm: Transition probabilities from active to 0. Size = [Nt x 1]
% -s2H: Variance of the channel coefficients. Size = [1 x L]
% -Z: State matrix. Size = [Nt x T]
% -nest: Number of jumps from 0/1 to 0/1. Size = [2 x 2 x Nt]
% -slice: Slice variable for the blocked sampling approach.
% -epAcc: Number of times that the MH algorithm accepts the EP proposal (only used if param.infer.symbolMethod='ep').
% 
%% CONFIGURATION PARAMETERS
% 
% param
% -flag0: Flag to determine if symbol 0 should be considered.
% -Nd: Number of devices.
% -D: Dimensionality of the observations.
% -T: Length of the symbol sequence.
% -Niter: Number of iterations of the sampler.
% -saveCycle: Save state of the sampler every saveCycle iterations.
% -storeIters: Number of iterations to be stored as local variables.
% -bcjr: Struct with configuration parameters for BCJR algorithm.
%   +p1: Probability of transitioning from 0 to 0.
%   +p2: Probability of transitioning from active to 0.
% -pgas: Struct with configuration parameters for PGAS algorithm.
%   +N_PF: Number of particles for the PF.
%   +N_PG: Number of particles for the PGAS.
%   +Niter: Number of iterations of the PGAS sampler.
%   +returnNsamples: Number of samples that the PGAS sampler should return.
%   +maxM: Maximum number of chains allowed. If exceeded, nothing happens but it will be slower
%   +particles: Matrix containing the particles. It has to be initialized due to efficient memory use. Size = [maxM x max(N_PF,N_PG) x T]
% -bnp: Struct with the configuration parameters for BNP inference.
%   +betaSlice1: Beta parameter to sample from the slice variable: s=betarnd(betaSlice1,betaSlice2)*c_min.
%   +betaSlice2: Beta parameter to sample from the slice variable: s=betarnd(betaSlice1,betaSlice2)*c_min.
%   +Mini: Initial number of parallel chains
%   +maxMnew: Maximum number of new chains allowed per iteration.
% -colGibbs: Struct with the configuration parameters for collapsed Gibbs algorithm.
%   +Niter: Number of iterations to run.
% -ffbs: Struct with the configuration parameters for FFBS algorithm.
%   +Niter: Number of iterations to run.
% -ep: Struct with the configuration parameters for EP algorithm.
%   +eps: Epsilon coefficient (minimum variance allowed).
%   +beta: Damping coefficient.
%   +Niter: Number of iterations for the EP algorithm.
% -infer: Struct with the options for the inference algorithm.
%   +symbolMethod: Selects the method to infer the symbols. It can be 'ep', 'pgas', 'ffbs' or 'gibbs'.
%   +sampleNoiseVar: Flag to indicate if the noise variace should be sampled.
%   +sampleChannel: Flag to indicate if the channel coefficients (and Nt) should be sampled.
%   +sampleVarH: Flag to indicate if the variance of the channel coefficients should be samples.
% 
%% HYPERPARAMETERS
% 
% hyper
% -muP: Mean for the prior on P.
% -s2P: Variance for the prior on P.
% -lambda: Define the mean value of the variance of the channel coefficients: E[s2H(r)]=s2h*exp(-lambda*(r-1)).
% -hyper.gamma: prior over the transition probabilities from x_t-1 to x_t follows a dirichlet with Q components and parameter gamma
% -alpha: Concentration parameter of the IBP.
% -gamma1: Hyperparameter for variables, bm~Beta(gamma1,gamma2).
% -gamma2: Hyperparameter for variables, bm~Beta(gamma1,gamma2).
% -nu: Hyperparameter for the noise variance, s2y~InvGamma(nu,tau).
% -tau: Hyperparameter for the noise variance, s2y~InvGamma(nu,tau).
% 

