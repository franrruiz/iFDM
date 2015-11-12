%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Description of variables and structs %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%% DATA
% 
% data
% -obs: Observations (including artificial noise). Size = [Nr x T]
% -obsWithoutNoise: Observations (without artificial noise). Size = [Nr x T]
% -artifNoise: Artificial noise. Size = [Nr x T]
% -channel: Channel coefficients. Size = [Nr x Nt x L]
% -symbols: Transmitted symbols. Size = [Nt x T]
% -seq: Transmitted indexes to symbols in the constellation. Size = [Nt x T]
% 

%% SAMPLES
% 
% samples
% -H: Channel coefficients. Size = [Nr x Nt x L]
% -s2y: Noise variance
% -d0: Reference distance
% -Ts: Sampling time
% -pathL: Path loss exponent
% -am: Transition probabilities from 0 to 0. Size = [Nt x 1]
% -bm: Transition probabilities from active to 0. Size = [Nt x 1]
% -s2H: Variance of the channel coefficients. Size = [1 x L]
% -Z: Transmitted symbols. Size = [Nt x T]
% -seq: Transmitted indexes to symbols in the constellation. Size = [Nt x T]
% -nest: Number of jumps from 0/1 to 0/1. Size = [2 x 2 x Nt]
% -slice: Slice variable for the blocked sampling approach.
% 
%% CONFIGURATION PARAMETERS
% 
% param
% -gen: Struct with the configuration parameters to generate data.
%   +s2n: Noise variance.
%   +L_true: Channel length for each transmitter. Size = [1 x Nt]
%   +varH: Variance of the channel coefficients.
%   +Nt: Number of transmitters.
%   +burstLength: Mean length of each burst of symbols. Size = [1 x Nt]
%   +burstLengthStdFactor: The std of the length of each burst is (burstLength/burstLengthStdFactor). If infinite, the length is exactly burstLength.
%   +symbol0: Symbol that is sent before (and after) transmission (typically 0).
%   +sparsityH: % of channel coefficients that are set to 0.
% -L: Channel length for inference.
% -flag0: Flag to determine if symbol 0 should be considered.
% -constellation: Vector with the constellation (must exclude 0!).
% -Nr: Number of receiving antennas.
% -T: Length of the symbol sequence.
% -Niter: Number of iterations of the sampler.
% -saveCycle: Save state of the sampler every saveCycle iterations.
% -storeIters: Number of iterations to be stored as local variables.
% -header: Vector containing the header. It can take any length or be empty. Its elements should be strictly positive integers pointing to elements of param.constellation.
% -onOffModel: If true (!=0), a transmitter that becomes inactive cannot start transmitting again.
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
%   +flagParallel: If this field exists and takes value 1, the parallelized version of the code is used (it makes use of OpenMP)
%   +blockNtSize: If a value is specified, transmitters are jointly sampled in blocks of this size. If unspecified, blockNtSize=inf
% -bnp: Struct with the configuration parameters for BNP inference.
%   +betaSlice1: Beta parameter to sample from the slice variable: s=betarnd(betaSlice1,betaSlice2)*c_min.
%   +betaSlice2: Beta parameter to sample from the slice variable: s=betarnd(betaSlice1,betaSlice2)*c_min.
%   +Mini: Initial number of parallel chains
%   +maxMnew: Maximum number of new chains allowed per iteration.
% -colGibbs: Struct with the configuration parameters for collapsed Gibbs algorithm.
%   +Niter: Number of iterations to run.
% -ffbs: Struct with the configuration parameters for FFBS algorithm.
%   +Niter: Number of iterations to run.
% -infer: Struct with the options for the inference algorithm.
%   +symbolMethod: Selects the method to infer the symbols. It can be 'pgas', 'ffbs' or 'gibbs'.
%   +sampleNoiseVar: Flag to indicate if the noise variace should be sampled.
%   +sampleChannel: Flag to indicate if the channel coefficients (and Nt) should be sampled.
%   +sampleVarH: Flag to indicate if the variance of the channel coefficients should be sampled.
%   +simulatedTempering: Flag to indicate if simulated tempering should be used
%   +addArtificialNoise: Flag to indicate if successive noise levels should be used during inference (by adding artificial noise)
% -artifNoise: Struct containing the configuration parameters for the artificial noise
%   +itCycle: How many iterations we should wait before decreasing the artificial noise variance
%   +stepDB: How many dB's we should decrease the artificial noise variance
%   +iniSNR: Initial SNR to be considered
%   +finalSNR: Final SNR to be considered (i.e., the true one)
% -temper: Struct containing the configuration parameters for simulated tempering
%   +pKeep: Probability of keeping the current temperature at each iteration
%   +pNext: Probability of proposing to decrease the noise level
%   +s2yValues: Vector with all the considered noise variances (temperature levels in decreasing order of temperature)
% 
%% HYPERPARAMETERS
% 
% hyper
% -s2h: Define the mean value of the variance of the channel coefficients: E[s2H(r)]=s2h*exp(-lambda*(r-1)).
% -lambda: Define the mean value of the variance of the channel coefficients: E[s2H(r)]=s2h*exp(-lambda*(r-1)).
% -kappa: Defome the variance of the variance of the channel coefficients: Std[s2H(r)]=kappa*E[s2H(r)].
% -alpha: Concentration parameter of the IBP.
% -gamma1: Hyperparameter for variables, bm~Beta(gamma1,gamma2).
% -gamma2: Hyperparameter for variables, bm~Beta(gamma1,gamma2).
% -nu: Hyperparameter for the noise variance, s2y~InvGamma(nu,tau).
% -tau: Hyperparameter for the noise variance, s2y~InvGamma(nu,tau).
% 
