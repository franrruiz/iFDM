function [Sest SeqEst] = pgas_main(data,samples,hyper,param)

%% Extract parameters from structs
Nt = size(samples.H,2);
M = param.pgas.Niter;
N_PF = param.pgas.N_PF;
N_PG = param.pgas.N_PG;
returnN = param.pgas.returnNsamples;

if(param.flag0)
    constellation = [0 param.constellation];
else
    constellation = [param.constellation];
end

%% Some configuration parameters
if(isempty(samples.seq))
    flagPG = 0;
else
    flagPG = 1;
end

flagParallel = 0;
if(isfield(param.pgas,'flagParallel'))
    if(param.pgas.flagParallel)
        flagParallel = 1;
    end
end

blockNtSize = inf;
if(isfield(param.pgas,'blockNtSize'))
    blockNtSize = param.pgas.blockNtSize;
end

%% Call the pgas function
if(Nt<=blockNtSize)
    % Option 1: All transmitters jointly sampled
    if(flagParallel)
        X_PG = pgas_C_parallel(Nt,param.Nr,N_PF,N_PG,M,param.T,param.L,length(constellation),constellation,data.obs,flagPG,samples.seq,samples.H,samples.s2y,samples.am,samples.bm,param.header,length(param.header),param.pgas.particles,param.onOffModel);
    else
        X_PG = pgas_C(Nt,param.Nr,N_PF,N_PG,M,param.T,param.L,length(constellation),constellation,data.obs,flagPG,samples.seq,samples.H,samples.s2y,samples.am,samples.bm,param.header,length(param.header),param.pgas.particles,param.onOffModel);
    end

    %% Return last obtained samples in a [Nt x T x M] matrix
    SeqEst = permute(X_PG(:,M-returnN+1:M,:),[1 3 2]);
    idxN0 = find(SeqEst~=0);
    Sest = zeros(Nt,param.T,returnN);
    Sest(idxN0) = constellation(SeqEst(idxN0)+param.flag0);
else
    % Option 2: Not all transmitters jointly sampled
    nRounds = ceil(Nt/blockNtSize);
    Sest = samples.Z;
    SeqEst = samples.seq;
    for m=1:M
        idxNt = randperm(Nt);
        for nn=1:nRounds
            consideredTx = idxNt((nn-1)*blockNtSize+1:min(nn*blockNtSize,length(idxNt)));
            notTx = setdiff(1:Nt,consideredTx);
            
            pseudoObs = data.obs;
            pseudoSeq = SeqEst(consideredTx,:);
            pseudoH = samples.H(:,consideredTx,:);
            pseudoam = samples.am(consideredTx);
            pseudobm = samples.bm(consideredTx);
            for ll=1:param.L
                pseudoObs = pseudoObs-samples.H(:,notTx,ll)*[zeros(length(notTx),ll-1) Sest(notTx,1:param.T-(ll-1))];
            end

            if(flagParallel)
                X_PG = pgas_C_parallel(length(consideredTx),param.Nr,N_PF,N_PG,1,param.T,param.L,length(constellation),constellation,pseudoObs,flagPG,pseudoSeq,pseudoH,samples.s2y,pseudoam,pseudobm,param.header,length(param.header),param.pgas.particles,param.onOffModel);
            else
                X_PG = pgas_C(length(consideredTx),param.Nr,N_PF,N_PG,1,param.T,param.L,length(constellation),constellation,pseudoObs,flagPG,pseudoSeq,pseudoH,samples.s2y,pseudoam,pseudobm,param.header,length(param.header),param.pgas.particles,param.onOffModel);
            end
            
            SeqEst(consideredTx,:) = permute(X_PG(:,1,:),[1 3 2]);
            idxN0 = find(SeqEst~=0);
            Sest = zeros(Nt,param.T);
            Sest(idxN0) = constellation(SeqEst(idxN0)+param.flag0);
        end
    end
end

