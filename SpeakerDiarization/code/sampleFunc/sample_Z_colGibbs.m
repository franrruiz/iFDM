function [Sest SeqEst nest] = sample_Z_colGibbs(data,samples,hyper,param)
% This function ignores samples.am, samples.bm and samples.H (instead, they
% are integrated out)

Sest = samples.Z;
SeqEst = samples.seq;
nest = samples.nest;

% Remove unused chains
idx0 = find(sum(SeqEst~=0,2)==0);
if(length(idx0)==size(Sest,1))
    idx0(1) = [];
end
Sest(idx0,:) = [];
SeqEst(idx0,:) = [];
nest(:,:,idx0) = [];

% Initialize some quantities
Mest = size(Sest,1);
T = param.T;
constellation = [0 param.constellation];
Q = length(constellation);

% Inference
for it=1:param.colGibbs.Niter
    for t=1:T
        %% Sample the states for instant t
        if(t==1)
            m=1;
            while(m<=Mest)
                % The previous state can be considered 0
                j = 0;
                nest(j+1,(SeqEst(m,t)~=0)+1,m) = nest(j+1,(SeqEst(m,t)~=0)+1,m)-1;
                nest((SeqEst(m,t)~=0)+1,(SeqEst(m,t+1)~=0)+1,m) = nest((SeqEst(m,t)~=0)+1,(SeqEst(m,t+1)~=0)+1,m)-1;
                  
                p = zeros(1,Q);
                fk = zeros(1,Q);
                % Compute transition log-probabilities
                if(SeqEst(m,t+1)==0)
                    % 0 --> 0 --> 0
                    p(1) = log(nest(1,1,m)+2)+log(nest(1,1,m)+1)-log(2+sum(nest(1,:,m)))-log(1+sum(nest(1,:,m)));
                    % 0 --> A --> 0
                    p(2:Q) = -log(Q-1)+log(nest(1,2,m))-log(1+sum(nest(1,:,m)))+log(hyper.gamma1+nest(2,1,m))-log(hyper.gamma1+hyper.gamma2+sum(nest(2,:,m)));
                else
                    % 0 --> 0 --> A
                    p(1) = log(1+nest(1,1,m))+log(nest(1,2,m))-log(1+sum(nest(1,:,m)))-log(2+sum(nest(1,:,m)));
                    % 0 --> A --> A
                    p(2:Q) = -log(Q-1)+log(nest(1,2,m))-log(1+sum(nest(1,:,m)))+log(hyper.gamma2+nest(2,2,m))-log(hyper.gamma1+hyper.gamma2+sum(nest(2,:,m)));
                end
                % Compute log-likelihood
                for k=1:Q
                    Sest(m,t) = constellation(k);
                    fk(k) = mLik_MIMO(data,Sest,samples.s2y,hyper,param);
                end
                % Normalize
                prob = fk+p;
                prob = exp(prob-max(prob));
                prob = prob/sum(prob);
                
                % Sample the current state and update variables
                SeqEst(m,t) = randmult2(prob)-1;
                Sest(m,t) = constellation(SeqEst(m,t)+1);
                nest(j+1,(SeqEst(m,t)~=0)+1,m) = nest(j+1,(SeqEst(m,t)~=0)+1,m)+1;
                nest((SeqEst(m,t)~=0)+1,(SeqEst(m,t+1)~=0)+1,m) = nest((SeqEst(m,t)~=0)+1,(SeqEst(m,t+1)~=0)+1,m)+1;
                
                % Remove empty chains
                if((sum(SeqEst(m,[1:t-1 t+1:end]))==0)&&(Mest>1))
                    SeqEst(m,:) = [];
                    nest(:,:,m) = [];
                    Sest(m,:) = [];
                    Mest = Mest-1;
                else
                    m = m+1;
                end
            end
        elseif(t==T)
            m=1;
            while(m<=Mest)
                % The previous state
                j = (SeqEst(m,t-1)~=0);
                nest(j+1,(SeqEst(m,t)~=0)+1,m) = nest(j+1,(SeqEst(m,t)~=0)+1,m)-1;
                  
                p = zeros(1,Q);
                fk = zeros(1,Q);
                % Compute transition log-probabilities
                if(j==0)
                    % 0 --> 0 --> x
                    p(1) = log(nest(1,1,m)+1)-log(1+sum(nest(1,:,m)));
                    % 0 --> A --> x
                    p(2:Q) = -log(Q-1)+log(nest(1,2,m))-log(1+sum(nest(1,:,m)));
                else
                    % A --> 0 --> x
                    p(1) = log(hyper.gamma1+nest(2,1,m))-log(hyper.gamma1+hyper.gamma2+sum(nest(2,:,m)));
                    % A --> A --> x
                    p(2:Q) = -log(Q-1)+log(hyper.gamma2+nest(2,2,m))-log(hyper.gamma1+hyper.gamma2+sum(nest(2,:,m)));
                end
                % Compute log-likelihood
                for k=1:Q
                    Sest(m,t) = constellation(k);
                    fk(k) = mLik_MIMO(data,Sest,samples.s2y,hyper,param);
                end
                % Normalize
                prob = fk+p;
                prob = exp(prob-max(prob));
                prob = prob/sum(prob);
                
                % Sample the current state and update variables
                SeqEst(m,t) = randmult2(prob)-1;
                Sest(m,t) = constellation(SeqEst(m,t)+1);
                nest(j+1,(SeqEst(m,t)~=0)+1,m) = nest(j+1,(SeqEst(m,t)~=0)+1,m)+1;
                
                % Remove empty chains
                if((sum(SeqEst(m,[1:t-1 t+1:end]))==0)&&(Mest>1))
                    SeqEst(m,:) = [];
                    nest(:,:,m) = [];
                    Sest(m,:) = [];
                    Mest = Mest-1;
                else
                    m = m+1;
                end
            end            
        else
            m=1;
            while(m<=Mest)
                % The previous state
                j = (SeqEst(m,t-1)~=0);
                nest(j+1,(SeqEst(m,t)~=0)+1,m) = nest(j+1,(SeqEst(m,t)~=0)+1,m)-1;
                nest((SeqEst(m,t)~=0)+1,(SeqEst(m,t+1)~=0)+1,m) = nest((SeqEst(m,t)~=0)+1,(SeqEst(m,t+1)~=0)+1,m)-1;
                  
                p = zeros(1,Q);
                fk = zeros(1,Q);
                % Compute transition log-probabilities
                if(j==0)
                    if(SeqEst(m,t+1)==0)
                        % 0 --> 0 --> 0
                        p(1) = log(nest(1,1,m)+2)+log(nest(1,1,m)+1)-log(2+sum(nest(1,:,m)))-log(1+sum(nest(1,:,m)));
                        % 0 --> A --> 0
                        p(2:Q) = -log(Q-1)+log(nest(1,2,m))-log(1+sum(nest(1,:,m)))+log(hyper.gamma1+nest(2,1,m))-log(hyper.gamma1+hyper.gamma2+sum(nest(2,:,m)));
                    else
                        % 0 --> 0 --> A
                        p(1) = log(1+nest(1,1,m))+log(nest(1,2,m))-log(1+sum(nest(1,:,m)))-log(2+sum(nest(1,:,m)));
                        % 0 --> A --> A
                        p(2:Q) = -log(Q-1)+log(nest(1,2,m))-log(1+sum(nest(1,:,m)))+log(hyper.gamma2+nest(2,2,m))-log(hyper.gamma1+hyper.gamma2+sum(nest(2,:,m)));
                    end
                else
                    if(SeqEst(m,t+1)==0)
                        % A --> 0 --> 0
                        p(1) = log(hyper.gamma1+nest(2,1,m))-log(hyper.gamma1+hyper.gamma2+sum(nest(2,:,m)))+log(1+nest(1,1,m))-log(1+sum(nest(1,:,m)));
                        % A --> A --> 0
                        p(2:Q) = -log(Q-1)+log(hyper.gamma1+nest(2,1,m))+log(hyper.gamma2+nest(2,2,m))-log(1+hyper.gamma1+hyper.gamma2+sum(nest(2,:,m)))-log(hyper.gamma1+hyper.gamma2+sum(nest(2,:,m)));
                    else
                        % A --> 0 --> A
                        p(1) = log(hyper.gamma1+nest(2,1,m))-log(hyper.gamma1+hyper.gamma2+sum(nest(2,:,m)))+log(nest(1,2,m))-log(1+sum(nest(1,:,m)));
                        % A --> A --> A
                        p(2:Q) = -log(Q-1)+log(1+hyper.gamma2+nest(2,2,m))+log(hyper.gamma2+nest(2,2,m))-log(1+hyper.gamma1+hyper.gamma2+sum(nest(2,:,m)))-log(hyper.gamma1+hyper.gamma2+sum(nest(2,:,m)));
                    end
                end
                % Compute log-likelihood
                for k=1:Q
                    Sest(m,t) = constellation(k);
                    fk(k) = mLik_MIMO(data,Sest,samples.s2y,hyper,param);
                end
                % Normalize
                prob = fk+p;
                prob = exp(prob-max(prob));
                prob = prob/sum(prob);
                
                % Sample the current state and update variables
                SeqEst(m,t) = randmult2(prob)-1;
                Sest(m,t) = constellation(SeqEst(m,t)+1);
                nest(j+1,(SeqEst(m,t)~=0)+1,m) = nest(j+1,(SeqEst(m,t)~=0)+1,m)+1;
                nest((SeqEst(m,t)~=0)+1,(SeqEst(m,t+1)~=0)+1,m) = nest((SeqEst(m,t)~=0)+1,(SeqEst(m,t+1)~=0)+1,m)+1;
                
                % Remove empty chains
                if((sum(SeqEst(m,[1:t-1 t+1:end]))==0)&&(Mest>1))
                    SeqEst(m,:) = [];
                    nest(:,:,m) = [];
                    Sest(m,:) = [];
                    Mest = Mest-1;
                else
                    m = m+1;
                end
            end
        end
        
        %% Add new chains for instant t
        % This is a truncated Gibbs scheme in which we try to
        % add either 0 or 1 chain. For 2 chains, we would have a complexity
        % increasing with the squared constellation order

        % Probability of adding 0 chains
        p0 = -hyper.alpha/T+mLik_MIMO(data,Sest,samples.s2y,hyper,param);
        % Probabilities of adding 1 chain with state ~=0
        pAct = zeros(1,Q-1);
        pAct(:) = -log(Q-1)-hyper.alpha/T+log(hyper.alpha/T);
        Saux = [Sest; zeros(1,T)];
        for q=2:Q
            Saux(Mest+1,t) = constellation(q);
            pAct(q-1) = pAct(q-1)+mLik_MIMO(data,Saux,samples.s2y,hyper,param);
        end
        
        prob = [p0 pAct];
        prob = exp(prob-max(prob));
        prob = prob/sum(prob);
        
        idxSampled = randmult2(prob);
        % Add a new chain
        if(idxSampled>1)
            Mnew = 1;
            Saux(Mest+1,t) = constellation(idxSampled);
            Sest = Saux;
            SeqEst = [SeqEst; zeros(1,T)];
            SeqEst(Mest+1,t) = idxSampled-1;
            nest = cat(3,nest,zeros(2,2,Mnew));
            if(t~=T)
                nest(1,1,Mest+1:end) = T-2;
                nest(1,2,Mest+1:end) = 1;
                nest(2,1,Mest+1:end) = 1;
            else
                nest(1,1,Mest+1:end) = T-1;
                nest(1,2,Mest+1:end) = 1;
            end
            Mest = Mest+Mnew;
        end
    end
end





