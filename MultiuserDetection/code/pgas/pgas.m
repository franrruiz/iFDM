function [X_PG] = pgas(Y,Z,Nt,C,Q,L,H,sy2,a,b,N_PF,N_PG,M)

[Nr T] = size(Y);

flagPG = 1;
if(isempty(Z))
    flagPG = 0;
    xc      =   zeros(Nt,T);      %   The particle that we condition on
else
    xc      =   Z;
end

X_PG    =   zeros(Nt,M,T);    %   Stores MCMC samples of trajectories

% Note that X_PG(:,m,t) is a vector that describes the transmitted symbols
% at time t, according to the m'th sample of the sequence of symbols. 

for m = 1 : M
    if ((m == 1) && (~flagPG))
        % The particle filter used to generate the first trajectory makes
        % use of more particles
        N   =   N_PF; 
    else
        % than the PGAS kernel
        N   =   N_PG;
    end
    % We now initialize three important variables:
    Xt      =   zeros(Nt,N,T);  %   Stores particles within PGAS
    a_ind   =   zeros(N,T);     %   Stores the ancestor indices
    W       =   zeros(N,T);     %   Stores the particle weights
    
    % In our first version we also store the last L elements 
    % in the particle trajectories. 
    X0_hist =   zeros(Nt,N,L); % Contains x_{t-L+1:t}^i, i=1,...,N
    X1_hist =   zeros(Nt,N,L); % Contains x_{t-L:t-1}^i, i=1,...,N
    
    % NOTE 1: Xt also contains particles but to find out how these are
    % connected across time we need to consider the ancestor indices a_ind.
    % NOTE 2: it is not necessary to store X0_hist and X1_hist since the same
    % information can be obtained from a_ind and Xt. We can probably speed
    % up the current implementation by not introducing the two variables. 
    % NOTE 3: within the PGAS algorithm, we use X0_hist to compute the
    % weights for particles 1 to N-1 and X1_hist to compute weights that we
    % use to sample the ancestors for particle N. 
    
    % In order to use built-in functions in Matlab it is convenient for us to 
    % produce matrices containing the am and bm coefficients
    A   =   repmat(a,1,N);
    B   =   repmat(b,1,N);
    
    for t = 1 : T 
        % To story away x_{t-L:t-1}^i for i = 1, 2, ..., N:
        X1_hist     =   X0_hist;
        if t == 1
            % At time t = 1 we sample the states from the prior at time 1.
            % We know that all transmitters were passive at time 0     [Line 1]   
            Xt(:,:,t)   =   binornd(ones(Nt,N),1-A).*unidrnd(Q,Nt,N);

            % Note that the particles do not yet have any ancestors that 
            % we can sample. However, if m > 1 the N'th particle should be the
            % particle that we condition on                            [Line 2]
            if ((m ~= 1) || (flagPG))
                Xt(:,N,t) = xc(:,t);
            end

            % For convenience, we store the particle trajectories here:
            X0_hist      =   cat(3,zeros(Nt,N,L-1),Xt(:,:,t));

        else
            % At later times we first sample the ancestors, using multinomial 
            % resampling, to obtain a(:,t) = ind                       [Line 5]
            [valnul,ind] = histc(rand(N,1), [0 cumsum(W(:,t-1)')]);

            % Then we propagate the selected particles from time t-1 to t, 
            % to obtain Xt(:,:,t)                                      [Line 5]
            Act         =   Xt(:,ind,t-1)>0;
            Xt(:,:,t)   =   (Act.*binornd(ones(Nt,N),1-B)+...
               (1-Act).*binornd(ones(Nt,N),1-A)).*unidrnd(Q,Nt,N);
           % Note: if we wish to improve the update rate for complex scenarios
           % we should probably try to improve how we propagate particles from
           % time t-1 to t. Specifically, it seems highly inefficient to
           % sample using "unidrnd(Q,Nt,N)" when Q and Nt are large. 

            % Within the PGAS algorithm, the ancestor index for the N'th
            % particle is sampled differently:
            if ((m ~= 1) || (flagPG))
                % Set the N'th particle to the particle that we condition on
                Xt(:,N,t) =   xc(:,t);                             %   [Line 6]

                % The ancestor probabilities (the weights) are obtained by
                % multiplying three factors. The objective here is to calculate
                % these factors in order to obtain what we denote w_a. [Line 7]


                % To compute the factor that depends on y, we currently use a 
                % double for-loop (certainly, this can be improved upon):
%                 logWY    =   zeros(N,L-1);
%                 for tau     =   t : min(t+L-2,T)
%                     % Symbols at different times contribute to the same measurement
%                     % and r defines the delay.
%                     % For small delays r, the receivers observe "symbols" 
%                     % transmitted after time t, thus described by the sequence 
%                     % that we condition on. For larger delays, the symbols 
%                     % are instead described by particles at time t-1.
%                     % The smallest delay for which the particle sequence in the 
%                     % interval 1 to t-1 has an impact is when tau-r=t-1.
%                     r1 = 0:tau-t;
%                     r2 = tau-t+1:min(L-1,tau-1);
%                     
%                     if(length(r1)==1)
%                         permV = [2 3 1];
%                     else
%                         permV = [1 3 2];
%                     end
%                     
%                     aux   = sum(mtimesx(H(:,:,r1+1),permute(C(xc(:,tau-r1)+1),permV)),3);
%                     Ypred = aux(:,ones(1,N)) + sum(mtimesx(H(:,:,r2+1),C(X1_hist(:,:,end+tau-r2-t+1)+1)),3);
% 
%                     Ydiff                =   Ypred - repmat(Y(:,tau),1,N);
%                     logWY(:,tau-t+1)     =   sum(-abs(Ydiff).^2/sy2,1);
%                 end
%                 logWY   =   sum(logWY,2);
                logWY = zeros(N,1);
                if(L>1)
                    Act_anc = zeros(N,L-1);
                    Act_anc(:,L-1) = 1:N;
                    r = L-2;
                    while((r>=1)&&(t+r-L>=0))
                         idxSurv = (Act_anc(:,r+1)>0);
                         aux = unique(a_ind(Act_anc(idxSurv,r+1),t+r-L+1));
                         Act_anc(1:length(aux),r) = aux;
                         r = r-1;
                    end
                    logWY = acc_loop(Nt,Nr,N,T,L,t,length(C),C,Y,Xt,xc,H,a_ind,Act_anc,sy2);
                end
                WY      =   exp(logWY-max(logWY));

                % Another factor represents the transition probabilities between
                % particle i at time t-1 and the particle that we condition on at
                % time t. To compute these probabilities we need to distinguish
                % between four events for each symbol and particle i:
                % 0    ->  0,   passive remains passive, Prob = am
                % 0    -> (>0), passive becomes active,  Prob = (1-am)/Q
                % (>0) -> (>0), active remains active,   Prob = (1-bm)/Q
                % (>0) -> 0,    active becomes passive,  Prob = bm
                WZ_mat  =   zeros(Nt,N);    % Transition probability for each element
                X1      =   Xt(:,:,t-1);    % The particles at time t-1
                Act1    =   X1>0;
                Act0    =   repmat(xc(:,t)>0,1,N);
                % We can now go through the four cases and compute transition
                % probabilities for each symbol and particle:
                WZ_mat  =   (1-Act1).*(1-Act0).*A + (1-Act1).*Act0.*(1-A)/Q ...
                    + Act1.*Act0.*(1-B)/Q +Act1.*(1-Act0).*B;
                logWZ   =   sum(log(WZ_mat),1); % Log-transition probabilities for each particle
                WZ      =   exp(logWZ-max(logWZ))';

                % Finally, we can compute the weights of interest
                w_a     =   W(:,t-1).*WY.*WZ;
                w_a     =   w_a/sum(w_a);
                if(sum(isnan(w_a)>0))
                    w_a = log(W(:,t-1))+logWY+logWZ.';
                    w_a = exp(w_a-max(w_a));
                    w_a = w_a/sum(w_a);
                end
                % from which we generate the N'th ancestor             [Line 8]
                ind(N) = find(rand(1) < cumsum(w_a),1); 
            end
            % We have now computed all the ancestor indices
            a_ind(:,t)  =   ind;        % Stores results of     [Lines 5 and 8]

            % We can also store away the particles recent histories,  ~[Line 9] 
            X0_hist     =   cat(3,X1_hist(:,ind,2:end),Xt(:,:,t));
        end


        % Compute the importance weights for all the particles, W(:,t).
        % First we compute the expected measurements for different particles
        % by considering the particles and their histories     [Lines 3 and 10]
        Ypred   =   zeros(Nr,N);
        for r = 0 : min(L-1,t-1)
             Ypred   =   Ypred   +   H(:,:,r+1)*C(X0_hist(:,:,end-r)+1);
        end        
        Ydiff       =   Ypred - repmat(Y(:,t),1,N);
        logW        =   sum(-abs(Ydiff).^2/sy2,1)';
        W(:,t)      =   exp(logW-max(logW));
        W(:,t)      =   W(:,t)/sum(W(:,t));
    end % This marks the end of the for-loop over time, t.


    % We can now compute the generated trajectories from the ancestor indices
    %                                                                  [Line 9]
    ind = a_ind(:,T);
    for(t = T-1:-1:1)
        Xt(:,:,t) = Xt(:,ind,t);
        ind = a_ind(ind,t);
    end
    % Two comments:
    %   --  Previously we only stored the most recent parts of the particle
    %   trajectories, which is more efficient as long L<<T. Here we construct
    %   the entire trajectories from Xt and a_ind.
    %   --  It may be slightly faster to first select one trajectory at random
    %   (see below) and then only construct that tractory, but this loop is
    %   anyway very fast so it probably does not make much difference. 

    % Finally we select one trajectory at random                      [Line 12]
    J = find(rand(1) < cumsum(W(:,T)),1);
    % that we can store away 
    X_PG(:,m,:)     =   Xt(:,J,:);
    xc              =   reshape(Xt(:,J,:),Nt,T);

    %%
end


