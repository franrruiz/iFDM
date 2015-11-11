function [X_PG] = pgas(Y,Z,Nt,C,Q,L,H,sy2, s2X,a,b,N_PF,N_PG,M)

[D T] = size(Y);

flagPG = 1;
if(isempty(Z))
    flagPG = 0;
    xc      =   zeros(Nt,T);      %   The particle that we condition on
else
    xc      =   double(Z~=0);
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
    An   =   repmat(1-a,1,N);
    Bn   =   repmat(1-b,1,N);
    tic
    for t = 1 : T 
        % To story away x_{t-L:t-1}^i for i = 1, 2, ..., N:
        X1_hist     =   X0_hist;
        if t == 1
            % At time t = 1 we sample the states from the prior at time 1.
            % We know that all transmitters were passive at time 0     [Line 1]   
            Xt(:,:,t)   =  (rand(Nt,N)<An);
        

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
            %slow part
            Act         =   Xt(:,ind,t-1)>0;
            Xt(:,:,t)   =   (Act.*(rand(Nt,N)<Bn)) +((1-Act).*(rand(Nt,N)<An)); 
           
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

                % This factor represents the transition probabilities between
                % particle i at time t-1 and the particle that we condition on at
                % time t. To compute these probabilities we need to distinguish
                % between four events for each symbol and particle i:
                % 0    ->  0,   passive remains passive, Prob = am
                % 0    -> (>0), passive becomes active,  Prob = (1-am)/Q
                % (>0) -> (>0), active remains active,   Prob = (1-bm)/Q
                % (>0) -> 0,    active becomes passive,  Prob = bm
                WZ_mat  =   zeros(Nt,N);    % Transition probability for each element
                X1      =   Xt(:,:,t-1);    % The particles at time t-1
                Act1    =   X1~=0;
                Act0    =   repmat(xc(:,t)~=0,1,N);
                % We can now go through the four cases and compute transition
                % probabilities for each symbol and particle:
                WZ_mat  =   (1-Act1).*(1-Act0).*A + (1-Act1).*Act0.*(An) ...
                    + Act1.*Act0.*(Bn) +Act1.*(1-Act0).*B;
                logWZ   =   sum(log(WZ_mat),1); % Log-transition probabilities for each particle
                WZ      =   exp(logWZ-max(logWZ))';

                % Finally, we can compute the weights of interest
                w_a     =   W(:,t-1).*WZ;
                w_a     =   w_a/sum(w_a);
                if(sum(isnan(w_a)>0))
                    w_a = log(W(:,t-1))+logWZ.';
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
        %Ypred   =   zeros(D,N);
        for mm = 1 : N
             Waux   =  squeeze(H(Xt(:,mm,t)~=0,:));
             Sy=Waux*Waux'+sy2/s2X*eye(size(Waux,1));
             Sy= eye(D)- Waux'*(Sy\Waux);
             logW(mm) = 1/2*log(det(Sy/sy2))-1/(2*sy2)*Y(:,t)'*Sy*Y(:,t);
        end        
        W(:,t)      =   exp(logW-max(logW));
        W(:,t)      =   W(:,t)/sum(W(:,t));
    end % This marks the end of the for-loop over time, t.


    % We can now compute the generated trajectories from the ancestor indices
    %                                                                  [Line 9]
    ind = a_ind(:,T);
    for t = T-1:-1:1
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

    %% Sampling X for the selected particle
    for t=1:T
        Waux   =  squeeze(H(X_PG(:,m,t)~=0,:));
        if ~isempty(Waux)
        Sx=Waux*Waux'+sy2/s2X*eye(size(Waux,1));
        LambdaX=Waux*Y(:,t);
        X_PG(X_PG(:,m,t)~=0,m,t)= mvnrnd2(Sx\LambdaX,sy2*inv(Sx),1);
        end
    end
end


