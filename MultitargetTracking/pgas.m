function [X_PG] = pgas(Y,sensors, Area, Gx, Gu, Ts, d0, pathL, Ptx, s2u,Z,Nt,L,sy2,a,b,N_PF,N_PG,M,s2vIni)

[Nr T] = size(Y);

flagPG = 1;
if(isempty(Z))
    flagPG = 0;
    xc      =   zeros(Nt,4,T);      %   The particle that we condition on
else
    xc      =   Z;
end

X_PG    =   zeros(Nt,M,4,T);    %   Stores MCMC samples of trajectories

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
    Xt      =   zeros(Nt,N,4,T);  %   Stores particles within PGAS
    a_ind   =   zeros(N,T);     %   Stores the ancestor indices
    W       =   zeros(N,T);     %   Stores the particle weights
    
    
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
        if t == 1
            % At time t = 1 we sample the states from the prior at time 1.
            % We know that all transmitters were passive at time 0     [Line 1]   
            Xt(:,:,:,t)   =  repmat(rand(Nt,N)<An,[1,1,4]).*cat(3,Area*rand([Nt N 2]), sqrt(s2vIni)*randn([Nt N 2]));
        

            % Note that the particles do not yet have any ancestors that 
            % we can sample. However, if m > 1 the N'th particle should be the
            % particle that we condition on                            [Line 2]
            if ((m ~= 1) || (flagPG))
                Xt(:,N,:,t) = xc(:,:,t);
            end

        else
            % At later times we first sample the ancestors, using multinomial 
            % resampling, to obtain a(:,t) = ind                       [Line 5]
            [valnul,ind] = histc(rand(N,1), [0 cumsum(W(:,t-1)')]);

            % Then we propagate the selected particles from time t-1 to t, 
            % to obtain Xt(:,:,t)                                      [Line 5]
            %slow part
            Act         =   double(Xt(:,ind,1,t-1)~=0);
            aux=zeros(Nt,N,4);
            for itm=1:Nt               
                aux(itm,:,:)=(Gx*squeeze(Xt(itm,ind,:,t-1))'+ sqrt(s2u)*Gu*randn([2 N]))';                
            end
            aux2=cat(3,Area*rand([Nt N 2]), sqrt(s2vIni)*randn([Nt N 2]));
            Xt(:,:,:,t)   =   repmat((Act.*(rand(Nt,N)<Bn)),[1,1,4]).*aux;
            Xt(:,:,:,t)   =   Xt(:,:,:,t)+repmat((1-Act).*(rand(Nt,N)<An),[1,1,4]).*aux2;
           
           % Note: if we wish to improve the update rate for complex scenarios
           % we should probably try to improve how we propagate particles from
           % time t-1 to t. Specifically, it seems highly inefficient to
           % sample using "unidrnd(Q,Nt,N)" when Q and Nt are large. 

            % Within the PGAS algorithm, the ancestor index for the N'th
            % particle is sampled differently:
            if ((m ~= 1) || (flagPG))
                % Set the N'th particle to the particle that we condition on
                Xt(:,N,:,t) =   xc(:,:,t);                             %   [Line 6]


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
                X1      =   Xt(:,:,:,t-1);    % The particles at time t-1
                Act1    =   X1(:,:,1)~=0;
                Act0    =   repmat(xc(:,1,t)~=0,1,N);
                % We can now go through the four cases and compute transition
                % probabilities for each symbol and particle:
                aux=zeros(Nt,N);
                for itm=1:Nt
                    Xaux=repmat(xc(itm,:,t),N,1)'-Gx*squeeze(X1(itm,:,:))';
                    aux(itm,:)=(-0.5*log(2*pi*s2u)-0.5*(Xaux(3,:)/Ts).^2/s2u)+ (-0.5*log(2*pi*s2u)-0.5*(Xaux(4,:)/Ts).^2/s2u);
                end
                aux2=log(1/Area*1/Area)+repmat((-0.5*log(2*pi*s2vIni)-0.5*(xc(:,3,t)).^2/s2vIni)+ (-0.5*log(2*pi*s2vIni)-0.5*(xc(:,4,t)).^2/s2vIni), 1,N);
                WZ_mat(((1-Act1).*(1-Act0))==1)=log(A(((1-Act1).*(1-Act0))==1)); 
                WZ_mat(((1-Act1).*Act0)==1)=log(An(((1-Act1).*Act0)==1))+aux2(((1-Act1).*Act0)==1);
                WZ_mat((Act1.*Act0)==1)=log(Bn((Act1.*Act0)==1))+aux((Act1.*Act0)==1); 
                WZ_mat((Act1.*(1-Act0))==1)=log(B((Act1.*(1-Act0))==1));
                logWZ   =   sum((WZ_mat),1); % Log-transition probabilities for each particle

                % Finally, we can compute the weights of interest
                %if(sum(isnan(w_a)>0))
                    w_a = log(W(:,t-1))+logWZ.';
                    w_a = exp(w_a-max(w_a));
                    w_a = w_a/sum(w_a);
                %end
                % from which we generate the N'th ancestor             [Line 8]
                ind(N) = find(rand(1) < cumsum(w_a),1); 
            end
            % We have now computed all the ancestor indices
            a_ind(:,t)  =   ind;        % Stores results of     [Lines 5 and 8]
            
        end


        % Compute the importance weights for all the particles, W(:,t).
        % First we compute the expected measurements for different particles
        % by considering the particles and their histories     [Lines 3 and 10]
        Ptot=zeros(Nr,N);
        d=zeros(Nr,N);
        for itm=1:Nt 
            idxX=find(Xt(itm,:,1,t)~=0);
            Xaux=Xt(itm,idxX,:,t);
            RR=length(idxX);
            d(:,idxX)= d(:,idxX)+1./((repmat(sensors(:,1),1,RR)-repmat(Xaux(:,:,1),Nr,1)).^2 +(repmat(sensors(:,2),1,RR)-repmat(Xaux(:,:,2),Nr,1)).^2).^(pathL/2);
        end  
        idxX=find(sum(d,1)~=0);
        Ptot(:,idxX)=   10^(Ptx/10)*d0^pathL*d(:,idxX);
        Ydiff       =   Ptot - repmat(Y(:,t),1,N);
        logW        =   sum(-abs(Ydiff).^2/sy2,1)';
        W(:,t)      =   exp(logW-max(logW));
        W(:,t)      =   W(:,t)/sum(W(:,t));
%W(:,t) = 1/N;
    end % This marks the end of the for-loop over time, t.


    % We can now compute the generated trajectories from the ancestor indices
    %                                                                  [Line 9]
    ind = a_ind(:,T);
    for(t = T-1:-1:1)
        Xt(:,:,:,t) = Xt(:,ind,:,t);
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
    X_PG(:,m,:,:)     =   Xt(:,J,:,:);
    %xc              =   reshape(Xt(:,J,:),Nt,T);

    %%
end


