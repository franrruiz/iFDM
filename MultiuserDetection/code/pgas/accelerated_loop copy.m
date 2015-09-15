%% (1) Build the matrix Act_anc (surely there is a more efficient way)
tic;
Act_anc = zeros(N,L-1);
Act_anc(:,L-1) = 1:N;
r = L-2;
while((r>=1)&&(t+r-L>=0))
     idxSurv = (Act_anc(:,r+1)>0);
     aux = unique(a_ind(Act_anc(idxSurv,r+1),t+r-L+1));
     Act_anc(1:length(aux),r) = aux;
     r = r-1;
end

logWY2 = acc_loop(Nt,Nr,N,T,L,t,length(C),C,Y,Xt,xc,H,a_ind,Act_anc,sy2);

toc
disp(num2str(sum(abs(logWY-logWY2))));


%% (2) Replacement of the for loop in L98-121 of the pgas.m file
% logWY2 corresponds to logWY in pgas.m (here it has a different name to compare both)
% logWY2 = zeros(N,1);
% for tau=t:min(t+L-2,T)
%     Ypred = Y(:,tau);   % Ypred is initialized to the tau-th observation
%     
%     % Remove the contribution of the fixed particle to instant tau
%     for r=0:tau-t
%         Ypred = Ypred-H(:,:,r+1)*C(xc(:,tau-r)+1).';
%     end
%     % Compute the contribution of ancestors to instant tau
%     minr = tau-t+1;
%     maxr = min(L-1,tau-1);
%     ancContrib = zeros(Nr,N,L-1);
%     for r=maxr:-1:minr
%         ii = 1;
%         while((ii<=N)&&(Act_anc(ii,tau-t+L-r)>0))
%             n = Act_anc(ii,tau-t+L-r);
%             ancContrib(:,n,tau-t+L-r) = H(:,:,r+1)*C(1+Xt(:,n,tau-r)).';
%             if(r<=maxr-1)
%                 nAnt = a_ind(n,tau-r);
%                 ancContrib(:,n,tau-t+L-r) = ancContrib(:,n,tau-t+L-r)+ancContrib(:,nAnt,tau-t+L-r-1);
%             end
%             ii = ii+1;
%         end
%     end
%     
%     % For each particle, remove the contribution of ancestors to instant
%     % tau and compute the likelihood (add it to logWY)
%     for n=1:N
%         logWY2(n) = logWY2(n)-sum(abs(Ypred-ancContrib(:,n,L-1)).^2)/sy2;
%     end
%     
%     % ------------- Sanity check (compare with the previous implementation)
%     % This part can be safely removed
%     r1 = 0:tau-t;
%     r2 = tau-t+1:min(L-1,tau-1);
%     if(length(r1)==1)
%         permV = [2 3 1];
%     else
%         permV = [1 3 2];
%     end
%     aux1 = sum(mtimesx(H(:,:,r1+1),permute(C(xc(:,tau-r1)+1),permV)),3);
%     aux2 = sum(mtimesx(H(:,:,r2+1),C(X1_hist(:,:,size(X1_hist,3)+tau-r2-t+1)+1)),3);
%     disp(['aux1: ' num2str(sum(abs(Ypred+(aux1-Y(:,tau)))))]);         % Should be 0
%     disp(['aux2: ' num2str(sum(sum(abs(ancContrib(:,:,L-1)-aux2))))]); % Should be 0
% end







