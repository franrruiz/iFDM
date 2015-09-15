function [Z]=FB(X,s2x,Z, Phi, A)
[D T]=size(X);
Q=size(A,1);
M=size(Z,1);
Xhat=cell(1,M);
fw= zeros(Q,T);
bw= zeros(Q,T);

 %% Forward-Backward

 
 for m=1:M
     suma=zeros(D,T);
     for q=1:Q-1
         Zq=(Z([1:m-1 m+1:end],:)==q);
         suma=suma+reshape(Phi(:,q,[1:m-1 m+1:end]),D,M-1)*Zq;
     end
     Xhat{m}=X-suma;
     
     %%Forward
     t=1;
     fw(:,t)=-1/(2*s2x)*diag((repmat(Xhat{m}(:,t),1,Q)-[zeros(D,1) Phi(:,:,m)])'*(repmat(Xhat{m}(:,t),1,Q)-[zeros(D,1) Phi(:,:,m)]))+log(A(1,:,m)');
     fw(:,t)=exp(fw(:,t)-max(fw(:,t)));
     fw(:,t)= fw(:,t)/sum( fw(:,t));
     for t=2:T
         fw(:,t)=-1/(2*s2x)*diag((repmat(Xhat{m}(:,t),1,Q)-[zeros(D,1) Phi(:,:,m)])'*(repmat(Xhat{m}(:,t),1,Q)-[zeros(D,1) Phi(:,:,m)]))+log(A(:,:,m)'*fw(:,t-1));
         fw(:,t)=exp(fw(:,t)-max(fw(:,t)));
         fw(:,t)= fw(:,t)/sum(fw(:,t));
     end
     
     %Backward
     t=T;
     Z(m,t)=mnrnd(1,fw(:,t))*[0:Q-1]';
     for t=T-1:-1:1
        bw(:,t)=fw(:,t).*A(:,Z(m,t+1)+1,m); 
        bw(:,t)= bw(:,t)/sum(bw(:,t));
        
        Z(m,t)=mnrnd(1,bw(:,t))*[0:Q-1]';
     end
     
 end
 
         
     
     