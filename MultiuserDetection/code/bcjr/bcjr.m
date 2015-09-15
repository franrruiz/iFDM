function qt = bcjr(Y,Haux,s2y,ptrans,A)

T = size(Y,2);

Yopt = Haux*A.';
Qtot = size(ptrans,1);

fw=zeros(Qtot,T);
bw=zeros(Qtot,T);
fw(:,1)=ptrans(1,:).*exp(-sum(abs(repmat(Y(:,1),1,Qtot)-Yopt).^2,1)/(s2y));
fw(:,1)=fw(:,1)/sum(fw(:,1));
bw(:,T)=1;
for t=2:T
    fw(:,t)=(ptrans'*fw(:,t-1)).*exp(-sum(abs(repmat(Y(:,t),1,Qtot)-Yopt).^2,1)/(s2y))';
    fw(:,t)=fw(:,t)/sum(fw(:,t));
    %aux2 = b_log(:,T-t+2)+log(bw(:,T-t+2));
    %bw(:,T-t+1)= a{l}*(b(:,T-t+2).*bw(:,T-t+2));
    bw(:,T-t+1)= ptrans*(bw(:,T-t+2).*exp(-sum(abs(repmat(Y(:,T-t+2),1,Qtot)-Yopt).^2,1)/(s2y))');
    bw(:,T-t+1)=bw(:,T-t+1)/sum(bw(:,T-t+1));
end
qt = fw.*bw;
qt=qt./repmat(sum(qt,1),Qtot,1);
