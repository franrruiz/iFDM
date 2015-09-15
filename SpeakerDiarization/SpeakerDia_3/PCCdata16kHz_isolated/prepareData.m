load('data/data4.mat','speakers');
speakers=100*squeeze(speakers(:,1,:));
[T aux1 aux2]=size(speakers);
W=rand(15,10);
s2y=0.3^2;
obs = (speakers*W)'+sqrt(s2y)*randn(10,T);

save('data/dataF4.mat','speakers','obs','W','s2y');
