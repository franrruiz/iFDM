function generate_data(Nt,s2y,T)
Ts = 0.4;
sy2 = s2y;
Ns = 25;
pathL = 2;
W = 800;
s2u = 1;
d0 = 1;
Pt = 50;
s2vIni = 0.01;
Gx = [1 0 Ts 0; 0 1 0 Ts; 0 0 1 0; 0 0 0 1];
Gu = [Ts^2/2 0; 0 Ts^2/2; Ts 0; 0 Ts];

pos = 0:round(W/(sqrt(Ns)-1)):W;
[a b] = meshgrid(pos);
sensors = [a(:) b(:)];

state=zeros(Nt,4,T);
for nt=1:Nt
    Tini = 1;%randi([1 round(T/2)],1,1);
    Tend = T;%min(Tini+round(T/2)-1,T);
    state(nt,:,Tini)= [W*rand(2,1); sqrt(s2vIni)*randn(2,1)];
    flagTxDone = 0;
    while(~flagTxDone)
        for t=Tini+1:Tend
            state(nt,:,t) = Gx*squeeze(state(nt,:,t-1))'+sqrt(s2u)*Gu*randn(2,1);
        end
        if(~(sum(state(nt,1,:)<-W/5)>0 || sum(state(nt,1,:)>W+W/5)>0 || sum(state(nt,2,:)<-W/5)>0 || sum(state(nt,2,:)>W+W/5)>0))
            flagTxDone = 1;
        end
    end
end

Ptx=zeros(Ns,T);
for t=1:T
    d=zeros(Ns,1);
    for nt=1:Nt
        if(state(nt,1,t)~=0)
        	d= d+1./(((sensors(:,1)-state (nt,1,t)).^2 +(sensors(:,2)-state (nt,2,t)).^2).^(pathL/2));
        end
    end
    Ptx(d>0,t)= 10^(Pt/10)*d0^pathL*d(d>0);
end

obs = Ptx+sqrt(s2y)*randn(Ns,T);

figure,plot(squeeze(state(:,1,:))', squeeze(state(:,2,:))','x')
hold on, plot(sensors(:,1),sensors(:,2),'o')
figure,plot((obs(:,:))')

save(['dataTracking_Nt' num2str(Nt) '_s2y' num2str(s2y) '_T' num2str(T) '.mat']);
