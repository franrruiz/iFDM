sload('/Users/ivalera/Desktop/uc3m/SourceSeparation/Tracking/Final.mat')
% figurapdf(5, 4)
% print -dpdf Tracking.pdf
Dtable=zeros(3,300);
jj=[3 2 1];
%jj=[1 2 3];
for ii=1:3
    Dtable(ii,data.states(ii,1,:)~=0)= sqrt((squeeze(data.states(ii,1,data.states(ii,1,:)~=0))-squeeze(samples.Z(jj(ii),1,data.states(ii,1,:)~=0))).^2+(squeeze(data.states(ii,2,data.states(ii,1,:)~=0))-squeeze(samples.Z(jj(ii),2,data.states(ii,1,:)~=0))).^2);
    
end
        
mean(Dtable(1,Dtable(2,:)~=0))
mean(Dtable(2,Dtable(2,:)~=0))
mean(Dtable(3,Dtable(3,:)~=0))

mean(Dtable(Dtable~=0))