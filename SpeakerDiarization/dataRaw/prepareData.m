close all

speakers=zeros(20*16000,2,5);

for sp=1:5
    t=0;
    for sen=1:4
        [data f]=wavread(['s' num2str(sp) '_' num2str(sen) '.wav']);
        t=t+poissrnd((5.*rand)*16000);
        if t>20*16000
            disp('missing sentence')
        end
        speakers(t:t+size(data,1)-1,:,sp)=data;
        t=t+size(data,1)-1;
    end
end

 plot(squeeze(speakers(1:500:end,1,:)))
 subplot(1,2,1);
 imshow(squeeze(speakers(1:1000:end,1,:)~=0))
 save('data3.mat','speakers')