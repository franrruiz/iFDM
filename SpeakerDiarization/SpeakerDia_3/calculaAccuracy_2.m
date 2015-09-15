function [acc, cadenas_ord, idx_Ord] = calculaAccuracy_2(cadenas,devices)
% Acc = 1 - (sum_t (sum_m |y^true_tm - y^est_tm|))) / 2*(sum_t y^obs_t)

Mest = size(cadenas,1);
nDev = size(devices,1);
M=max(Mest,nDev);

Mcomb=zeros(M, M);

if(Mest>nDev)
   devices=[devices; zeros(Mest-nDev,size(cadenas,2))];
elseif(Mest<nDev)
   cadenas=[cadenas; zeros(nDev-Mest,size(cadenas,2))]; 
end
accP=zeros(1,M);
idx_Ord=zeros(1,M);
for i=1:M
    
    for j=1:M
        if sum(sum(devices(i,:)))~=0
            Mcomb(i,j) = 1-(sum(sum(abs(cadenas(j,:)~=devices(i,:)))))/(2*sum(sum(devices(i,:))));
        else
            Mcomb(i,j) = 1-(sum(sum(abs(cadenas(j,:)~=devices(i,:)))));
        end
    end
    [accP(i), idx_Ord(i)]=max(Mcomb(i,:));
    ii=i;
    idxDif=1:M;
    while sum(idx_Ord(ii)==idx_Ord([1:ii-1 ii+1:end]))>0
        idxjj=find(idx_Ord==idx_Ord(ii));
        jj=setdiff(idxjj,ii);
        [val]=max(accP(ii),accP(jj));
        if val==accP(jj)
           idxDif=setdiff(idxDif,idx_Ord(jj));
           accP(ii) = max(Mcomb(ii,idxDif));
           idx_Ord(ii)=find(Mcomb(ii,:)==accP(ii));
        else
           idxDif=setdiff(idxDif,idx_Ord(ii));
           accP(jj) = max(Mcomb(jj,idxDif));
           try
           idx_Ord(jj)=find(Mcomb(jj,:)==accP(jj));
           catch
               %disp('prueba')
               aux=find(Mcomb(jj,:)==accP(jj));
               idx_Ord(jj)=aux(2);
           end
           ii=jj;
        end
    end
       
end

cadenas_ord = cadenas(idx_Ord,:);
acc=mean(accP);