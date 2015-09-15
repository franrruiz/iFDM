%function ACC_REDD(H,Q)
Nd=5;
T=960;
Tsubsample=1000;
BASEDIR1=['PCCdata16kHz_isolated/resultsPGAS/S' num2str(Nd) '_T' num2str(T) '_Tsub' num2str(Tsubsample)];
load([BASEDIR1 '/Final.mat'],'data','init','samples','samplesAll', 'LLH', 'M_EST');
Mmax=max(M_EST);
cad_ord=cell(1,2000);
devices=double(data.speakers~=0)';
for it=1:2000
    cadenas=double(samplesAll{it}.Z~=0);
    
end

[ACC cad_ord]= calculaAccuracy(cadenas,devices);
BASEDIR1=['PCCdata16kHz_isolated/resultsPGAS'];
save([BASEDIR1 '/ADER2_S' num2str(Nd) '_T' num2str(T) '_Tsub' num2str(Tsubsample) '.mat'],'ACC','cad_ord','devices');