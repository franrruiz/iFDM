function ACC_REDD(H,Q)
Nd=6;
Tini= 1;
Tend= 2880;
BASEDIR1=['REDD/resultsPGAS/House' num2str(H) '_M' num2str(Nd) '_T' num2str(Tini) '_' num2str(Tend) '_Q' num2str(Q)];
load([BASEDIR1 '/Final.mat'],'data','init','samples','samplesAll', 'LLH', 'M_EST');

ACC=zeros(1,2000);
cad_ord=cell(1,2000);
devices=data.devices;
for it=1:2000
    [ACC(it) cad_ord{it}]= calculaAccuracy(samplesAll{it}.Z,devices);
end

BASEDIR1=['REDD/resultsPGAS'];
save([BASEDIR1 '/Acc_House' num2str(H) '_Q' num2str(Q) '.mat'],'ACC','cad_ord','devices');