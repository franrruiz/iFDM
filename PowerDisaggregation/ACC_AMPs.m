function ACC_AMPs(Q, Tini, Tend)
Nd=8;
BASEDIR1=['AMPs/resultsPGAS/M' num2str(Nd) '_T' num2str(Tini) '_' num2str(Tend) '_Q' num2str(Q)];
load([BASEDIR1 '/Final.mat'],'data','init','samples','samplesAll', 'LLH', 'M_EST');

ACC=zeros(1,2000);
cad_ord=cell(1,2000);
devices=data.devices;
for it=1:2000
    [ACC(it) cad_ord{it}]= calculaAccuracy(samplesAll{it}.Z,devices);
end

BASEDIR1=['AMPs/resultsPGAS'];
save([BASEDIR1 '/Acc_T' num2str(Tini) '_' num2str(Tend) '_Q' num2str(Q) '.mat'],'ACC','cad_ord','devices');