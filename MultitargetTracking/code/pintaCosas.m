close all;

figure;
plot(squeeze(data.states(:,1,:))')
xlabel('t');
ylabel('x_1');

figure;
plot(squeeze(samples.Z(:,1,:))')
xlabel('t');
ylabel('Inferred x_1');

colores = 'rgcbmkyrgcbmky';
figure;
for nt=1:size(data.states,1)
    idxNZ = find(squeeze(data.states(nt,1,:))~=0);
    plot(squeeze(data.states(nt,1,idxNZ))',squeeze(data.states(nt,2,idxNZ))','Color',colores(nt));
    hold on;
end
plot(data.sensors(:,1),data.sensors(:,2),'bo');
xlabel('x_1');
ylabel('x_2');

figure;
for nt=1:size(samples.Z,1)
    idxNZ = find(squeeze(samples.Z(nt,1,:))~=0);
    plot(squeeze(samples.Z(nt,1,idxNZ))',squeeze(samples.Z(nt,2,idxNZ))','Color',colores(nt));
    hold on;
end
plot(data.sensors(:,1),data.sensors(:,2),'bo');
xlabel('x_1');
ylabel('x_2');

