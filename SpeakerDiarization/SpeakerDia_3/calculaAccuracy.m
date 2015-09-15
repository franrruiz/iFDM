function [acc cadenas_ord] = calculaAccuracy(cadenas,devices)
% Acc = 1 - (sum_t (sum_m |y^true_tm - y^est_tm|))) / 2*(sum_t y^obs_t)

Mest = size(cadenas,1);
nDev = size(devices,1);

if(Mest==nDev)
    % Hay Mest! permutaciones
    permutac = perms(1:Mest);
    acc = -Inf;
    for i=1:factorial(Mest)
        acc_aux = 1-(sum(sum(abs(cadenas(permutac(i,:),:)~=devices))))/(2*sum(sum(devices)));
        if acc_aux>acc
            cadenas_ord = cadenas(permutac(i,:),:);
            acc = acc_aux;
        end
    end
elseif(Mest<nDev)
    cadenas = [cadenas; zeros(nDev-Mest,size(cadenas,2))];
    Mest = nDev;
    permutac = perms(1:Mest);
    acc = -Inf;
    for i=1:factorial(Mest)
        acc_aux = 1-(sum(sum(abs(cadenas(permutac(i,:),:)~=devices))))/(2*sum(sum(devices)));
        if acc_aux>acc
            cadenas_ord = cadenas(permutac(i,:),:);
            acc = acc_aux;
        end
    end
elseif(Mest>nDev)
    combinac = combnk(1:Mest,nDev);
    acc = -Inf;
    for j=1:size(combinac,1)
        permutac = perms(combinac(j,:));
        for i=1:factorial(nDev)
            acc_aux = 1-(sum(sum(abs(cadenas(permutac(i,:),:)~=devices)))+sum(sum(abs(cadenas(setdiff(1:Mest,permutac(i,:)),:)))))/(2*sum(sum(devices)));
            if acc_aux>acc
                cadenas_ord = [cadenas(permutac(i,:),:); cadenas(setdiff(1:Mest,permutac(i,:)),:)];
                acc = acc_aux;
            end
        end
    end
else
    error(['Mest y nDev no son comparables: Mest=' num2str(Mest) ' y nDev=' num2str(nDev)]);
end
