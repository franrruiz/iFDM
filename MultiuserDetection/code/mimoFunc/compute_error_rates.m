function [ADER SER_ALL SER_ACT MMSE vec_ord rot ADER_indiv SER_ALL_indiv SER_ACT_indiv MMSE_indiv] = compute_error_rates(data,samples,hyper,param,flagExhaustiveOrder,flagExhaustiveRot)
% Returns:
% -ADER: Activity detection error rate (i.e., prob that a trasmitter is
%        idle when it is said to be trasnmitting or vice-versa)
% -SER_ALL: Symbol error rate. It is computed as if 0 were another symbol of
%           the constellation
% -SER_ACT: Symbol error rate, conditioned on the actual trasmitter being
%           active
% -MMSE: Mean Square Error (w.r.t. the channel coefficients)
% -vec_ord: Vector containing the order of the estimated transmitters
%           needed to match the true ones
% -rot: Vector containing the coefficients to rotate the constellation in
%       order to match the true symbols
%
% The variables ending with "_indiv" correspond to the individual ADER, SER
% or MSE (for each transmitter)
% 

if(nargin==4)
    flagExhaustiveOrder = 1;
    flagExhaustiveRot = 1;
elseif(nargin==5)
    flagExhaustiveRot = 1;
elseif(nargin<4)
    error('Wrong number of inputs');
end

[Mest T] = size(samples.Z);
Nt = size(data.symbols,1);

thr = min(abs(param.constellation))/10;

if(Mest<Nt)
    samples.H = cat(2,samples.H,zeros(param.Nr,Nt-Mest,param.L));
	samples.Z = [samples.Z; zeros(Nt-Mest,size(samples.Z,2))];
    samples.seq = [samples.seq; zeros(Nt-Mest,size(samples.Z,2))];
	Mest = Nt;
end

if(param.L>size(data.channel,3))
    data.channel = cat(3,data.channel,zeros(param.Nr,Nt,param.L-size(data.channel,3)));
elseif(param.L<size(data.channel,3))
    warning('Cannot compute performance due to wrong choice of param.L');
    ADER = NaN;
    SER_ALL = NaN;
    SER_ACT = NaN;
    MMSE = NaN;
    vec_ord = NaN;
    rot = NaN;
    ADER_indiv = NaN;
    SER_ALL_indiv = NaN;
    SER_ACT_indiv = NaN;
    MMSE_indiv = NaN;
    return;
end

idxNZ = find(data.symbols~=0);
idxNZ_indiv = cell(1,Nt);
for m=1:Nt
    idxNZ_indiv{m} = find(data.symbols(m,:)~=0);
end

if(Mest==Nt)
    % Hay Mest! permutaciones
    permutac = perms(1:Mest);
    ADER = Inf;
    maxI = factorial(Nt);
    if(~flagExhaustiveOrder)
        maxI = 1;
    end
    for i=1:maxI
        ADER_aux = sum(sum((samples.seq(permutac(i,:),:)~=0)~=(data.seq~=0)))/(Nt*T);
        if(ADER_aux<ADER)
            vec_ord = permutac(i,:);
            ADER = ADER_aux;
            SER_ALL = Inf;
            maxM = 4^Nt-1;
            if(~flagExhaustiveRot)
                maxM = 0;
            end
            for m=0:maxM
                rot_aux = exp(1i*pi*de2bi(m,Nt,4)/2);
                Zaux = samples.Z(permutac(i,:),:).*repmat(rot_aux.',1,T);
                Haux = samples.H(:,permutac(i,:),:).*repmat(1./rot_aux,[param.Nr,1,param.L]);
                SER_ALL_aux = sum(abs(Zaux(:)-data.symbols(:))>thr)/(Nt*T);
                if(SER_ALL_aux<SER_ALL)
                    SER_ALL = SER_ALL_aux;
                    SER_ACT = sum(abs(Zaux(idxNZ)-data.symbols(idxNZ))>thr)/length(idxNZ);
                    rot = rot_aux;
                    MMSE = sum(abs(Haux(:)-data.channel(:)).^2)/numel(data.channel);
                    
                    ADER_indiv = sum((Zaux~=0)~=(data.seq~=0),2)/T;
                    SER_ACT_indiv = zeros(Nt,1);
                    SER_ALL_indiv = zeros(Nt,1);
                    MMSE_indiv = zeros(Nt,1);
                    for nt=1:Nt
                        SER_ACT_indiv(nt) = sum(abs(Zaux(nt,idxNZ_indiv{nt})-data.symbols(nt,idxNZ_indiv{nt}))>thr)/length(idxNZ_indiv{nt});
                        SER_ALL_indiv(nt) = sum(abs(Zaux(nt,:)-data.symbols(nt,:))>thr)/T;
                        MMSE_indiv(nt) = sum(sum(abs(Haux(:,nt,:)-data.channel(:,nt,:)).^2))/numel(data.channel(:,nt,:));
                    end
                end
            end
        end
    end
elseif(Mest>Nt)
    combinac = combnk(1:Mest,Nt);
    ADER = Inf;
    for j=1:size(combinac,1)
        permutac = perms(combinac(j,:));
        idxNotUsed = setdiff(1:Mest,combinac(j,:));
        for i=1:factorial(Nt)
            ADER_aux = (sum(sum((samples.seq(permutac(i,:),:)~=0)~=(data.seq~=0)))+sum(sum(samples.seq(idxNotUsed,:)~=0)))/(Nt*T);
            if(ADER_aux<ADER)
                vec_ord = permutac(i,:);
                ADER = ADER_aux;
                SER_ALL = Inf;
                maxM = 4^Nt-1;
                if(~flagExhaustiveRot)
                    maxM = 0;
                end
                for m=0:maxM
                    rot_aux = exp(1i*pi*de2bi(m,Nt,4)/2);
                    Zaux = samples.Z(permutac(i,:),:).*repmat(rot_aux.',1,T);
                    Haux = samples.H(:,permutac(i,:),:).*repmat(1./rot_aux,[param.Nr,1,param.L]);
                    SER_ALL_aux = (sum(abs(Zaux(:)-data.symbols(:))>thr)+sum(sum(samples.seq(idxNotUsed,:)~=0)))/(Nt*T);
                    if(SER_ALL_aux<SER_ALL)
                        SER_ALL = SER_ALL_aux;
                        SER_ACT = sum(abs(Zaux(idxNZ)-data.symbols(idxNZ))>thr)/length(idxNZ);
                        rot = rot_aux;
                        MMSE = sum(abs(Haux(:)-data.channel(:)).^2)/numel(data.channel);
                        
                        ADER_indiv = sum((Zaux~=0)~=(data.seq~=0),2)/T;
                        SER_ACT_indiv = zeros(Nt,1);
                        SER_ALL_indiv = zeros(Nt,1);
                        MMSE_indiv = zeros(Nt,1);
                        for nt=1:Nt
                            SER_ACT_indiv(nt) = sum(abs(Zaux(nt,idxNZ_indiv{nt})-data.symbols(nt,idxNZ_indiv{nt}))>thr)/length(idxNZ_indiv{nt});
                            SER_ALL_indiv(nt) = sum(abs(Zaux(nt,:)-data.symbols(nt,:))>thr)/T;
                            MMSE_indiv(nt) = sum(sum(abs(Haux(:,nt,:)-data.channel(:,nt,:)).^2))/numel(data.channel(:,nt,:));
                        end
                    end
                end
            end
        end
    end
else
    error(['Mest y Nt no son comparables: Mest=' num2str(Mest) ' y Nt=' num2str(Nt)]);
end

