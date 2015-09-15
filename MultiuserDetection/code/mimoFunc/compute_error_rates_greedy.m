function [ADER SER_ALL SER_ACT MMSE vec_ord rot desp ADER_indiv SER_ALL_indiv SER_ACT_indiv MMSE_indiv] = compute_error_rates_greedy(data,samples,hyper,param)
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
% -desp: Vector containing the delay of the inferred symbols w.r.t. the true ones
%
% The variables ending with "_indiv" correspond to the individual ADER, SER
% or MSE (for each transmitter)
% 

[Mest T] = size(samples.Z);
Nt = size(data.symbols,1);
Nr = size(data.obs,1);
thr = min(abs(param.constellation))/10;

if(Mest==Nt)
    [vec_ord rot desp] = compute_error_rates_greedy_NtEqMest(data,samples,hyper,param);
elseif(Mest<Nt)
    samples.H = cat(2,samples.H,zeros(param.Nr,Nt-Mest,param.L));
	samples.Z = [samples.Z; zeros(Nt-Mest,T)];
    samples.seq = [samples.seq; zeros(Nt-Mest,T)];
	Mest = Nt;
    [vec_ord rot desp] = compute_error_rates_greedy_NtEqMest(data,samples,hyper,param);
elseif(Mest>Nt)
    data.channel = cat(2,data.channel,zeros(param.Nr,Mest-Nt,size(data.channel,3)));
	data.symbols = [data.symbols; zeros(Mest-Nt,T)];
    data.seq = [data.seq; zeros(Mest-Nt,T)];
	Nt = Mest;
    [vec_ord rot desp] = compute_error_rates_greedy_NtEqMest(data,samples,hyper,param);
end

if(param.L>size(data.channel,3))
    data.channel = cat(3,data.channel,zeros(param.Nr,Nt,param.L-size(data.channel,3)));
elseif(param.L<size(data.channel,3))
    samples.H = cat(3,samples.H,zeros(param.Nr,Mest,size(data.channel,3)-param.L));
    param.L = size(data.channel,3);
end

ADER_indiv = zeros(Nt,1);
SER_ALL_indiv = zeros(Nt,1);
SER_ACT_indiv = zeros(Nt,1);
MMSE_indiv = zeros(Nt,1);

badSymbols = 0;
countSymbols = 0;
for n=1:Nt
    m = vec_ord(n);  % m refers to the inferred chain and n to the actual Tx
    r = rot(n);
    ll = desp(n);
    Zm = [zeros(1,max(0,ll)) r*samples.Z(m,max(1-ll,1):min(T-ll,T)) zeros(1,max(0,-ll))];
    
    ADER_indiv(n) = sum((Zm~=0)~=(data.symbols(n,:)~=0))/T;
    SER_ALL_indiv(n) = sum(abs(Zm-data.symbols(n,:))>thr)/T;
    idxAct = (data.seq(n,:)~=0);
    badSymbols = badSymbols+sum(abs(Zm(idxAct)-data.symbols(n,idxAct))>thr);
    countSymbols = countSymbols+sum(idxAct);
    SER_ACT_indiv(n) = sum(abs(Zm(idxAct)-data.symbols(n,idxAct))>thr)/sum(idxAct);
    
    if(ll==0)
        MMSE_indiv(n) = sum(sum(abs((1/r)*samples.H(:,m,:)-data.channel(:,n,:)).^2))/numel(data.channel(:,n,:));
    elseif(ll>0)
        trueH = cat(3,zeros(Nr,1,ll),data.channel(:,n,:));
        auxH = (1/r)*samples.H(:,m,:);
        countCut = 0;
        for j=1:ll
            if(sum(sum(abs(trueH(:,:,end))))==0)
                trueH(:,:,end) = [];
                countCut = countCut+1;
            end
        end
        if(countCut<ll)
            auxH = cat(3,auxH,zeros(Nr,1,ll-countCut));
        end
        MMSE_indiv(n) = sum(sum(abs(auxH-trueH).^2))/numel(trueH);
    elseif(ll<0)
        auxH = cat(3,zeros(Nr,1,-ll),(1/r)*samples.H(:,m,:));
        trueH = data.channel(:,n,:);
        countCut = 0;
        for j=1:-ll
            if(sum(sum(abs(auxH(:,:,end))))==0)
                auxH(:,:,end) = [];
                countCut = countCut+1;
            end
        end
        if(countCut<-ll)
            trueH = cat(3,trueH,zeros(Nr,1,-ll-countCut));
        end
        MMSE_indiv(n) = sum(sum(abs(auxH-trueH).^2))/numel(trueH);
    end
end

ADER = mean(ADER_indiv);
SER_ALL = mean(SER_ALL_indiv);
SER_ACT = badSymbols/countSymbols;
MMSE = mean(MMSE_indiv);

