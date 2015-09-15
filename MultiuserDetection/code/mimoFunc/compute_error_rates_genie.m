function [ADER SER_ALL SER_ACT MMSE vec_ord rot ADER_indiv SER_ALL_indiv SER_ACT_indiv MMSE_indiv] = compute_error_rates_genie(data,samples,hyper,param)
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

[Mest T] = size(samples.Z);
Nt = size(data.symbols,1);

thr = inf;
for q1=1:length(param.constellation)
    for q2=1:length(param.constellation)
        if(q1~=q2)
            auxThr = abs(param.constellation(q1)-param.constellation(q2));
            if(auxThr<thr)
                thr = auxThr;
            end
        end
    end
end
thr = thr/10;

if(Mest~=Nt)
    error('This is a genie-aided function: The size of samples.Z should coincide with the size of data.symbols');
end

% Set MSE to 0
MMSE = 0;

% Set vec_ord and rot
vec_ord = 1:Nt;
rot = ones(1,Nt);

% Compute SER_ALL, SER_ACT, ADER
idxNZ = find(data.symbols~=0);
SER_ALL = sum(abs(samples.Z(:)-data.symbols(:))>thr)/(Nt*T);
SER_ACT = sum(abs(samples.Z(idxNZ)-data.symbols(idxNZ))>thr)/length(idxNZ);
ADER = sum((samples.seq(:)~=0)~=(data.seq(:)~=0))/(Nt*T);

% Compute tx-individual performance
ADER_indiv = sum((samples.Z~=0)~=(data.seq~=0),2)/T;
SER_ACT_indiv = zeros(Nt,1);
SER_ALL_indiv = zeros(Nt,1);
MMSE_indiv = zeros(Nt,1);
for nt=1:Nt
    idxNZ_indiv = find(data.symbols(nt,:)~=0);
    
    SER_ACT_indiv(nt) = sum(abs(samples.Z(nt,idxNZ_indiv)-data.symbols(nt,idxNZ_indiv))>thr)/length(idxNZ_indiv);
    SER_ALL_indiv(nt) = sum(abs(samples.Z(nt,:)-data.symbols(nt,:))>thr)/T;
end
