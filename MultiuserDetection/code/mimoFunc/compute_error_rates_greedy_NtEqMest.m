function [vec_ord rot desp] = compute_error_rates_greedy_NtEqMest(data,samples,hyper,param)
% Returns:
% -vec_ord: Vector containing the order of the estimated transmitters
%           needed to match the true ones
% -rot: Vector containing the coefficients to rotate the constellation in
%       order to match the true symbols
% -desp: Vector containing the delay of the inferred symbols w.r.t. the true ones
%

[Mest T] = size(samples.Z);
Nt = size(data.symbols,1);
despMax = max(param.L,20);
thr = min(abs(param.constellation))/10;
thrSER = 0.1;

if(Mest~=Nt)
    error('This function only works for Mest==Nt');
end

%% Choose the order

vec_ord = zeros(1,Nt);
rot = zeros(1,Nt)+1i*zeros(1,Nt);
desp = zeros(1,Nt);
alreadyChosen = [];
notChosen = 1:Nt;

vec_ord_inf = zeros(1,Mest);
rot_inf = zeros(1,Mest);
desp_inf = zeros(1,Mest);

% We match the inferred transmitters with the true ones
for m=1:Mest
    flagChosen = 0;
    min_ser_m = inf;
    ll=-despMax;
    while(ll<=despMax && ~flagChosen)
        % Shifted replica of Z(m,:)
        Zm = [zeros(1,max(0,ll)) samples.Z(m,max(1-ll,1):min(T-ll,T)) zeros(1,max(0,-ll))];
        r = 0;
        while(r<=3 && ~flagChosen)
            % Rotate the constellation by a factor of pi*r/2
            Zm_rot = Zm*exp(1i*r*pi/2);
            % Compute the SER vs all transmitters
            ser_vs_all = sum(abs(repmat(Zm_rot,Nt,1)-data.symbols)>thr,2)/T;
            % If the SER of the m-th Tx has improved
            if(min(ser_vs_all(notChosen))<min_ser_m)
                % Find those tx's for which the SER is below threshold
                idx_all = find(ser_vs_all(notChosen)<thrSER);
                if(isempty(idx_all))
                    [min_ser_m idxSelected] = min(ser_vs_all(notChosen));
                    vec_ord_inf(m) = notChosen(idxSelected);
                    rot_inf(m) = exp(1i*r*pi/2);
                    desp_inf(m) = ll;
                else
                    [min_ser_m idxSelected] = min(ser_vs_all(notChosen));
                    vec_ord_inf(m) = notChosen(idxSelected);
                    rot_inf(m) = exp(1i*r*pi/2);
                    desp_inf(m) = ll;
                    alreadyChosen = [alreadyChosen notChosen(idxSelected)];
                    notChosen(notChosen==notChosen(idxSelected)) = [];
                    flagChosen = 1;
                end
            end
            r = r+1;
        end
        ll = ll+1;
    end
end

% In case some transmitters have not been yet assigned, perform a more greedy approach
% -First:
selectedOnes = unique(vec_ord_inf);
for ii=1:length(selectedOnes)
    m = selectedOnes(ii);
    % If exactly one inferred chain has been assigned to the m-th Tx, add m to the list
    if(sum(vec_ord_inf==m)==1)
        alreadyChosen = union(alreadyChosen,m);
        notChosen(notChosen==m) = [];
    % If >1 inferred chains have been assigned to the m-th Tx and the error is below thrSER
    % for one of them, then de-assign the rest of them
    elseif(sum(vec_ord_inf==m)>1 && sum(alreadyChosen==m)==1)
        mask_inf = (vec_ord_inf==m);
        for maux=1:Mest
            Zmaux = [zeros(1,max(0,desp_inf(maux))) samples.Z(maux,max(1-desp_inf(maux),1):min(T-desp_inf(maux),T)) zeros(1,max(0,-desp_inf(maux)))];
            if(sum(abs(Zmaux-data.symbols(vec_ord_inf(maux),:))>thr)/T<thrSER)
                mask_inf(maux) = 0;
            end
        end
        vec_ord_inf(mask_inf) = 0;
    % If >1 inferred chains have been assigned to the m-th Tx and the error is above thrSER
    % for all of them, add the best one to the list and de-assign the rest of them
    elseif(sum(vec_ord_inf==m)>1 && sum(alreadyChosen==m)==0)
        mask_inf = (vec_ord_inf==m);
        list_maux = find(vec_ord_inf==m);
        seraux = inf;
        winner = -1;
        for maux=list_maux
            Zmaux = [zeros(1,max(0,desp_inf(maux))) samples.Z(maux,max(1-desp_inf(maux),1):min(T-desp_inf(maux),T)) zeros(1,max(0,-desp_inf(maux)))];
            if(sum(abs(Zmaux-data.symbols(m,:))>thr)/T<seraux)
                winner = maux;
            end
        end
        mask_inf(winner) = 0;
        vec_ord_inf(mask_inf) = 0;
        alreadyChosen = union(alreadyChosen,m);
        notChosen(notChosen==m) = [];
    end
end
% -Second: repeat the initial step more greedily for the still-not-assigned Tx
for m=1:Mest
    flagChosen = (vec_ord_inf(m)~=0);
    min_ser_m = inf;
    ll=-despMax;
    while(ll<=despMax && ~flagChosen)
        % Shifted replica of Z(m,:)
        Zm = [zeros(1,max(0,ll)) samples.Z(m,max(1-ll,1):min(T-ll,T)) zeros(1,max(0,-ll))];
        r = 0;
        while(r<=3 && ~flagChosen)
            % Rotate the constellation by a factor of pi*r/2
            Zm_rot = Zm*exp(1i*r*pi/2);
            % Compute the SER vs all transmitters
            ser_vs_all = sum(abs(repmat(Zm_rot,Nt,1)-data.symbols)>thr,2)/T;
            % If the SER of the m-th Tx has improved
            if(min(ser_vs_all(notChosen))<min_ser_m)
                [min_ser_m idxSelected] = min(ser_vs_all(notChosen));
                vec_ord_inf(m) = notChosen(idxSelected);
                rot_inf(m) = exp(1i*r*pi/2);
                alreadyChosen = [alreadyChosen notChosen(idxSelected)];
                notChosen(notChosen==notChosen(idxSelected)) = [];
                flagChosen = 1;
            end
            r = r+1;
        end
        ll = ll+1;
    end
end

%% Now that we have computed the right order for the inferred Tx, compute vec_ord and rot
for m=1:Nt
    auxIdx = find(vec_ord_inf==m);
    if(length(auxIdx)~=1)
        error('This should not happen');
    end
    vec_ord(m) = auxIdx;
    rot(m) = rot_inf(auxIdx);
    desp(m) = desp_inf(auxIdx);
end    
