function samples = sample_newsticks(data,samples,hyper,param)

if(~param.infer.sampleChannel)
    return;
end

%% (1) Sample the sticks
flagEmpty = 0;
if((size(samples.Z,1)==1)&&(sum(sum(samples.seq))==0))
    c = 1;
    flagEmpty = 1;
else
    c = 1-samples.am;
end
cmin = min(c);

varargin = cell(1,3);
varargin{1} = hyper.alpha;
varargin{2} = param.T;
options = optimset('Display','none','TolX',1e-12);

itSlice = 1;
flagA = 1;
Aaux = cmin;
Puntos = [];
A0no0 = [];
while((flagA==1) && (itSlice<=param.bnp.maxMnew))
    domain = [-inf,log(Aaux)];
    varargin{3} = Aaux;
    % Find the maximum of the log-pdf in the log-domain
    a = fminsearchbnd(@(z)lc_lpdf_neg(z,varargin),log(Aaux)-10,-Inf,log(Aaux),options);
    % Sample new stick length in the log-domain using ARS
    b = log(Aaux);
    try
        Puntos = Puntos(find(Puntos<b));
        [Aaux Puntos] = ars(@lc_lpdf, a-10, b, domain, 1, Puntos, varargin);
        Aaux = exp(Aaux);
    catch err
        disp(err);
        exit;
    end
    
    if(Aaux>samples.slice)
        A0no0 = [A0no0; Aaux];
    else
        flagA = 0;
    end
    itSlice = itSlice+1;
end

Mnew = length(A0no0);
if(Mnew==0)
    return;
end

% Add the new sticks to the existing ones
samples.am = [samples.am; 1-A0no0];
if(flagEmpty)
    samples.am(1) = [];
end

%% (2) Sample new variables bm from prior
samples.bm = [samples.bm; betarnd(hyper.gamma1,hyper.gamma2,Mnew,1)];
if(flagEmpty)
    samples.bm(1) = [];
end

%% (3) Sample new channel coefficients H from prior
% First, detect wheter the constellation is complex-valued
flagComplex = 1;
if(sum(abs(imag(param.constellation)))<1e-5*sum(abs(real(param.constellation))))
    flagComplex = 0;
end
Hnew = zeros(param.Nr,Mnew,param.L);
for ll=1:param.L
    Hnew(:,:,ll) = sqrt(samples.s2H(ll))*(randn(param.Nr,Mnew,1)+1i*randn(param.Nr,Mnew,1));
end
samples.H = cat(2,samples.H,Hnew);
if(flagEmpty)
    samples.H(:,1,:) = [];
end
if(~flagComplex)
    samples.H = real(samples.H);
end

%% (4) Extend the representation of the symbol matrix Z and update nest
samples.Z = [samples.Z; zeros(Mnew,param.T)];
samples.seq = [samples.seq; zeros(Mnew,param.T)];
samples.nest = cat(3,samples.nest,zeros(2,2,Mnew));
samples.nest(1,1,end-Mnew+1:end) = param.T;
if(flagEmpty)
    samples.Z(1,:) = [];
    samples.seq(1,:) = [];
    samples.nest(:,:,1) = [];
end
