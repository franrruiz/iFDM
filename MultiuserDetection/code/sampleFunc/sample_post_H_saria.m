function [H post_mean] = sample_post_H_saria(data,samples,hyper,param)

if(~param.infer.sampleChannel)
    H = samples.H;
    post_mean = samples.H;
    return;
end

% Obtain parameters from the structs
Nr = param.Nr;
[Nt T] = size(samples.Z);
L = param.L;

% Detect wheter the constellation is complex-valued
flagComplex = 1;
if(sum(abs(imag(param.constellation)))<1e-5*sum(abs(real(param.constellation))))
    flagComplex = 0;
end

% Build matrix S containing the symbols and their shifted replicas
S = zeros(T,Nt*L);
for ll=1:L
    S(ll:T,ll:L:L*Nt) = samples.Z(:,1:T-ll+1).';
end

% Build the posterior covariance matrix in complex and real forms
v_1 = 1./samples.s2H;  % These are the inverse prior variances
v_1 = repmat(v_1,1,Nt);

Gam = (diag(v_1)+(1/samples.s2y)*(S'*S))\eye(Nt*L);  % Posterior covariance
Vxx = 0.5*real(Gam);   % Covariance for the real and imaginary parts (they are both the same)
Vxy = -0.5*imag(Gam);  % Covariance between the real and the imaginary part

% Compute the posterior mean
post_mean = Gam*S'*(data.obs.')/samples.s2y;  % Size = [Nt*L x Nr]

% Sample the channel coefficients
if(flagComplex)
    % Covariance matrix
    SS = [Vxx Vxy; Vxy' Vxx];
    SS = 0.5*(SS+SS.');
    try
        U = chol(SS);
    catch e
        [E,Lambda] = eig(SS);
        if(min(diag(Lambda))<0)
            error('Sigma must be positive semi-definite.');
        end
        U = sqrt(Lambda)*E';
    end
    % Sample from the multivariate Gaussian distribution
    muestr = randn(Nr,2*Nt*L)*U + [real(post_mean.') imag(post_mean.')];
    h = muestr(:,1:Nt*L)+1i*muestr(:,Nt*L+1:2*Nt*L);  % Size = [Nr x Nt*L]
else
    % Covariance matrix
    SS = Vxx;
    SS = 0.5*(SS+SS.');
    try
        U = chol(SS);
    catch e
        [E,Lambda] = eig(SS);
        if(min(diag(Lambda))<0)
            error('Sigma must be positive semi-definite.');
        end
        U = sqrt(Lambda)*E';
    end
    % Sample from the multivariate Gaussian distribution
    muestr = randn(Nr,Nt*L)*U + real(post_mean.');
    h = muestr;   % Size = [Nr x Nt*L]
end

% Reshape to fit matrix H
H = permute(reshape(h,[Nr,L,Nt]),[1 3 2]);  % Size = [Nr x Nt x L]
post_mean = permute(reshape(post_mean.',[Nr,L,Nt]),[1 3 2]);  % Size = [Nr x Nt x L]
