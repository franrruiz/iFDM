function R = mvnrnd2(mu,sigma,N)
% 
% Returns the same as:
%   R = mvnrnd(mu,sigma,N)
% but performs an additional test to ensure that the covariance matrix is
% symmetric
%

cte = max(max(abs(sigma-sigma')./abs(sigma)));
if((cte>0) )%&& (cte<1e-2)) % should be 0
    sigma = diag(diag(sigma))+(triu(sigma,1)'+triu(sigma,1)+tril(sigma,-1)'+tril(sigma,-1))/2;
elseif(cte>=1e-2) % matrix is not symmetric
    save('sigma.mat','sigma');
    error('mvnrnd2:SigmaNotSymmetric','The covariance matrix should be symmetric');
end

R = mvnrnd(mu,sigma,N);
