function theta = dirichletrnd(alpha, N)
% SAMPLE_DIRICHLET Sample N vectors from Dir(alpha(1), ..., alpha(k))
% theta = sample_dirichlet(alpha, N)
% theta(i,j) = i'th sample of theta_j, where theta ~ Dir

% We use the method from p. 482 of "Bayesian Data Analysis", Gelman etal.

k = size(alpha,1);
scale = 1; % arbitrary



if (k == 1)
    theta = gamrnd(repmat(alpha,N,1), scale*ones(N,size(alpha,2)), N, size(alpha,2)); % size [N*k]
    S = sum(theta,2);
    theta = theta ./ repmat(S, 1, size(alpha,2));   
elseif (N == k)
    % we want as many samples as rows
    theta = gamrnd(alpha, scale*ones(size(alpha)));
    S = sum(theta,2);
    theta = theta ./ repmat(S,1,size(alpha,2));
else
    error('Check')
end

% k = length(alpha);
% theta = zeros(N, k);
% scale = 1; % arbitrary?
% for i=1:k
%   theta(:,i) = gamrnd(alpha(i), scale, N, 1);
% end
% S = sum(theta,2);
% theta = theta ./ repmat(S, 1, k);