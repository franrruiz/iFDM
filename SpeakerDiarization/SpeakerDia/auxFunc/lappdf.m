function y  = lappdf(x, mu, sigma)
%LAPRND generate i.i.d. laplacian random number drawn from laplacian distribution
%   with mean mu and standard deviation sigma. 
%   mu      : mean
%   sigma   : standard deviation


% Generate Laplacian noise
b = sigma / sqrt(2);
y = 1/(2*b)*exp(-abs(x-mu)/b);
