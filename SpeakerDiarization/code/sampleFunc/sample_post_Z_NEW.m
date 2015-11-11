function [Sest SeqEst nest out] = sample_post_Z(data,samples,hyper,param)

%% Choose the desired method to sample Z:

% For PGAS
if(strcmp(param.infer.symbolMethod,'pgas'))
    [Sest, SeqEst] = pgas_main_matlab(data,samples,hyper,param);
    % Update the number of transitions, nest
    nest = zeros(2,2,size(Sest,1));
    for m=1:size(Sest,1)
        % From 0 to 0
        nest(1,1,m) = sum([0 SeqEst(m,1:end-1)]==0 & SeqEst(m,:)==0);
        % From 0 to active
        nest(1,2,m) = sum([0 SeqEst(m,1:end-1)]==0 & SeqEst(m,:)~=0);
        % From active to 0
        nest(2,1,m) = sum([0 SeqEst(m,1:end-1)]~=0 & SeqEst(m,:)==0);
        % From active to active
        nest(2,2,m) = sum([0 SeqEst(m,1:end-1)]~=0 & SeqEst(m,:)~=0);
    end
    out = 0;
    
% For EP
elseif(strcmp(param.infer.symbolMethod,'ep'))
    [Sest SeqEst nest out] = ep_main(data,samples,hyper,param);
    
% For Collapsed Gibbs sampling
elseif(strcmp(param.infer.symbolMethod,'gibbs'))
    [Sest SeqEst nest] = sample_Z_colGibbs(data,samples,hyper,param);
    out = 0;
    
% For FFBS
%  - Gaussian prior over x_tm
elseif(strcmp(param.infer.symbolMethod,'ffbs_gauss'))
    [Sest SeqEst nest] = sample_Z_FFBS_gauss(data,samples,hyper,param);
    out = 0;
%  - Laplace prior over x_tm
elseif(strcmp(param.infer.symbolMethod,'ffbs_laplace'))
    [Sest SeqEst nest] = sample_Z_FFBS_laplace(data,samples,hyper,param);
    out = 0;
    
% Else, give an error
else
    error(['Invalid sampling method: ' param.infer.symbolMethod]);
end

