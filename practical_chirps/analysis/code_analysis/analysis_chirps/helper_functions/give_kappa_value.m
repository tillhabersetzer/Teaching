function kappa = give_kappa_value(cov,label,sensors)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (computes rank of matrix with svd)
% the following lines detect the location of the first large 'cliff' in the 
% singular value spectrum of the grads and mags
% kappa  : number of spatial components in covariance matrix (rank)
%
% cov    :   covariance matrix
% label  : label vector of covariance matrix
% sensors: cell array for specific sensors e.g. {'megmag','megplanar'}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N     = length(sensors);
kappa = zeros(1,N);

for n = 1:N
    idx      = ft_chantype(label,sensors{n}); % logical for specific sensors
    [~,s,~]  = svd(cov(idx,idx));
    d        = -diff(log10(diag(s))); 
    d        = d./std(d);
    kappa(n) = find(d>4,1,'first');   
end

end