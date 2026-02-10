function xt = gauss_trf(x)
% function xt = gauss_trf(x)

[~,~,ranks] = unique(x);  % ranks go from 1 to n
u = (ranks-0.5)/numel(x); % center ranks in interval [0,n]
xt = norminv(u,0,1);    % guassian values (mean 0, std 1)

end