function pn = wei_smp(pr)
% function pn = wei_smp(pr)
%
% Inputs:
% . pr : parameter ranges
%
% Outputs : 
% . pn : parameter value sampled

fs = fieldnames(pr);
pn = struct();
for ctr = 1:numel(fs)
    fn = fs{ctr};
    range = pr.(fn);
    pn.(fn) = exp(rand()*log(range(2)/range(1))+log(range(1))); %uniform sampling
end
end