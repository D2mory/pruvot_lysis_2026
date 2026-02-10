function pn = vr_smp(pr)
% function pn = vr_smp(pr)
%
% This function uniformally samples parameters values within defined ranges
%
% Inputs:
% . pr : structure of parameter ranges
% Outputs:
% . pn : one set of sampled parameters 

fs = fieldnames(pr);
pn = struct();
for ctr = 1:numel(fs)
    fn = fs{ctr};
    rngv = pr.(fn);
    pn.(fn) = exp(rand()*log(rngv(2)/rngv(1))+log(rngv(1))); %uniform sampling across the range
end

end