function pn = model_WEI_smp(pr)
% function pn = model_WEI_smp(pr)
%
% Randomly samples a parameter set from given parameter value ranges
%
% Inputs:
% . pr : parameter value ranges
%
% Outputs: 
% . pn : sampled parameter set

fs = fieldnames(pr);
pn = struct();
for ctr = 1:numel(fs)
    fn = fs{ctr};
    range = pr.(fn);
    pn.(fn) = exp(rand()*log(range(2)/range(1))+log(range(1)));
end

end
