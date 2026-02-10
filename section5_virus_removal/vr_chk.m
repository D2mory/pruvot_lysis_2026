function chk = vr_chk(x,rng)
% function chk = vr_chk(x,rng)
%
% Function insuring that all variables steady states are within realistic ranges
%
% Inputs:
% . x : variables steady state values
% . rng : ranges of realistic values
% Outputs : 
% . chk : logical (==0 if values are not within desired ranges)

x=x(:);
chk=all([x>rng(:,1); x<rng(:,2)]);

end