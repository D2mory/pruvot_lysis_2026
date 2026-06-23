function chk = virusremoval_model_chk(x,p,rng)
% function chk = virusremoval_model_chk(x,p,rng)
%
% Checks whether variables are within specified ranges
% and whether parameters are consistent
%
% Inputs:
% . x : variables as state vector (vector of length 6)
% . p : parameters (structure with 26 fields)
% . rng : ranges of variable values (matrix with dimensions 6x2)
%
% Outputs : 
% . chk : true if all checks passed

% check whether variables are within specified ranges
x = x(:);
chk1 = all([x>rng(:,1); x<rng(:,2)]);

% check whether parameters are consistent
chk2 = (p.qP>p.betaP*p.qV) && (p.qB>p.betaB*p.qV) && (p.pG_assim+p.pG_org<1);

chk = chk1 && chk2;

end
