function [dxc,phc] = virusremoval_model_rhsC(xc,p)
% function [dxc,phc] = virusremoval_model_rhsC(xc,p)
%
% Evaluates the right-hand side of the dynamical equations of
% the minixmal marine microbial ecosystem model
% All variables are expressed in carbon units, which makes
% the integration of the differential equations more stable
%
% Inputs:
% . xc : system state (vector of length 6)
% . p : parameters (stucture with 26 fields)
%
% Outputs:
% . dxc : right-hand side of dynamical equations
% . phc : right-hand side written as flux matrix
%
% Notes:
% . Variables and fluxes are expressed in carbon units
%   parameters are expressed in original units

xc = xc(:);
convC = [p.qP p.qB p.qG p.qV p.qV 1]';
x = xc./convC;
[dx,ph] = virusremoval_model_rhs(x,p);
dxc = convC.*dx;
phc = diag(convC)*ph;

end
