function x = virusremoval_model_eqs(p)
% function x = virusremoval_model_eqs(p)
%
% Evaluates the equilibrium state of the minimal marine microbial ecosystem
%
% Inputs:
% . p : parameters (structure with 26 fields)
%
% Outputs:
% . x : equilibrium state
%
% Notes:
% . The order of the variables in the state vector is P, B, G, VP, VB, D

% P and B
P = p.lVP/(p.betaP*p.phiP);
B = p.lVB/(p.betaB*p.phiB);

% G
G = 1/p.lG_cons * (p.pG_assim* ...
    (p.psiPG*P*p.qP/p.qG + p.psiBG*B*p.qB/p.qG) ...
    - p.lG_inorg);

% total losses
ellP = p.lP_inorg + p.lP_org;
ellB = p.lB_inorg + p.lB_org;

% VP
VP = 1/p.phiP * (p.muP_max*(1-P/p.KP) - p.psiPG*G - ellP);

% auxiliary ratios
cPB = p.qP*P/(p.qB*B);
cGB = p.qG*G/(p.qB*B);

% x1
x1 = cPB*(p.phiP*VP + p.pG_org*p.psiPG*G + p.lP_org ) ...
    + p.pG_org*p.psiBG*G + p.lB_org + cGB*p.lG_cons*G;

% VB
VB = (p.pDB_assim*x1 - p.psiBG*G - ellB) / ...
     (p.phiB*(1-p.pDB_assim));

% x2
x2 = p.phiB*VB + p.psiBG*G + ellB;

% D
D = p.KDB*x2/(p.muB_max - x2);

x = [P;B;G;VP;VB;D];

end
