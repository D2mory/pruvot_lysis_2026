function x = vr_eqs(p)
% function x = vr_eqs(p)
% 
% This function calculates the steady state values for the different variables 
%
% Inputs:
% . p : parameter set as a structure
% Outputs:
% . x : vector of steady state values for the variables in order P,B,G,VP,VB,D

% Phytoplankton
P=p.mVP/(p.beP*p.phP);

% Heterotrophic bacteria
B=p.mVB/(p.beB*p.phB);

% Grazers
G=(p.pi*p.psP*P*p.qP/p.qG+p.pi*p.psB*B*p.qB/p.qG-p.r1G)/p.m2G;

% Phytopolankton viruses
VP=(p.muP*(1-P/p.KP)-p.psP*G-p.m1P-p.r1P)/p.phP;

% Heterotrophic bacteria viruses
VB=(p.m1P*P*p.qP+p.m1B*B*p.qB+p.m2G*G^2*p.qG ...
    +p.qP*p.phP*P*VP-p.qB/p.eB*p.psB*G*B ...
    +p.pe*(p.psP*P*G*p.qP+p.psB*B*G*p.qB) ...
    -p.qB/p.eB*(p.m1B+p.r1B)*B) ...
    /(p.qB*(1-p.eB)/p.eB*p.phB*B);

% DIssolved organic carbon (DOC)
D=p.KB*(p.phB*VB+p.psB*G+p.m1B+p.r1B) ...
   /(p.muB-(p.phB*VB+p.psB*G+p.m1B+p.r1B));

x=[P;B;VP;VB;G;D];

end