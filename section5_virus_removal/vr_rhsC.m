function [dxc,phc] = vr_rhsC(xc,pc)
% function [dxc,phc] = vr_rhsC(xc,pc)
%
% Dynamical equations in carbon units
%
% Inputs:
% . xc : variables steady state values in carbon
% . pc : parameters in carbon
% Outputs: 
% . dxc : variables dynamical equations
% . phc : flux matrix (dimension 6x8)
%        . column 7: virus-driven DOC
%        . column 8: loss flux

Pc=xc(1);
Bc=xc(2);
VPc=xc(3);
VBc=xc(4);
Gc=xc(5);
Dc=xc(6);

phc = zeros(7,8);
% Phytoplankton (P) dynamics
phc(1,:) = [pc.muPc*Pc*(1-Pc/pc.KPc) 0 ...
           -pc.phPc*Pc*VPc*pc.bePc 0 ...
           -pc.pic*pc.psPc*Pc*Gc -pc.m1Pc*Pc-pc.pec*pc.psPc*Pc*Gc ...
           -pc.phPc*Pc*VPc*(1-pc.bePc) ...
           -pc.r1Pc*Pc-(1-pc.pic-pc.pec)*pc.psPc*Pc*Gc ];
dPc = sum(phc(1,:));

% Heterotrophic bacteria (B) dynamics 
phc(2,:) = [0 pc.muBc*Dc*Bc/(pc.KBc+Dc) ...
           0 -pc.phBc*Bc*VBc*pc.beBc ...
           -pc.pic*pc.psBc*Bc*Gc -pc.m1Bc*Bc-pc.pec*pc.psBc*Bc*Gc ...
           -pc.phBc*Bc*VBc*(1-pc.beBc) ...
           -pc.r1Bc*Bc-(1-pc.pic-pc.pec)*pc.psBc*Bc*Gc];
dBc = sum(phc(2,:));

% Phytoplakton viruses (VP) dynamics
phc(3,:) = [0 0 pc.bePc*pc.phPc*Pc*VPc 0 0 -pc.mVPc*VPc 0 0];
dVPc= sum(phc(3,:));


% Heterotrophic bacteria viruses (VB) dynamics
phc(4,:) = [0 0 0 pc.beBc*pc.phBc*Bc*VBc 0 -pc.mVBc*VBc 0 0];
dVBc= sum(phc(4,:));

% Grazers (G) dynamics
phc(5,:) = [0 0 0 0 pc.pic*pc.psPc*Gc*Pc + pc.pic*pc.psBc*Gc*Bc ...
           -pc.m2Gc*Gc^2 0 -pc.r1Gc*Gc];
dGc = sum(phc(5,:));

% Dissolved organic carbon (DOC) dynamics
phc(6,:) = [0 -pc.muBc*Dc*Bc/(pc.KBc+Dc) 0 0 0 ...
         pc.m1Pc*Pc + pc.m1Bc*Bc + pc.m2Gc*Gc^2 + ...
         pc.mVPc*VPc + pc.mVBc*VBc + ...
         (1-pc.bePc)*pc.phPc*Pc*VPc + (1-pc.beBc)*pc.phBc*Bc*VBc + ...
         pc.pec*pc.psPc*Pc*Gc + pc.pec*pc.psBc*Bc*Gc 0 ...
         -(1/pc.eBc-1)*pc.muBc*Dc*Bc/(pc.KBc+Dc)];
dDc = sum(phc(6,:));

% Fluxes for virus-driven DOC
phc(7,:) = [0 0 0 0 0 -1 1 0] * ...
    (pc.phPc*Pc*VPc*(1-pc.bePc)+pc.phBc*Bc*VBc*(1-pc.beBc));

dxc=[dPc;dBc;dVPc;dVBc;dGc;dDc];

end