function x = wei_eqs(p)
% function x = wei_eqs(p)
%
% Equilibrium expression from Weitz et al. 2015
%
% Inputs:
% . p : parameters
%
% Outputs:
% . x : system's equilibra in order H, C, E, Z, VH, VC, VE, xo, xi

% Heterotrophs
H=p.mVH/(p.beH*p.phH);

% Cyanobacteria
C=p.mVC/(p.beC*p.phC);

% Eukaryotic autotrophs
E=p.mVE/(p.beE*p.phE);

% Zooplankton
Z=(p.pg/p.qZ*(p.qH*p.psH*H+p.qC*p.psC*C+p.qE*p.psE*E)-p.mZ)/p.mZP;

% Inorganic nitrogen
xi=p.xsb-p.qZ*p.mZP*Z^2/p.omg;

% Viruses of cyanobacteria
VC=(p.muC*xi/(xi+p.KiC)-p.psC*Z-p.miC-p.moC)/p.phC;

% Viruses of eukaryotic autotrophs
VE=(p.muE*xi/(xi+p.KiE)-p.psE*Z-p.miE-p.moE)/p.phE;

% Viruses of heterotrophs
NN=-p.qH*H/p.eH*(p.psH*Z+p.miH+p.moH)+p.qV*p.mVC*VC+p.qV*p.mVE*VE ...
   +(p.qC-p.qV*p.beC)*p.phC*C*VC+(p.qE-p.qV*p.beE)*p.phE*E*VE ...
   +p.po*p.qH*p.psH*H*Z+p.po*p.qC*p.psC*C*Z+p.po*p.qE*p.psE*E*Z ...
   +p.qH*p.moH*H+p.qC*p.moC*C+p.qE*p.moE*E;

VH=NN/(p.qH*p.phH*H/p.eH-p.qV*p.mVH-(p.qH-p.qV*p.beH)*p.phH*H);

% Organic nitrogen
xo=p.Ko*(p.phH*VH+p.psH*Z+p.miH+p.moH)/(p.muH-p.phH*VH-p.psH*Z-p.miH-p.moH);


x=[H;C;E;Z;VH;VC;VE;xo;xi];
end