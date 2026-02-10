function F_out = vr_flx_trf(F_ini)
% function F_out = vr_flx_trf(F_ini)
%
% Transform the matrix into the common format 

% 7x7 (P,B,VP,VB,G,D,DV)
% common: 5x5 (P,B,G,D,DV)
F_out=F_ini([1:2 5:7],[1:2 5:7]); %taking the initial matrix without VP and VB
F_out(1,5)=F_out(1,5)+F_ini(1,3); %VP flux into DV
F_out(2,5)=F_out(2,5)+F_ini(2,4); %VB flux into DV
F_out(5,4)=F_out(5,4)+sum(F_ini(3:4,6)); %fluxes exiting VP and VB, from DV to D
F_out(5,5)=F_out(5,5)+F_ini(3,3)+F_ini(4,4); %the through-flow flux of DV also include the trhough-flow flux of VP and VB

end