function F_out = virusremoval_model_std(F_ini)
% function F_out = virusremoval_model_std(F_ini)
%
% Transforms the flux matrix obtained from the minimal marine microbial
% ecosystem into the standard format used by Pruvôt et al (2026)
%
% Inputs:
% . F_ini: flux matrix from ecosystem model
%
% Outputs:
% . F_out: flux matrix in standard format

% recall order of variables
% in F_ini: P,B,G,VP,VB,D,DV
% in F_out: P,B,G,D,DV

F_out = zeros(5);
F_out(1:4,:) = F_ini([1:3 6],[1:3 6:7]);
F_out(1,5) = F_out(1,5)+F_ini(1,4);
F_out(2,5) = F_out(2,5)+F_ini(2,5);
F_out(5,4) = sum(F_out(1:2,5));
F_out(5,5) = -sum(F_out(1:2,5));

end