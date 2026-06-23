function F = model_MOJ_flx(variant)
% function F = model_MOJ_flx(variant)
%
% Transforms the flux networks reported in Mojica et al. (2015)
% into the format of flux matrices used in Pruvôt et al. (2026)
%
% Inputs:
%   variant: specifies model variant
%   . variant=='N': Northern region
%   . variant=='S': Southern region
%
% Outputs:
%   F : flux matrix

variant = char(variant);
switch variant

    case 'N'   % ---------- model MOJ-N ----------
        A = [ ...
            47.2  0   -11.6  -16.3  -9.9;   % Phytoplankton
             0   10.4  -3.0   -1.7  -5.7;   % Heterotrophic bacteria
             0    0    14.6   -3.3   0;     % Grazers
             0  -10.4   0     36.9   0;     % DOC
             0    0     0    -15.6  15.6];  % DOC_V

    case 'S'   % ---------- model MOJ-S ----------
        A = [ ...
            10.6  0   -2.7  -2.1  -3.7;   % Phytoplankton
             0   9.0  -3.8  -0.4  -4.8;   % Heterotrophic bacteria
             0    0   6.5  -1.4   0;      % Grazers
             0  -9.0   0   12.4   0;      % DOC
             0    0    0   -8.5   8.5];   % DOC_V
end

% Normalise relative to primary production
F = A / A(1,1) * 100;

end
