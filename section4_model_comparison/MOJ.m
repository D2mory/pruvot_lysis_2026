function F = MOJ(region)
% function F = MOJ(region)
%
% MOJ adapts the model from Mojica et al. 2015 to a single format
%
% Inputs:
%   region : 'N' -> Northern region
%            'S' -> Southern region
%
% Outputs:
%   F : flux matrix

if nargin == 0
    error('MOJ requires an input: ''N'' or ''S''');
end

switch upper(region)

    case 'N'   % ---------- MOJ_N ----------
        A = [ ...
            47.2  0   -11.6  -16.3  -9.9;   % Phytoplankton
             0   10.4  -3.0   -1.7  -5.7;   % Heterotrophic bacteria
             0    0    14.6   -3.3   0;     % Grazers
             0  -10.4   0     36.9   0;     % DOC
             0    0     0    -15.6  15.6];  % DOC_V

    case 'S'   % ---------- MOJ_S ----------
        A = [ ...
            10.6  0   -2.7  -2.1  -3.7;   % Phytoplankton
             0   9.0  -3.8  -0.4  -4.8;   % Heterotrophic bacteria
             0    0   6.5  -1.4   0;      % Grazers
             0  -9.0   0   12.4   0;      % DOC
             0    0    0   -8.5   8.5];   % DOC_V

    otherwise
        error('Unknown MOJ region. Use ''N'' or ''S''.');
end

%% Flux matrix (normalized to phytoplankton production)
F = A / A(1,1) * 100;

end
