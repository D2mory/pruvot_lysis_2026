function F = model_FUH_flx(variant)
% function F = model_FUH_flx(variant)
%
% Transforms the flux networks reported in Fuhrman (1992, 1999)
% into the format of flux matrices used in Pruvôt et al. (2026)
%
% Inputs:
%   variant: specifies model variant
%   . variant=='B': flux network from Fuhrman (1992);
%     includes only viruses infecting heterotrophic bacteria
%   . variant=='BP': flux network from Fuhrman (1999);
%     includes both viruses infecting phytoplankton
%     and viruses infecting bacteria
%
% Outputs:
%   F : flux matrix

%% Original flux network
A = zeros(9);

variant = char(variant);
switch variant

    case 'B'   % ---------- model FUH-B ----------
        %1 Phytoplankton
        A(1, 3:9) = [20 20 30 0 30 0 0];
        %2 Bacteria
        A(2, 3:9) = [19 0 0 19 0 38 0];
        %3 Nanozooplankton
        A(3, 4:9) = [15 0 0 8 16 0];
        %4 Microzooplankton
        A(4, 5:9) = [13 0 9 13 0];
        %5 Macrozooplankton
        A(5, 7:9) = [10 19 14];
        %6 Viruses
        A(6, 7) = 19;
        %7 DOM
        A(7, 2) = 76;
        %8 CO2
        A(8, 1) = 100;
        %9 EXT
        A(9, 8) = 14;

    case 'BP'  % ---------- model FUH-BP ----------
        %1 Phytoplankton
        A(1, 3:9) = [18 18 27 3 34 0 0];
        %2 Bacteria
        A(2, 3:9) = [20 0 0 20 0 40 0];
        %3 Nanozooplankton
        A(3, 4:9) = [16 0 0 9 16 0];
        %4 Microzooplankton
        A(4, 5:9) = [14 0 8 12 0];
        %5 Macrozooplankton
        A(5, 7:9) = [9 18 14];
        %6 Viruses
        A(6, 7) = 23;
        %7 DOM
        A(7, 2) = 80;
        %8 CO2
        A(8, 1) = 100;
        %9 EXT
        A(9, 8) = 14;
end

A = A / 100;

%% Compute partitioning ratios
PR = [ 0 0 sum(A(1,3:5)) A(1,7) A(1,6); ...
       0 0 sum(A(2,3:5))/sum(A(:,2)) ...
           sum(A(2,7))/sum(A(:,2)) ...
           A(2,6)/sum(A(:,2)); ...
       0 0 0 sum(A(3:5,7)) / ...
               (sum(A([1:2,6:end],3)) + ...
                sum(A([1:2,6:end],4)) + ...
                sum(A([1:2,6:end],5))) 0; ...
       0 A(7,2)/sum(A(:,7)) 0 0 0; ...
       0 0 0 1 0 ];

%% Compute flux matrix
F = diag([100 0 0 0 0] / (eye(5) - PR)) * (eye(5) - PR);

end
