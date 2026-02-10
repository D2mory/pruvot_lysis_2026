function F = FUH(host)
% function F = FUH(host)
%
% FUH adapts Fuhrman (1992 / 1999) models to a single format
%
% Inputs:
%   type : 'B'  -> Fuhrman 1992 (bacteria only)
%          'PB' -> Fuhrman 1999 (phyto + bacteria)
%
% Outputs:
%   F : flux matrix

if nargin == 0
    error('FUH requires an input: ''B'' or ''PB''');
end

%% Original flux matrix ratio
A = zeros(9);

switch upper(host)

    case 'B'   % ---------- FUH_B ----------
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

    case 'PB'  % ---------- FUH_PB ----------
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

    otherwise
        error('Unknown FUH type. Use ''B'' or ''PB''.');
end

A = A / 100;

%% Partitioning ratios matrix
PR = @(A)[ ...
    0 0 sum(A(1,3:5)) A(1,7) A(1,6); ...
    0 0 sum(A(2,3:5))/sum(A(:,2)) ...
          sum(A(2,7))/sum(A(:,2)) ...
          A(2,6)/sum(A(:,2)); ...
    0 0 0 ...
      sum(A(3:5,7),'all') / ...
      (sum(A([1:2,6:end],3)) + ...
       sum(A([1:2,6:end],4)) + ...
       sum(A([1:2,6:end],5))) ...
      0; ...
    0 A(7,2)/sum(A(:,7)) 0 0 0; ...
    0 0 0 1 0 ];

PR_mat = PR(A);

%% Resolution
solveF = @(PR) diag([100 0 0 0 0] / (eye(5) - PR)) * (eye(5) - PR);

F = solveF(PR_mat);

end
