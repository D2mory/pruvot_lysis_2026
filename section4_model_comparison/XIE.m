function F = XIE(site)
% function F = XIE(site)
%
% XIE adapts the model from Xie et al. (2022) to a single format
%
% Inputs:
%   site : 'A' -> Arabian Sea (1995)
%          'H' -> HOT site (2002)
%
% Outputs:
%   F : flux matrix

if nargin == 0
    error('XIE requires an input: ''A'' or ''H''');
end

%% Original flux matrix ratio
A = zeros(9);

switch upper(site)

    case 'A'   % ---------- XIE_A (Arabian Sea) ----------

        %1 PHY
        A(1, 3:9) = [ ...
            0.551 ...   % PRT
            0 ...       % MZ
            0 ...       % Virus
            0.143 ...   % LDOC (modified for equilibrium, original value =0.043 (+0.1))
            0.069 ...   % SDOC
            0.221 ...   % DET (modified for equilibrium, original value =0.121 (+0.1))
            0];         % EXT

        %2 BA
        A(2, 3:9) = [ ...
            0.168 ...   % PRT
            0 ...       % MZ
            0.009 ...   % Virus
            0.175 ...   % LDOC
            0 ...       % SDOC
            0 ...       % DET
            0.371 + 0.005]; % EXT

        %3 PRT
        A(3, 4:9) = [ ...
            0.283 ...   % MZ
            0 ...       % Virus
            0.088 ...   % LDOC
            0.040 ...   % SDOC
            0.020 ...   % DET
            0.271];     % EXT

        %4 MZ
        A(4, 6:9) = [ ...
            0.051 ...   % LDOC
            0.020 ...   % SDOC
            0.044 ...   % DET
            0.168];     % EXT

        %5 Virus
        A(5,6) = 0.009;

        %6 LDOC
        A(6,2) = 0.472;

        %7 SDOC
        A(7,2) = 0.242;

        %8 DET
        A(8,7) = 0.117;
        A(8,9) = 0.162;

        %9 EXT
        A(9,1) = 1.0;

    case 'H'   % ---------- XIE_H (HOT) ----------

        %1 PHY
        A(1, 3:9) = [ ...
            0.511 + 0.178 ... % PRT
            0.111 ...         % MZ
            0 ...             % Virus
            0.082 + 0.034 + 0.011 ... % LDOC
            0.041 + 0.004 + 0.002 ... % SDOC
            0.007 + 0.002 + 0.017 ... % DET
            0];               % EXT

        %2 BA
        A(2, 3:9) = [ ...
            0.128 ...
            0 ...
            0.009 ...
            0.115 ...
            0 ...
            0 ...
            0.002 + 0.306];

        %3 PRT
        A(3, 4:9) = [ ...
            0.013 ...
            0 ...
            0.134 ...
            0.029 ...
            0.022 ...
            0.619];

        %4 MZ
        A(4, 6:9) = [ ...
            0.028 ...
            0.013 ...
            0.023 ...
            0.060];

        %5 Virus
        A(5,6) = 0.009;

        %6 LDOC
        A(6,2) = 0.402;

        %7 SDOC
        A(7,2) = 0.146;

        %8 DET
        A(8,7) = 0.028;
        A(8,9) = 0.043;

        %9 EXT
        A(9,1) = 1.0;

    otherwise
        error('Unknown XIE site. Use ''A'' or ''H''.');
end

%% Partitioning ratios matrix (identical for both)
PR = @(A)[ ...
    0 0 sum(A(1,3:4)) ...
      sum(A(1,6:7)) + A(1,8)*(A(8,7)/(A(8,7)+A(8,9))) ...
      A(1,5); ...
    0 0 sum(A(2,3:4))/sum(A(:,2)) ...
      sum(A(2,6:7))/sum(A(:,2)) ...
      A(2,5)/sum(A(:,2)); ...
    0 0 0 ...
      (sum(A(3:4,6:7),'all')) / ...
      (sum(A([1:2,5:end],3)) + sum(A([1:2,5:end],4))) + ...
      (A(3,8)+A(4,8))*(A(8,7)/(A(8,7)+A(8,9))) ...
      0; ...
    0 round(sum(A(6:7,2)) / ...
      (sum(A([1:4 6:end],6)) + sum(A([1:4 6:end],7)))) ...
      0 0 0; ...
    0 0 0 1 0 ];

PR_mat = PR(A);

%% Resolution
solveF = @(PR) diag([100 0 0 0 0] / (eye(5) - PR)) * (eye(5) - PR);

F = solveF(PR_mat);

end
