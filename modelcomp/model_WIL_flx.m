function [F_mean, F_real] = model_WIL_flx()
% function [F_mean, F_real] = model_WIL_flx()
% 
% Transforms the flux network reported in Wilhelm & Suttle (1999)
% into the format of flux matrices used in Pruvôt et al. (2026)
% To account for the variability reported in Wilhelm & Suttle (1999),
% we generate 1000 model realisations
%
% Inputs: none
%
% Outputs:
% . F_real: flux matrices for the 1000 model realisations
% . F_mean: mean flux matrix of the 1000 model realisations

%% Mean flux matrix

% Matrix reported in Wilhelm & Suttle (1999)
A = [.00 .00 .84 .00 .10 .06 .00;   % Phytoplankton
     .00 .00 .17 .00 .00 .25 .58;   % Heterotrophic bacteria
     .00 .00 .00 .31 .28 .01 .40;   % Grazers
     .00 .00 .00 .00 .18 .00 .82;   % Carnivores
     .00 .55 .00 .00 .00 .00 .45;   % DOC
     .00 .00 .00 .00 1.0 .00 .00;   % DOC_V (viral lysis)
     .00 .00 .00 .00 .00 .00 .00];  % Losses

% Compute flux matrix
c = [100 0 0 0 0 0 0]/(eye(7)-A);
P1 = diag(c)*(eye(7)-A);

% Reduce to 5x5 matrix
F_mean = P1([1:3 5:6],[1:3 5:6]);
F_mean(1:2,3) = sum(P1(1:2,3:4),2);
F_mean(3,4:5) = sum(P1(3:4,5:6));

%% Sampling within fluxes ranges

nReal = 1000;                          % number of realisations
F_real = struct('mat', cell(1,nReal)); % preallocate struct array

% Uncertainty for selected fluxes
a = [.08 .1 .26 .16 .15 .1];

for k = 1:nReal

    % Random perturbations in the range [-0.5, 0.5]
    rnd = rand(1,6) - .5;
    
    % Construct perturbed flux network
    A1 = A;

    % Phytoplankton
    A1(1,[3 6])   = A(1,[3 6])   + a(1)*rnd(1)*[1 -1];
    % Heterotrophic bacteria
    A1(2,[6 7])   = A(2,[6 7])   + a(2)*rnd(2)*[1 -1];
    % Grazers
    A1(3,[4 5 7]) = A(3,[4 5 7]) + a(3)*rnd(3)*[1 0 -1] ...
                                + a(4)*rnd(4)*[0 1 -1];
    % DOC
    A1(4,[5 7])   = A(4,[5 7])   + a(5)*rnd(5)*[1 -1];
    % DOC_V
    A1(5,[2 7])   = A(5,[2 7])   + a(6)*rnd(6)*[1 -1];
    
    % Compute flux matrix
    c = [100 0 0 0 0 0 0]/(eye(7)-A1);
    P1 = diag(c)*(eye(7)-A1);

    % Reduce to 5x5 matrix
    F_ind = P1([1:3 5:6],[1:3 5:6]);
    F_ind(1:2,3) = sum(P1(1:2,3:4),2);
    F_ind(3,4:5) = sum(P1(3:4,5:6));

    F_real(k).mat = F_ind;
end

end
