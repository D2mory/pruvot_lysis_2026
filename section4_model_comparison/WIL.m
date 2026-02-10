function [F_mean, F] = WIL()
% function F = WIL(mode)
% Computes flux matrix for the modle of Wilhelm and Suttle 1999
%
% Inputs: none
%
% Outputs: flux matrix F

%% Mean flux matrix
% Wilhelm & Suttle (P,B,G1,G2,D,VD,L)
A=[.00 .00 .84 .00 .10 .06 .00; % Phytoplankton
   .00 .00 .17 .00 .00 .25 .58; % Heterotrophic bacteria
   .00 .00 .00 .31 .28 .01 .40; % Grazers
   .00 .00 .00 .00 .18 .00 .82; % Carnivores
   .00 .55 .00 .00 .00 .00 .45; % DOC
   .00 .00 .00 .00 1.0 .00 .00; % DOC_V (viral lysis)
   .00 .00 .00 .00 .00 .00 .00]; % Losses

% Input vector: 100% of primary production
b=[100 0 0 0 0 0 0];

% Resolution
c=b/(eye(7)-A);
P1=diag(c)*(eye(7)-A);

% Reduce to our common 5x5 system
F_mean = P1([1:3 5:6],[1:3 5:6]);
F_mean(1:2,3) = sum(P1(1:2,3:4),2);
F_mean(3,4:5) = sum(P1(3:4,5:6));

%% Sampling within fluxes ranges
nReal = 1000; %number of realisations
F = struct('mat', cell(1,nReal));   % preallocate struct array

for k = 1:nReal
    % Wilhelm & Suttle (P,B,G1,G2,D,VD,L)
    A=[.00 .00 .84 .00 .10 .06 .00; % Phytoplankton
       .00 .00 .17 .00 .00 .25 .58; % Heterorophic bacteria
       .00 .00 .00 .31 .28 .01 .40; % Grazers
       .00 .00 .00 .00 .18 .00 .82; % Carnivores
       .00 .55 .00 .00 .00 .00 .45; % DOC
       .00 .00 .00 .00 1.0 .00 .00; % DOC_V (viral lysis)
       .00 .00 .00 .00 .00 .00 .00]; % Losses
    
    % Relative uncertainty amplitudes for selected fluxes
    a=[.08 .1 .26 .16 .15 .1];

    % Random perturbations in the range [-0.5, 0.5]
    rnd = rand(1,6) - .5;
    
    % Phytoplankton
    A(1,[3 6])   = A(1,[3 6])   + a(1)*rnd(1)*[1 -1];

    % Heterotrophic bacteria
    A(2,[6 7])   = A(2,[6 7])   + a(2)*rnd(2)*[1 -1];

    % Grazers
    A(3,[4 5 7]) = A(3,[4 5 7]) + a(3)*rnd(3)*[1 0 -1] ...
                                + a(4)*rnd(4)*[0 1 -1];

    % DOC
    A(4,[5 7])   = A(4,[5 7])   + a(5)*rnd(5)*[1 -1];

    % DOC_V
    A(5,[2 7])   = A(5,[2 7])   + a(6)*rnd(6)*[1 -1];
    
    % Input vector: 100% of primary production
    b = [100 0 0 0 0 0 0];

    % Resolution
    c = b/(eye(7)-A);
    P1 = diag(c)*(eye(7)-A);

    % Reduce to our common 5x5 system
    F_ind = P1([1:3 5:6],[1:3 5:6]);
    F_ind(1:2,3) = sum(P1(1:2,3:4),2);
    F_ind(3,4:5) = sum(P1(3:4,5:6));

    F(k).mat = F_ind;
end

end