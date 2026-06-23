% script modelcomp_data.m
%
% Generates data for the comparison of previous model predictions
% on the role of marine viruses on carbon cycling;
% see Figs. 3 and S1 in Pruvôt et al. (2026)
%
% Data is saved in modelcomp_data.mat

%% Compute flux matrices

% for Wilhelm & Suttle (1999)
% flx_WIL with 1000 realisations, fmn_WIL with average flux matrix
[fmn_WIL, flx_WIL] = model_WIL_flx();

% for Weitz et al. (2015)
% flx_WEI with 1000 realisations, fmn_WEI with average flux matrix
[fmn_WEI, flx_WEI] = model_WEI_flx();

% for Fuhrman (1992, 1999)
flx_FUH_B = model_FUH_flx('B');    % model FUH-B
flx_FUH_PB = model_FUH_flx('BP');  % model FUH-BP

% for Mojica (2015)
flx_MOJ_S = model_MOJ_flx('S');    % model MOJ-S
flx_MOJ_N = model_MOJ_flx('N');    % model MOJ-N

% for Xie et al. (2022)
flx_XIE_H = model_XIE_flx('H');    % model XIE-H
flx_XIE_A = model_XIE_flx('A');    % model XIE-A

%% Store flux matrices in structure refs

res = struct( ...
        'nms', { ...                         % model identifier
           'WIL', 'WEI', 'FUH', 'MOJ', 'XIE'}, ...
        'nml', { ...                         % paper reference for model
           'Wilhelm and Suttle (1999)', ...
           'Weitz et al. (2015)', ...
           'Fuhrman (1992, 1999)', ...
           'Mojica (2015)', ...
           'Xie et al. (2022)'}, ...
        'fmn', { ...                         % mean flux matrix
           fmn_WIL, fmn_WEI, [], [], []}, ...
        'flx', { ...                         % individual flux matrices
           flx_WIL, flx_WEI, [], [], []});

res(3).flx(1).mat = flx_FUH_B;
res(3).flx(2).mat = flx_FUH_PB;
res(4).flx(1).mat = flx_MOJ_S;
res(4).flx(2).mat = flx_MOJ_N;
res(5).flx(1).mat = flx_XIE_H;
res(5).flx(2).mat = flx_XIE_A;

%% Save data
save('modelcomp_data.mat', 'res')
