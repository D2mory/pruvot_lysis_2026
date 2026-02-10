% This code generate the data used for the 4th section of the article:
% "Common ground to compare modelling approaches"
%
% The data is save as a .mat file called "data_model_comparison.mat"

%% Flux matrices (flx) and mean flux matrices (fmn) for the 5 selected models
% Each are transformed into a common model format for comparison in their associated functions 

% Wilhelm & Suttle (1999)
[fmn_WIL, flx_WIL] = WIL();

% Weitz et al. (2015)
[fmn_WEI, flx_WEI] = WEI();

% Fuhrman (1992, 1999)
flx_FUH_B = FUH('B'); %viruses of bacteria only
flx_FUH_PB = FUH('PB'); %viruses of phytoplankton and bacteria 

% Mojica (2015)
flx_MOJ_S = MOJ('S'); %south
flx_MOJ_N = MOJ('N'); %north

% Xie et al. (2022)
flx_XIE_H = XIE('H'); %HOT sampling site
flx_XIE_A = XIE('A'); %AS sampling site

%%
refs = struct( ...
    'nms', {'WIL','WEI','FUH','MOJ','XIE'}, ... %codes
    'nml', { ... %full names
        'Wilhelm and Suttle (1999)', ...
        'Weitz et al. (2015)', ...
        'Fuhrman (1992, 1999)', ...
        'Mojica (2015)', ...
        'Xie et al. (2022)'}, ...
    'fmn', { ... %mean flux matrices
        fmn_WIL, ...
        fmn_WEI, ...
        [], ...
        [], ...
        [] });        
res(1:numel(refs)) = refs;

% All flux matrices
res(1).flx = flx_WIL;
res(2).flx = flx_WEI;
res(3).flx(1).mat = flx_FUH_B;
res(3).flx(2).mat = flx_FUH_PB;
res(4).flx(1).mat = flx_MOJ_S;
res(4).flx(2).mat = flx_MOJ_N;
res(5).flx(1).mat = flx_XIE_H;
res(5).flx(2).mat = flx_XIE_A;

%% Saving the data
save("data_model_comparison.mat", "res")