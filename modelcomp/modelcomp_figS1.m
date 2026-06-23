% script modelcomp_figS1.m
%
% Generates Fig. S1 of Pruvôt et al. (2026) using modelcomp_data.mat
% This data file is provided, or can be generated using modelcomp_data.m

%% Load data

load modelcomp_data.mat

%% Extract data matrix

% each of the 5 datasets is weighted equally
% . for datasets WIL and WEI, we take 1000 realisations
% . for datasets FUH, MOJ and XIE, we take 500 copies of each model
%   variant, so that there are also 1000 realisations in total

% X_data1 contains relevant fluxes without equal weighting
X_data1 = zeros(0,10);
for cc=1:5
    x1 = zeros(length(res(cc).flx),10);
    x1(:,1) = arrayfun(@(s) s.mat(3,3)/s.mat(1,1), res(cc).flx);    % fGtot
    x1(:,2) = arrayfun(@(s) s.mat(4,4)/s.mat(1,1), res(cc).flx);    % fDtot
    x1(:,3) = arrayfun(@(s) s.mat(5,5)/s.mat(1,1), res(cc).flx);    % fVtot
    x1(:,4) = arrayfun(@(s) -s.mat(1,3)/s.mat(1,1), res(cc).flx);   % aPG
    x1(:,5) = arrayfun(@(s) -s.mat(1,4)/s.mat(1,1), res(cc).flx);   % aPD^0
    x1(:,6) = arrayfun(@(s) -s.mat(1,5)/s.mat(1,1), res(cc).flx);   % aPD^V
    x1(:,7) = arrayfun(@(s) -s.mat(2,3)/s.mat(4,4), res(cc).flx);   % aBG
    x1(:,8) = arrayfun(@(s) -s.mat(2,4)/s.mat(4,4), res(cc).flx);   % aBD^0
    x1(:,9) = arrayfun(@(s) -s.mat(2,5)/s.mat(4,4), res(cc).flx);   % aBD^V
    x1(:,10) = arrayfun(@(s) -s.mat(3,4)/s.mat(3,3), res(cc).flx);  % aGD
    X_data1 = [X_data1; x1];
end

rngs = [1 1000; 1001 2000; 2001 2002; 2003 2004; 2005 2006];

%% Plot figure

cols = [244 162 97; 42 157 143; 231 111 81; 233 196 106;0 109 119]/255;

% variable labels
lab_flx = {'f_G^{tot}', 'f_D^{tot}', 'f_V^{tot}', ...
          'a_{PG}', 'a_{PD}^{0}', 'a_{PD}^{V}', ...
          'a_{BG}', 'a_{BD}^0', 'a_{BD}^V', 'a_{GD}'};
% model labels
lab_mod = {'WIL','WEI','FUH-B','FUH-BP','MOJ-S','MOJ-N','XIE-H','XIE-A'};

% limits and ticks on x-axes
xlims={[0 1],[0 1.25],[0 1.25],[0 1],[0 .5],[0 1],[0 .5],[0 .5],[0 .5],[0 .5]}; 
xtcks={0:.2:1,0:.25:1.25,0:.25:1.25,0:.2:1,0:.1:.5,0:.2:1,0:.1:.5,0:.1:.5,0:.1:.5,0:.1:.5};

splts=[1 2 3 5 6 7 9 10 11 12];
for cc = 1:10 
    subplot(3,4,splts(cc))

    % for datasets WIL and WEI draw histograms
    for m = 1:2
        data = X_data1(rngs(m,1):rngs(m,2), cc);
        xlm = xlims{cc};
        bins = linspace(xlm(1),xlm(2),21);
        hh(m) = histogram(data,'BinEdges',bins, ...
            'Normalization','probability', ...
            'EdgeColor','k','FaceColor',cols(m,:),'FaceAlpha',0.5);
        hold on
    end
    
    % for datasets FUH, MOJ and XIE draw vertical lines
    for m = 3:5
        data_points1 = X_data1(rngs(m,1), cc);
        data_points2 = X_data1(rngs(m,2), cc);
        hh(2*m-3) = plot(data_points1*[1 1], [0 .8], 'Color', cols(m,:), ...
            'LineStyle','-', 'LineWidth', 2);
        hh(2*m-2) = plot(data_points2*[1 1], [0 .8], 'Color', cols(m,:), ...
            'LineStyle',':', 'LineWidth', 2);
    end
    
    title(lab_flx{cc},'fontweight','normal')
    axis([xlims{cc} 0 .8])
    set(gca,'xtick',xtcks{cc})
    set(gca,'ytick',0:.2:.8)
    set(gca, 'Box','on', 'FontSize', 12)
end

% construct legend
hl=legend(hh,lab_mod);
set(hl,'Location','north','Box','off')
aux=get(hl,'position');
set(hl,'position',aux+[0 .59 0 0])
