%% Figure SI: distribution of the different fluxes, for the studied models

load('data_model_comparison.mat');

%% Data extraction

% Construct data matrix - equal weighting
X_data1 = zeros(0,10);
for cc=1:5
    x1 = zeros(length(res(cc).flx),10);
    x1(:,1) = arrayfun(@(s) s.mat(3,3), res(cc).flx); %fGtot
    x1(:,2) = arrayfun(@(s) s.mat(4,4), res(cc).flx); %fDtot
    x1(:,3) = arrayfun(@(s) s.mat(5,5), res(cc).flx); %fVtot
    x1(:,4) = arrayfun(@(s) -s.mat(1,3)/s.mat(1,1), res(cc).flx); %aPG
    x1(:,5) = arrayfun(@(s) -s.mat(1,4)/s.mat(1,1), res(cc).flx); %a0PD
    x1(:,6) = arrayfun(@(s) -s.mat(1,5)/s.mat(1,1), res(cc).flx); %aVPD
    x1(:,7) = arrayfun(@(s) -s.mat(2,3)/s.mat(4,4), res(cc).flx); %aDBG
    x1(:,8) = arrayfun(@(s) -s.mat(2,4)/s.mat(4,4), res(cc).flx); %a0DBD
    x1(:,9) = arrayfun(@(s) -s.mat(2,5)/s.mat(4,4), res(cc).flx); %aVDBD    
    x1(:,10) = arrayfun(@(s) -s.mat(3,4)/s.mat(3,3), res(cc).flx); %aGD
    X_data1 = [X_data1; x1];
end

X_data2 = zeros(5000,10);
X_data2(1:1000,:) = X_data1(1:1000,:);
X_data2(1001:2000,:) = X_data1(1001:2000,:);
X_data2(2001:3000,:) = repmat(X_data1(2001:2002,:), 500, 1);
X_data2(3001:4000,:) = repmat(X_data1(2003:2004,:), 500, 1);
X_data2(4001:5000,:) = repmat(X_data1(2005:2006,:), 500, 1);

X_data = X_data2;
rngs = [
    1     1000;   
    1001  2000;   
    2001  3000;   
    3001  4000;   
    4001  5000];  

%% Plotting figure S1

% Colors
cols = [
    244 162 97;  
    42 157 143;  
    0 109 119;   
    233 196 106; 
    231 111 81;  
] / 255;

% Variables labels
labflx = {'f_G^{tot}', 'f_D^{tot}', 'f_V^{tot}', ...
          'a_{PG}', 'a_{PD}^{0}', 'a_{PD}^{V}', ...
          'a_{BG}', 'a_{BD}^0', 'a_{BD}^V', ...
          'a_{GD}'};

% Models labels
modelNames = {'WIL', 'WEI', {'FUH-B', 'FUH-PB'}, {'MOJ-S','MOJ-N'}, {'XIE-H', 'XIE-A'}};

% X-axis limits
xlim_f = [0 130];
xlim_a1 = [0 1];
xlim_a05 = [0 .5];

% Number of bins for the histograms
nBins = 60;

% BinEdges
binEdges_f = linspace(xlim_f(1), xlim_f(2), nBins + 1);
binEdges_a1 = linspace(xlim_a1(1), xlim_a1(2), nBins + 1);
binEdges_a05 = linspace(xlim_a05(1), xlim_a05(2), nBins + 1);

% Compute bin edges per variable
nVar = size(X_data,2);


%%% Plot %%%
clf
for cc = 1:size(X_data,2)
    sp_idx = cc;
    if cc == 10
        sp_idx = 11; % move a_GD to 11
    end

    ax(cc) = subplot(4,3,sp_idx);
    hold on; grid on
    
    % Adapting x limits and ticks depending on the variable
    n_xIntervals = 5;
    if cc <= 3
        binEdges = binEdges_f;
        xlim(xlim_f)
        xticks_vals = xlim_f(1):20:xlim_f(2); 
        xticks(xticks_vals)
    elseif cc==4 || cc==6 %a_PG or a_VPD
        binEdges = binEdges_a1;
        xlim(xlim_a1)
        xticks_vals = linspace(xlim_a1(1), xlim_a1(2), n_xIntervals + 1);
        xticks(xticks_vals)
    else 
        binEdges = binEdges_a05;
        xlim(xlim_a05)
        xticks_vals = linspace(xlim_a05(1), xlim_a05(2), n_xIntervals + 1);
        xticks(xticks_vals)
    end

    % WIL & WEI
    yMax_hist = 0;
    iden_lines = gobjects(0);
    for m = 1:2
        data_m = X_data(rngs(m,1):rngs(m,2), cc);
    
        if all(data_m == data_m(1))   % case with identical values (no incertainties)
            % draw temporary line (fixed later)
            p = plot([data_m(1) data_m(1)], [0 1], 'Color', cols(m,:), 'LineWidth', 2);
            iden_lines(end+1) = p;  % store handle
            hCur = p;
        else
            counts = histcounts(data_m, binEdges);
            yMax_hist = max(yMax_hist, max(counts));
    
            h = histogram(data_m, ...
                'BinEdges', binEdges, ...
                'FaceColor', cols(m,:), ...
                'FaceAlpha', 0.5, ...
                'EdgeColor', 'none');
            hCur = h;
        end
    
        if cc == 1
            hLegObj(m) = hCur;
        end
    end
    yMax = yMax_hist;
    
    % Order of magnitude
    order = 10^floor(log10(yMax));
    yMax = ceil(yMax / order) * order;
    
    % 4 intervals
    n_yIntervals = 4;
    yticks_vals = linspace(0, yMax, n_yIntervals + 1);
    ylim([0 yMax])
    yticks(yticks_vals)

    for k = 1:numel(iden_lines)
        iden_lines(k).YData = [0 yMax];
    end

    % FUH, MOJ, XIE → xlines
    for m = 3:5
        data_points1 = X_data(rngs(m,1), cc);
        data_points2 = X_data(rngs(m,2), cc);

        l1 = xline(data_points1, 'Color', cols(m,:), 'LineStyle','-', 'LineWidth', 2);
        l2 = xline(data_points2, 'Color', cols(m,:), 'LineStyle',':', 'LineWidth', 2);

        % store one handle per model for legend (do it once)
        if cc == 1
            hLegObj(m) = l1; % or l2, choice is arbitrary, just color/style
        end
    end

    title(labflx{cc})
    set(gca, 'Box','on','XGrid','on','YGrid','on', ...
        'XLimMode', 'manual', ...
        'FontSize', 12)
end


%%% Legend in a single block %%%
% Create handles
h_FUH_B  = plot(NaN, NaN, 'Color', cols(3,:), 'LineStyle','-',  'LineWidth',2);
h_FUH_PB = plot(NaN, NaN, 'Color', cols(3,:), 'LineStyle',':',  'LineWidth',2);
h_MOJ_S  = plot(NaN, NaN, 'Color', cols(4,:), 'LineStyle','-',  'LineWidth',2);
h_MOJ_N  = plot(NaN, NaN, 'Color', cols(4,:), 'LineStyle',':',  'LineWidth',2);
h_XIE_H  = plot(NaN, NaN, 'Color', cols(5,:), 'LineStyle','-',  'LineWidth',2);
h_XIE_A  = plot(NaN, NaN, 'Color', cols(5,:), 'LineStyle',':',  'LineWidth',2);

% Collect handles and labels
hLegObj_legend = [hLegObj(1), hLegObj(2), h_FUH_B, h_FUH_PB, h_MOJ_S, h_MOJ_N, h_XIE_H, h_XIE_A];
legStrings = {'WIL', 'WEI', 'FUH-B', 'FUH-PB', 'MOJ-S', 'MOJ-N', 'XIE-H', 'XIE-A'};

% Create legend subplot
ax_leg = subplot(4,3,10);
axis(ax_leg, 'off');
lgd = legend(ax_leg, hLegObj_legend, legStrings, 'Location', 'north');
lgd.NumColumns = 2;
lgd.FontSize = 12;
lgd.Position(3:4) = [0.08 0.12]; % width height of subplot
lgd.Position(1:2) = [0.2 0.12]; % reposition to fit: x, y
lgd.Units = 'normalized';

%% Save figure

exportgraphics(gcf, 'figS1.png', 'Resolution',600)