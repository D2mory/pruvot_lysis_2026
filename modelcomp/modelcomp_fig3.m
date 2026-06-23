% script modelcomp_fig3.m
%
% Generates Fig. 3 of Pruvôt et al. (2026) using modelcomp_data.mat
% This data file is provided, or can be generated using modelcomp_data.m

%% Load data

load modelcomp_data.mat

%% Construct data matrix

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
    x1(:,5) = arrayfun(@(s) -s.mat(2,3)/s.mat(4,4), res(cc).flx);   % aBG
    x1(:,6) = arrayfun(@(s) -s.mat(1,4)/s.mat(1,1), res(cc).flx);   % aPD^0
    x1(:,7) = arrayfun(@(s) -s.mat(2,4)/s.mat(4,4), res(cc).flx);   % aBD^0
    x1(:,8) = arrayfun(@(s) -s.mat(1,5)/s.mat(1,1), res(cc).flx);   % aPD^V
    x1(:,9) = arrayfun(@(s) -s.mat(2,5)/s.mat(4,4), res(cc).flx);   % aBD^V
    x1(:,10) = arrayfun(@(s) -s.mat(3,4)/s.mat(3,3), res(cc).flx);  % aGD
    X_data1 = [X_data1; x1];
end

% X_data contains relevant fluxes with equal weighting
X_data = zeros(5000,10);
X_data(1:1000,:) = X_data1(1:1000,:);
X_data(1001:2000,:) = X_data1(1001:2000,:);
X_data(2001:3000,:) = repmat(X_data1(2001:2002,:), 500, 1);
X_data(3001:4000,:) = repmat(X_data1(2003:2004,:), 500, 1);
X_data(4001:5000,:) = repmat(X_data1(2005:2006,:), 500, 1);

rngs = [1 1000; 1001 2000; 2001 3000; 3001 4000; 4001 5000]; 

%% Construct PCA

X_mean = mean(X_data, 1);
X_ctrd = X_data - X_mean;
[U, S, V] = svd(X_ctrd, 'econ');
scores = X_ctrd * V(:, 1:2);
loadings = V(:, 1:2);

% Enforce orientation of PCA axes
% large positive PC1 should correspond to strong viral lysis of phytoplankton
if dot(scores(:,1), X_ctrd(:,8)) < 0
    V(:,1)        = -V(:,1);
    scores(:,1)   = -scores(:,1);
    loadings(:,1) = -loadings(:,1);
end
% large positive PC2 should correspond to strong viral lysis of bacteria
if dot(scores(:,2), X_ctrd(:,9)) < 0
    V(:,2)        = -V(:,2);
    scores(:,2)   = -scores(:,2);
    loadings(:,2) = -loadings(:,2);
end

% Compute graphical representation
eigvals = diag(S).^2;
explained = eigvals / sum(eigvals) * 100;
n = size(X_ctrd,1);
lambda = eigvals / (n-1);

idx = [1 2 3 4 8 9];                        % selected variables
      % 1:fGtot 2:fDtot 3:fVtot 4:aPG 8:aPD^V 9:aBD^V
not_idx = setdiff(1:size(V,1), idx);        % non-selected variables

coords = V(:,1:2) .* sqrt(lambda(1:2)');    % projections on PC1 and PC2
arrow_lengths = sqrt(sum(coords.^2, 2));
missing_lengths = arrow_lengths(not_idx);
max_missing_length = max(missing_lengths);  % longest non-selected variable


%%  Plot first subplot

cols = [244 162 97; 42 157 143; 231 111 81; 233 196 106; 0 109 119]/255;
set(gcf,'Units','centimeters')
set(gcf,'Position',[2 2 50 18])  

axes('Position',[0.08 0.15 0.55 0.75]); 
hold on

% model WIL
iModel = 1;
x = scores(rngs(iModel,1):rngs(iModel,2),1);
y = scores(rngs(iModel,1):rngs(iModel,2),2);
% determine KDE grid
xrng = [min(x) max(x)];
xrng = xrng + .2*diff(xrng)*[-1 1];
xlin = linspace(xrng(1), xrng(2), 100);
yrng = [min(y) max(y)];
yrng = yrng + .2*diff(yrng)*[-1 1];
ylin = linspace(yrng(1), yrng(2), 100);
[Xgrid, Ygrid] = meshgrid(xlin, ylin);
% compute KDE function
ww = .02;
zz = zeros(size(Xgrid));
for ctr = 1:length(x)
    zz = zz+exp(-(Xgrid-x(ctr)).^2/2/ww^2 ...
                -(Ygrid-y(ctr)).^2/2/ww^2);
end
% determine mode of KDE function
[~, idxMax] = max(zz(:));
x_mode = Xgrid(idxMax);
y_mode = Ygrid(idxMax);
% plot density lines
contour(Xgrid, Ygrid, zz, 5, ...
    'LineColor', cols(iModel,:), 'LineWidth', 1.5);
scatter(x_mode, y_mode, 150, cols(iModel,:), 'x', ...
    'LineWidth', 3, 'MarkerEdgeAlpha', 1);

% model WEI
iModel = 2;
x = scores(rngs(iModel,1):rngs(iModel,2),1);
y = scores(rngs(iModel,1):rngs(iModel,2),2);        
% determine KDE grid
xrng = [min(x) max(x)];
xrng = xrng + .2*diff(xrng)*[-1 1];
xlin = linspace(xrng(1), xrng(2), 100);
yrng = [min(y) max(y)];
yrng = yrng + .2*diff(yrng)*[-1 1];
ylin = linspace(yrng(1), yrng(2), 100);
[Xgrid, Ygrid] = meshgrid(xlin, ylin);
% compute KDE function
ww = .03;
zz = zeros(size(Xgrid));
for ctr = 1:length(x)
    zz = zz+exp(-(Xgrid-x(ctr)).^2/2/ww^2 ...
                -(Ygrid-y(ctr)).^2/2/ww^2);
end
% determine mode of KDE function
[~, idxMax] = max(zz(:));
x_mode = Xgrid(idxMax);
y_mode = Ygrid(idxMax);
% plot density lines
contour(Xgrid, Ygrid, zz, 5, ...
    'LineColor', cols(iModel,:), 'LineWidth', 1.5);
scatter(x_mode, y_mode, 150, cols(iModel,:), 'filled');

% model FUH
iModel = 3;
x = scores(rngs(iModel,1):rngs(iModel,2),1);
y = scores(rngs(iModel,1):rngs(iModel,2),2);
scatter(x, y, 150, cols(iModel,:), 'x', ...
    'MarkerEdgeColor', cols(iModel,:), ...
    'MarkerFaceAlpha', 0.1, 'LineWidth', 3);

% mode MOJ
iModel = 4;
x = scores(rngs(iModel,1):rngs(iModel,2),1);
y = scores(rngs(iModel,1):rngs(iModel,2),2);
scatter(x, y, 150, cols(iModel,:), 'x', ...
    'MarkerEdgeColor', cols(iModel,:), ...
    'MarkerFaceAlpha', 0.1, 'LineWidth', 3);

% model XIE
iModel = 5;
x = scores(rngs(iModel,1):rngs(iModel,2),1);
y = scores(rngs(iModel,1):rngs(iModel,2),2);
scatter(x, y, 150, cols(iModel,:), 'filled');

% variables: arrows and labels
labs = {'f^{tot}_G','f^{tot}_D','f^{tot}_V', ...
    'a_{PG}','a_{BG}', 'a^{0}_{PD}','a^{0}_{BD}', ...
    'a^{V}_{PD}','a^{V}_{BD}', 'a_{GD}'};
for i = idx
    quiver(0, 0, loadings(i,1), loadings(i,2), ...
        'Color',[.3 .3 .3], 'LineWidth',1.5);
    if i == 1
        text(loadings(i,1), loadings(i,2)+0.05, labs{i}, ...
            'Color',[.3 .3 .3], 'FontSize',16);
    elseif i == 4
        text(loadings(i,1), loadings(i,2)-0.05, labs{i}, ...
            'Color',[.3 .3 .3], 'FontSize',16);
    elseif i == 7
        text(loadings(i,1), loadings(i,2)-0.05, labs{i}, ...
            'Color',[.3 .3 .3], 'FontSize',16);
    else
        text(loadings(i,1), loadings(i,2), labs{i}, ...
            'Color',[.3 .3 .3], 'FontSize',16);
    end
end

% model labels
% model WIL
text(scores(rngs(1,1),1)-0.2, scores(rngs(1,1),2)+0.07, 'WIL', ...
    'Color', cols(1,:), 'FontSize', 13, 'FontWeight', 'bold');
% model WEI
text(scores(rngs(2,1),1)-0.1, scores(rngs(2,1),2)+0.14, 'WEI', ...
    'Color', cols(2,:), 'FontSize', 13, 'FontWeight', 'bold');
% model FUH
text(scores(rngs(3,1),1)-0.2, scores(rngs(3,1),2)+0.01, 'FUH-B', ...
    'Color', cols(3,:), 'FontSize', 13, 'FontWeight', 'bold');
text(scores(rngs(3,2),1)-0.2, scores(rngs(3,2),2)+0.06, 'FUH-BP', ...
    'Color', cols(3,:), 'FontSize', 13, 'FontWeight', 'bold');
% model MOJ
text(scores(rngs(4,1),1)+0.03, scores(rngs(4,1),2)+0.05, 'MOJ-S', ...
    'Color', cols(4,:), 'FontSize', 13, 'FontWeight', 'bold');
text(scores(rngs(4,2),1)-0.1, scores(rngs(4,2),2)-0.07, 'MOJ-N', ...
    'Color', cols(4,:), 'FontSize', 13, 'FontWeight', 'bold');
% model XIE
text(scores(rngs(5,1),1)-0.12, scores(rngs(5,1),2)-0.05, 'XIE-H', ...
    'Color', cols(5,:), 'FontSize', 13, 'FontWeight', 'bold');
text(scores(rngs(5,2),1)-0.13, scores(rngs(5,2),2)-0.05, 'XIE-A', ...
    'Color', cols(5,:), 'FontSize', 13, 'FontWeight', 'bold');

% plot axes
plot([-1 2],[0 0], ':', 'Color','k', 'LineWidth', 1)
plot([0 0],[-1 1], ':', 'Color','k', 'LineWidth', 1)

% plot circle of non-selected variables
theta = linspace(0, 2*pi, 200);
r = max_missing_length;  % radius = longest arrow length
x_circ = r * cos(theta);
y_circ = r * sin(theta);
fill(r * cos(theta), r * sin(theta), [0.8 0.8 0.8], ...
        'FaceAlpha', 0.4, 'EdgeColor', 'none');
plot(x_circ, y_circ, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);

xlim([-1 1.5])
ylim([-0.5 0.75])
xticks(-1:0.5:1.5)
yticks(-1.5:0.5:1.5)

xlabel(sprintf('PC1 (%.0f%%)', explained(1)));
ylabel(sprintf('PC2 (%.0f%%)', explained(2)));
set(gca,'Box','on','XGrid','on','YGrid','on','FontSize',14)

text(0.025,0.94,'A','Units','normalized', ...
    'FontSize',24,'FontWeight','bold')

%% Plot second subplot

% extract viral lysis fluxes
lysis_flx = 100*X_data1(:,3);

axes('Position',[0.70 0.15 0.25 0.75]);  
hold on

% model WIL
m=1; x=2;
data = lysis_flx(1:1000);
[f, xi] = ksdensity(data);
f = f / max(f) * 0.4;
fill([x - f, fliplr(x + f)], [xi, fliplr(xi)], cols(m,:), 'FaceAlpha',0.5,'EdgeColor','none');
% find mode for 1D KDE
[~, idxMax] = max(f);
mode_y = xi(idxMax);
if m==1
    scatter(x, mode_y, 150, cols(m,:), 'x', 'LineWidth', 3, 'MarkerEdgeAlpha', 1);
elseif m==2
    scatter(x, mode_y, 150, cols(m,:), 'filled');
end

% model WEI
m=2; x=5;
data = lysis_flx(1001:2000);
[f, xi] = ksdensity(data);
f = f / max(f) * 0.4;
fill([x - f, fliplr(x + f)], [xi, fliplr(xi)], cols(m,:), 'FaceAlpha',0.5,'EdgeColor','none');
% find mode for 1D KDE
[~, idxMax] = max(f);
mode_y = xi(idxMax);
if m==1
    scatter(x, mode_y, 150, cols(m,:), 'x', 'LineWidth', 3, 'MarkerEdgeAlpha', 1);
elseif m==2
    scatter(x, mode_y, 150, cols(m,:), 'filled');
end

% model FUH
m=3; x=3;
scatter(x,lysis_flx(2001),150,'x','LineWidth',3, ...
    'MarkerEdgeColor',cols(3,:),'MarkerFaceAlpha',0.1);
scatter(x,lysis_flx(2002),150,'x','LineWidth',3, ...
    'MarkerEdgeColor',cols(3,:),'MarkerFaceAlpha',0.1);

% model MOJ
m=4; x=4;
scatter(x,lysis_flx(2003),150,'x','LineWidth',3, ...
    'MarkerEdgeColor',cols(4,:),'MarkerFaceAlpha',0.1);
scatter(x,lysis_flx(2004),150,'x','LineWidth',3, ...
    'MarkerEdgeColor',cols(4,:),'MarkerFaceAlpha',0.1);

% model XIE
m=5; x=1;
scatter(x,lysis_flx(2005),150,cols(5,:),'filled');
scatter(x,lysis_flx(2006),150,cols(5,:),'filled');

% model labels
text(1-0.3, lysis_flx(2005)+9, 'XIE-H', 'Color', cols(5,:), ...
    'FontSize', 13, 'FontWeight', 'bold')
text(1-0.3, lysis_flx(2006)+5, 'XIE-A', 'Color', cols(5,:), ...
    'FontSize', 13, 'FontWeight', 'bold')
text(3-0.3, lysis_flx(2001)-5, 'FUH-B', 'Color', cols(3,:), ...
    'FontSize', 13, 'FontWeight', 'bold')
text(3-0.3, lysis_flx(2002)+5, 'FUH-BP', 'Color', cols(3,:), ...
    'FontSize', 13, 'FontWeight', 'bold')
text(4-0.3, lysis_flx(2003)+5, 'MOJ-S', 'Color', cols(4,:), ...
    'FontSize', 13, 'FontWeight', 'bold')
text(4-0.3, lysis_flx(2004)+5, 'MOJ-N', 'Color', cols(4,:), ...
    'FontSize', 13, 'FontWeight', 'bold')

xticks(1:5)
xticklabels({'XIE','WIL','FUH','MOJ','WEI'})
xlim([0.5 5.5])
ylim([0 125])

ylabel('Lysis percentage $\mathcal{L}^{\mathrm{flux}}$', ...
    'Interpreter','latex','fontsize',18)
set(gca,'FontSize',14,'Box','on','XGrid','on','YGrid','on')

text(0.05, 0.94, 'B', 'Units','normalized', ...
    'FontSize', 24, 'FontWeight','bold')
