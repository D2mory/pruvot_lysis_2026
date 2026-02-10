% analysis 1: Effects of virus removal using a simple microbial ecosystem model

%% Load data
load('data_virus_removal.mat'); %from 'data_genetation_virus_removal.m'

% Explanatory variables
x1 = arrayfun(@(s) s.with(5,5), flxs); % lysis flux
x2 = arrayfun(@(s) s.with(5,5) ./ s.with(1,1), flxs); % lysis flux / primary production

% Response variables
y1 = arrayfun(@(s) s.wout(1,1) - s.with(1,1), flxs); % Delta(primprod)
y2 = arrayfun(@(s) (s.wout(1,1) - s.with(1,1)) ./ s.with(1,1), flxs); % delta(primprod)

xx = [x1' x2'];
yy = [y1' y2'];

%% Data transformation
% apply gaussianisation

x1t=gauss_trf(x1);
x2t=gauss_trf(x2);
y1t=gauss_trf(y1);
y2t=gauss_trf(y2);

xxt=[gauss_trf(x1') gauss_trf(x2')];
yyt=[gauss_trf(y1') gauss_trf(y2')];

%% Figure 4
% determine ticks for x/y-axes
clf

fig = gcf;
set(fig,'Units','centimeters')
set(fig,'Position',[2 2 18 18])   % [left bottom width height]


tcks={[0 .15 1.2 8 30], ...                 % x1
      [0 .15 .8 1.1 1.2], ...               % x2
      [-10 0 5 50 500], ...                 % y1
      [-.9 0 3 30 150]};                    % y2

dats={[x1 x1t],[x2 x2t], ...
      [y1 y1t],[y2 y2t]};

% ticks for transformed variables
for cc=1:4
    [~,idx]=sort(dats{cc}(:,1));
    tcks{cc}(2,:)=0; % ticks
    for c1=1:size(tcks{cc},2)
        yy=interp1(dats{cc}(idx,1),dats{cc}(idx,2), ...
            tcks{cc}(1,c1),'linear','extrap');
        yy=max([yy -4]);
        yy=min([yy 4]);
        tcks{cc}(2,c1)=yy;
    end
end

tcls={{'0','0.15','1.2','8','30'}, ...      % x1
      {'0','0.15','0.8','1.1','1.2'}, ...   % x2
      {'-10','0','5','50','500'}, ...       % y1
      {'-0.9','0','3','30','150'}};         % y2

% determine labels for x/y-axis

xlabs={'F_V^{tot} (\mumol-C/L/day)','f_V^{tot} (-)'};
ylabs={'{\Delta}F_P^{tot} (\mumol-C/L/day)','{\delta}F_P^{tot} (-)'};

% generate figure

maxCounts = zeros(2,2);
for cx=1:2
    for cy=1:2
        subplot(2,2,cx+2*(cy-1))
        [counts, ~] = hexbin_plot(xxt(:,cx),yyt(:,cy),30,[-4 4 -4 4]);
        maxCounts(cy,cx) = max(counts);

        text(0.03, 0.98, sprintf('n_{max}=%d', maxCounts(cy,cx)), ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'fontsize',11);

        hold on

        axis([-4 4 -4 4])
        axis equal
        set(gca,'PositionConstraint','innerposition')

        yz=tcks{2+cy}(2,tcks{2+cy}(1,:)==0);
        plot([-4 4],yz*[1 1],'k--')
        set(gca,'xtick',tcks{cx}(2,:))
        set(gca,'XTickLabel',tcls{cx})
        set(gca,'ytick',tcks{2+cy}(2,:))
        set(gca,'YTickLabel',tcls{2+cy})
        set(gca,'fontsize',11)
    end
end

maxCounts % maximum hexagon density observed for each subplot

% variable labels

for cx=1:2
    subplot(2,2,cx+2)
    xlabel(xlabs{cx},'fontsize',14)
end
for cy=1:2
    subplot(2,2,1+2*(cy-1))
    ylabel(ylabs{cy},'fontsize',14)
end

% print figure

set(gcf,'color','w')
aux=get(gcf,'position');

cmap = colormap;

for k = [3 4]
    ax = subplot(2,2,k);
    pos = get(ax,'Position');
    pos(2) = pos(2) + 0.05;
    set(ax,'Position',pos);
end

fig = gcf;
% position of the legend bar
ax_leg = axes('Parent',fig,'Position',[0.4 0.07 0.25 0.02]);  % [left bottom width height]

% draw horizontal color strip

ncol = size(cmap,1);
for k = 1:ncol
    patch(ax_leg, [k k+1 k+1 k], [0 0 1 1], cmap(k,:), 'EdgeColor','none');
    hold(ax_leg,'on')
end

% set axis limits and appearance

xlim(ax_leg, [1 100])
ylim(ax_leg, [0 1])
set(ax_leg,'ytick',[],'Box','on','Layer','top')

% set ticks and label of the legend bar

aux=get(gca,'position');
set(gca,'position',[aux(1)-.01 aux(2)+aux(4)-.02 aux(3) .02])
xlim([-10 100])
set(gca,'xtick',[-10 1 50 100],'ytick',[])
set(gca,'xticklabel',{'0','1','n_{max}/2','n_{max}'}, 'FontSize', 12)
set(gca,'box','on','layer','top','ticklength',[.03 .03])
set(gca, 'XTickLabelRotation', 0)
xlabel('Number of model realisations', 'FontSize',13)



%% Save figure

exportgraphics(gcf, 'fig_4_remove1.png', 'Resolution',600)