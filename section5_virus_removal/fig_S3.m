% analysis 1 (Supplementary Information): Effects of virus removal using a simple microbial ecosystem model

%% Load data
load('data_virus_removal.mat'); % from 'data_genetation_virus_removal.m'

% Explanatory variables
x1 = arrayfun(@(s) s.with(5,5), flxs); % lysis rate
x2 = arrayfun(@(s) s.with(5,5) ./ s.with(1,1), flxs); % lysis flux / primary production
x3 = arrayfun(@(s) s.with(5,5) ./ s.with(4,4), flxs); % lysis flux / DOC production

% Response variables
y1 = arrayfun(@(s) s.wout(1,1) - s.with(1,1), flxs); % Delta(primprod)
y2 = arrayfun(@(s) (s.wout(1,1) - s.with(1,1)) ./ s.with(1,1), flxs); % delta(primprod)
y3 = arrayfun(@(s) s.wout(4,4) - s.with(4,4), flxs); % Delta(DOCprod)
y4 = arrayfun(@(s) (s.wout(4,4) - s.with(4,4)) ./ s.with(4,4), flxs); % delta(DOCprod)

xx = [x1' x2' x3'];
yy = [y1' y2' y3' y4'];

%% Data transformation
% apply gaussianisation

x1t=gauss_trf(x1);
x2t=gauss_trf(x2);
x3t=gauss_trf(x3);
y1t=gauss_trf(y1);
y2t=gauss_trf(y2);
y3t=gauss_trf(y3);
y4t=gauss_trf(y4);

xxt=[gauss_trf(x1') gauss_trf(x2') gauss_trf(x3')];
yyt=[gauss_trf(y1') gauss_trf(y2') gauss_trf(y3') gauss_trf(y4')];

%% Figure S2
% determine ticks for x/y-axes
clf 

fig = gcf;
set(fig,'Units','centimeters')
set(fig,'Position',[2 2 18 18])   % [left bottom width height]

tcks={[0 .15 1.2 8 30], ...                 % x1
      [0 .15 .8 1.1 1.2], ...               % x2
      [0 .2 .8 .97 .99], ...                % x3
      [-10 0 5 50 500], ...                 % y1
      [-.9 0 3 30 150], ...                 % y2
      [-20 -5 0 10 300], ...                % y3
      [-.9 0 3 20 100]};                    % y4

dats={[x1 x1t],[x2 x2t],[x3 x3t], ...
      [y1 y1t],[y2 y2t],[y3 y3t],[y4 y4t]};

for cc=1:7
    [~,idx]=sort(dats{cc}(:,1));
    tcks{cc}(2,:)=0;    % ticks for transformed variables
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
      {'0','0.2','0.8','0.97','0.99'}, ...  % x3
      {'-10','0','5','50','500'}, ...       % y1
      {'-0.9','0','3','30','150'}, ...      % y2
      {'-20','-5','0','10','300'}, ...      % y3
      {'-0.9','0','3','20','100'}};         % y4

% determine labels for x/y-axis

xlabs={'F_V^{tot}','F_V^{tot}/F_P^{tot}','F_V^{tot}/F_D^{tot}'};
ylabs={'{\Delta}F_P^{tot}','{\delta}F_P^{tot}', ...
       '{\Delta}F_D^{tot}','{\delta}F_D^{tot}'};

% generate figure

maxCounts = zeros(4,3);
for cx=1:3
    for cy=1:4
        subplot(4,3,cx+3*(cy-1))
        [counts, ~] = hexbin_plot(xxt(:,cx),yyt(:,cy),30,[-4 4 -4 4]);
        maxCounts(cy,cx) = max(counts);
        hold on
        yz=tcks{3+cy}(2,tcks{3+cy}(1,:)==0);
        plot([-4 4],yz*[1 1],'k--')         % draw line y=0
        set(gca,'xtick',tcks{cx}(2,:))      % set ticks x-axis
        set(gca,'XTickLabel',tcls{cx})      % set ticklabels x-axis
        set(gca,'ytick',tcks{3+cy}(2,:))    % set ticks y-axis
        set(gca,'YTickLabel',tcls{3+cy})    % set ticklabels y-axis
        set(gca,'fontsize',9)
        text(0.03, 0.98, sprintf('n_{max}=%d', maxCounts(cy, cx)), ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top', ...
            'fontsize',11);
        axis square
    end
end
maxCounts % maximum hexagon density observed for each subplot

for cx=1:3
    subplot(4,3,cx+9)
    xlabel(xlabs{cx},'fontsize',14)
end
for cy=1:4
    subplot(4,3,1+3*(cy-1))
    ylabel(ylabs{cy},'fontsize',14)
end

% print figure

set(gcf,'color','w')
aux=get(gcf,'position');
set(gcf,'position',[aux(1:2) 750 960])

cmap = colormap;

% shifting slightly the subplots

for k = 10:12
    ax = subplot(4,3,k);
    pos = get(ax,'Position'); 
    pos(2) = pos(2) + 0.04;
    set(ax,'Position',pos);
end
for k = 4:9 
    ax = subplot(4,3,k);
    pos = get(ax,'Position');
    pos(2) = pos(2) + 0.02;  
    set(ax,'Position',pos);
end

fig = gcf;
ax_leg = axes('Parent',fig,'Position',[0.4 0.06 0.23 0.02]); % [left bottom width height]

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

% set ticks

aux=get(gca,'position');
set(gca,'position',[aux(1)-.01 aux(2)+aux(4)-.02 aux(3) .02])
xlim([-10 100])
set(gca,'xtick',[-10 1 50 100],'ytick',[])
set(gca,'xticklabel',{'0','1','n_{max}/2','n_{max}'}, 'FontSize', 12)
set(gca,'box','on','layer','top','ticklength',[.03 .03])
set(gca, 'XTickLabelRotation', 0)
xlabel('Number of model realisations', 'FontSize',13)

%% Save figure

exportgraphics(gcf, 'figS3.png', 'Resolution',600)