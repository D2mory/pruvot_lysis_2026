% script virusremoval_fig4.m
%
% Generates Fig. 4 of Pruvôt et al. (2026) using virusremoval_data.mat

%% Load data
load('virusremoval_data.mat')

% Explanatory variables
x1 = arrayfun(@(s) s.with(5,5), flxs); % lysis flux
x2 = arrayfun(@(s) s.with(5,5) ./ s.with(1,1), flxs); % lysis flux / primary production

% Response variables
y1 = arrayfun(@(s) s.wout(1,1) - s.with(1,1), flxs); % Delta(primprod)
y2 = arrayfun(@(s) (s.wout(1,1) - s.with(1,1)) ./ s.with(1,1), flxs); % delta(primprod)

xx = [x1' x2'];
yy = [y1' y2'];

%% Transform data
% apply gaussianisation

x1t=gauss_trf(x1);
x2t=gauss_trf(x2);
y1t=gauss_trf(y1);
y2t=gauss_trf(y2);

xxt=[gauss_trf(x1') gauss_trf(x2')];
yyt=[gauss_trf(y1') gauss_trf(y2')];

%% Plot figure

% determine ticks for x/y-axes

tcks={[0 .15 1.2 8 30], ...                 % x1
      [0 .15 .8 1.1 1.2], ...               % x2
      [-10 0 5 50 500], ...                 % y1
      [-.9 0 3 30 150]};                    % y2
dats={[x1 x1t],[x2 x2t], ...
      [y1 y1t],[y2 y2t]};

for cc=1:4
    [~,idx]=sort(dats{cc}(:,1));
    tcks{cc}(2,:)=0;
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

% generate figure

maxCounts = zeros(2,2);
for cx=1:2
    for cy=1:2
        subplot(2,2,cx+2*(cy-1))
        [counts, ~] = hexbin_plot(xxt(:,cx),yyt(:,cy),30,[-4 4 -4 4]);
        maxCounts(cy,cx) = max(counts);
        hold on

        text(0.03, 0.98, sprintf('n_{max}=%d', maxCounts(cy,cx)), ...
            'Units', 'normalized', 'fontsize', 11, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
        axis([-4 4 -4 4])

        yz=tcks{2+cy}(2,tcks{2+cy}(1,:)==0);
        plot([-4 4],yz*[1 1],'k--')
        set(gca,'xtick',tcks{cx}(2,:))
        set(gca,'XTickLabel',tcls{cx})
        set(gca,'ytick',tcks{2+cy}(2,:))
        set(gca,'YTickLabel',tcls{2+cy})
        set(gca,'fontsize',11)
    end
end

% variable labels

xlabs={'F_V^{tot} (\mumol-C/L/day)','f_V^{tot} (-)'};
ylabs={'{\Delta}F_P^{tot} (\mumol-C/L/day)','{\delta}F_P^{tot} (-)'};

for cx=1:2
    subplot(2,2,cx+2)
    xlabel(xlabs{cx},'fontsize',14)
end
for cy=1:2
    subplot(2,2,1+2*(cy-1))
    ylabel(ylabs{cy},'fontsize',14)
end
