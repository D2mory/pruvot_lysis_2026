% script virusremoval_figS2.m
%
% Generates Fig. S2 of Pruvôt et al. (2026) using virusremoval_data.mat

%% Load data
load('virusremoval_data.mat')

%% Collect results

% total viral flux
x1 = arrayfun(@(s) s.with(5,5), flxs);
% viral flux normalised by total phytoplankton flux
x2 = arrayfun(@(s) s.with(5,5)./s.with(1,1), flxs);
% viral flux normalised by total DOC flux
x3 = arrayfun(@(s) s.with(5,5)./s.with(4,4), flxs);

% absolute change in phytoplankton flux
y1 = arrayfun(@(s) s.wout(1,1)-s.with(1,1), flxs);
% relative change in phytoplankton flux
y2 = arrayfun(@(s) (s.wout(1,1)-s.with(1,1)) ...
        /s.with(1,1), flxs);
% absolute change in DOC flux
y3 = arrayfun(@(s) s.wout(4,4)-s.with(4,4), flxs);
% relative change in DOC flux
y4 = arrayfun(@(s) (s.wout(4,4)-s.with(4,4)) ...
        /s.with(4,4), flxs);

% assemble variables for plotting
xx = [x1 x2 x3];
yy = [y1 y2 y3 y4];


%% Plot figure

col=[.9 .2 .2];

%x labels
labels={'F_V^{tot} ({\mu}mol-C/L/day)', ...
        'F_V^{tot}/F_P^{tot} (—)', ...
        'F_V^{tot}/F_D^{tot} (—)'};

for cc=1:3
    % Main histogram axis
    ah(cc,1)=subplot(2,4,cc);
    histogram(xx(:,cc),40,'normalization','probability', ...
        'facecolor',col,'edgecolor','w','facealpha',.8)
    xl=xlim; % store x-axis limits
    set(gca,'xtick',[],'fontsize',8)

    % Add y-label only to first panel
    if cc==1
        ylabel('Probability','fontsize',10)
    end

    % Shrink main axis vertically to make space for rug plot
    aux=get(gca,'position');
    set(gca,'position',[aux(1) aux(2)+aux(4)*.2 aux(3) aux(4)*.8])

    % Small lower axis (rug plot)
    ah(cc,2)=axes('position',aux.*[1 1 1 .15],'box','on');
    line([xx(:,cc) xx(:,cc)]', ...
        [zeros(size(xx(:,cc))) ones(size(xx(:,cc)))]', ...
        'color',brighten(col,.2))
    xlim(xl);
    set(gca,'fontsize',8)
    xlabel(labels{cc},'fontsize',12)
    set(gca,'ytick',[]);
end

labels={'\Delta{}F_P^{tot} ({\mu}mol-C/L/day)', ...
      '\delta{}F_P^{tot} (—)', ...
      '\Delta{}F_D^{tot} ({\mu}mol-C/L/day)', ...
      '\delta{}F_D^{tot} (—)'};

for cc=1:4
    bh(cc,1)=subplot(2,4,4+cc);
    histogram(yy(:,cc),40,'normalization','probability', ...
        'facecolor',col(1,:),'edgecolor','w','facealpha',.8)
    xl=xlim;
    set(gca,'xtick',[],'fontsize',8)
    if cc==1
        ylabel('Probability','fontsize',10)
    end
    aux=get(gca,'position');
    set(gca,'position',[aux(1) aux(2)+aux(4)*.2 aux(3) aux(4)*.8])

    bh(cc,2)=axes('position',aux.*[1 1 1 .15],'box','on');
    line([yy(:,cc) yy(:,cc)]', ...
        [zeros(size(yy(:,cc))) ones(size(yy(:,cc)))]', ...
        'color',brighten(col,.2))
    xlim(xl);
    set(gca,'fontsize',8)
    xlabel(labels{cc},'fontsize',12)
    set(gca,'ytick',[]);
    hold on
end
