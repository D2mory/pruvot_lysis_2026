%%% FIGURE S2: distributions of variables from virus removal study

%% Load data

% Loads:
%   - flxs : struct array containing flux matrices
%           flxs(s).with  = fluxes with viruses
%           flxs(s).wout  = fluxes without viruses
%   - nrms : normalized residuals used to filter valid solutions

load fig45.mat

%%

% Keep only realizations with very small residuals (good convergence)
idxs = find(nrms<1e-6);

% Total viral flux (virus-derived DOM to itself)
x1 = arrayfun(@(s) flxs(s).with(5,5), idxs);

% Viral flux normalized by total phytoplankton flux
x2 = arrayfun(@(s) flxs(s).with(5,5)./flxs(s).with(1,1), idxs);

% Viral flux normalized by total detrital flux
x3 = arrayfun(@(s) flxs(s).with(5,5)./flxs(s).with(4,4), idxs);

% Absolute change in phytoplankton flux
y1 = arrayfun(@(s) flxs(s).wout(1,1)-flxs(s).with(1,1), idxs);

% Relative change in phytoplankton flux
y2 = arrayfun(@(s) (flxs(s).wout(1,1)-flxs(s).with(1,1)) ...
        /flxs(s).with(1,1), idxs);

% Absolute change in detrital flux
y3 = arrayfun(@(s) flxs(s).wout(4,4)-flxs(s).with(4,4), idxs);

% Relative change in detrital flux
y4 = arrayfun(@(s) (flxs(s).wout(4,4)-flxs(s).with(4,4)) ...
        /flxs(s).with(4,4), idxs);

% Assemble variables for plotting
xx = [x1' x2' x3'];
yy = [y1' y2' y3' y4'];


%% Plot distributions of viral flux metrics

col=[.9 .2 .2];

%x labels
tits={'F_V^{tot} ({\mu}mol-C/L/day)', ...
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
    xlabel(tits{cc},'fontsize',12)
    set(gca,'ytick',[]);
end
tits={'\Delta{}F_P^{tot} ({\mu}mol-C/L/day)', ...
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
    xlabel(tits{cc},'fontsize',12)
    set(gca,'ytick',[]);
    hold on
end

set(gcf,'color','w')
aux=get(gcf,'position');
set(gcf,'position',[aux(1:2) 920 500])

%% Save figure

exportgraphics(gcf, 'figS2.png', 'Resolution',600)
