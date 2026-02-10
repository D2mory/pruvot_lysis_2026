% analysis 2: compare fluxes without viruses, true vs reconstructed

%% Load data
load('data_virus_removal.mat'); %from 'data_genetation_virus_removal.m'

% create struct array with results
nSim = numel(flxs);

recs(length(flxs)) = struct('P1',[],'s1',[],'P2',[],'s2',[]);
for c1=1:nSim
    pn = psel(c1); % selected 8651 parameter sets
    P1 = flxs(c1).with; 
    P1 = diag(1./diag(P1))*P1; % carbon flux network with viruses
    s1 = [100 0 0 0 0 ]/P1;
    P1a = flxs(c1).wout; 
    P1a(5,4:5) = [-1 1];
    P1a = diag(1./diag(P1a))*P1a; % carbon flux network without viruses constructed from the dynamical model
    s1a = [100 0 0 0 0 ]/P1a;
    P2 = P1;
    P2(1,3:4) = P2(1,3:4)/(1+P2(1,5));
    P2(1,5) = 0;
    P2(2,3:4) = P2(2,3:4)/(1+P2(2,5));
    P2(2,5) = 0;
    % P1a : removal of viruses in the DS model, then transforamtion into the flux network format
    % P2 : removal of viruses in the FN model directly
    s2 = [100 0 0 0 0 ]/P2;
    recs(c1) = struct('P1',P1a,'s1',s1a,'P2',P2,'s2',s2);
end

%% Figure 5
% generate figure

clf

fig = gcf;
set(fig,'Units','centimeters')
set(fig,'Position',[2 2 18 18])   % [left bottom width height]

titleFS = 12; %title FontSize
maxCounts = zeros(3,3);

subplot(331) %%% f_Btot
x = arrayfun(@(s) s.s1(2), recs)/100;
y = arrayfun(@(s) s.s2(2), recs)/100;
[counts, ~] = hexbin_plot(x,y,30,[0 .25 0 .25]);
maxCounts(1,1) = max(counts);
hold on
plot([0 .25],[0 .25],'k--')
axis([0 .25 0 .25])
set(gca,'xtick',0:.05:.25,'ytick',0:.05:.25)
set(gca,'box','on')
title('Through-flow f_B^{tot}','FontWeight','normal', 'FontSize', titleFS)
ylabel('Flux in FN')

subplot(332) %%% f_Gtot
x = arrayfun(@(s) s.s1(3), recs)/100;
y = arrayfun(@(s) s.s2(3), recs)/100;
[counts, ~] = hexbin_plot(x,y,30,[0 .6 0 .6]);
maxCounts(1,2) = max(counts);
hold on
plot([0 .6],[0 .6],'k--')
axis([0 .6 0 .6])
set(gca,'xtick',0:.2:.6,'ytick',0:.2:.6)
set(gca,'box','on')
title('Through-flow f_G^{tot}','FontWeight','normal', 'FontSize', titleFS)

subplot(333) %%% f_Dtot
x = arrayfun(@(s) s.s1(4), recs)/100;
y = arrayfun(@(s) s.s2(4), recs)/100;
[counts, ~] = hexbin_plot(x,y,30,[0 1.2 0 1.2]);
maxCounts(1,3) = max(counts);
hold on
plot([0 1.2],[0 1.2],'k--')
axis([0 1.2 0 1.2])
set(gca,'xtick',0:.4:1.2,'ytick',0:.4:1.2)
set(gca,'box','on')
title('Through-flow f_D^{tot}','FontWeight','normal', 'FontSize', titleFS)

subplot(334) %%% a_PG
x = arrayfun(@(s) -s.P1(1,3), recs);
y = arrayfun(@(s) -s.P2(1,3), recs);
[counts, ~] = hexbin_plot(x,y,30,[0 .6 0 .6]);
maxCounts(2,1) = max(counts);
hold on
plot([0 .6],[0 .6],'k--')
axis([0 .6 0 .6])
set(gca,'xtick',0:.2:.6,'ytick',0:.2:.6)
set(gca,'box','on')
title('Flux partitioning a_{PG}','FontWeight','normal', 'FontSize', titleFS)
ylabel('Flux in FN')

subplot(335) %%% a_BG
x = arrayfun(@(s) -s.P1(2,3), recs);
y = arrayfun(@(s) -s.P2(2,3), recs);
[counts, ~] = hexbin_plot(x,y,30,[0 .6 0 .6]);
maxCounts(2,2) = max(counts);
hold on
plot([0 .6],[0 .6],'k--')
axis([0 .6 0 .6])
set(gca,'xtick',0:.2:.6,'ytick',0:.2:.6)
set(gca,'box','on')
title('Flux partitioning a_{BG}','FontWeight','normal', 'FontSize', titleFS)

subplot(336); %%% legend bar
col=[.9 .2 .2];
col(2,:)=.7+.3*col;
itp=linspace(1,0,100).^3;
cmap=ones(100,1)*col(1,:)+itp'*(col(2,:)-col(1,:));
for cc=1:100
    patch(cc+[0 1 1 0],[0 0 1 1],cmap(cc,:),'edgecolor','none')
    hold on
end
aux=get(gca,'position');
set(gca,'position',[aux(1)+0.04 aux(2)+aux(4)/2-.015 aux(3)-0.075 .03])
xlim([-10 100])
set(gca,'xtick',[-10 1 50 100],'ytick',[])
set(gca,'xticklabel',{'0','1','n_{max}/2','n_{max}'}, 'FontSize', 10)
set(gca, 'XTickLabelRotation', 0)
set(gca,'box','on','layer','top')
xlabel('Number of model realisations', 'FontSize', 12)

subplot(337) %%% a_PD
x = arrayfun(@(s) -s.P1(1,4), recs);
y = arrayfun(@(s) -s.P2(1,4), recs);
[counts, ~] = hexbin_plot(x,y,30,[0 1 0 1]);
maxCounts(3,1) = max(counts);
hold on
plot([0 1],[0 1],'k--')
axis([0 1 0 1])
set(gca,'xtick',0:.2:1,'ytick',0:.2:1)
set(gca,'box','on')
title('Flux partitioning a_{PD}','FontWeight','normal', 'FontSize', titleFS)
xlabel('Flux in DS')
ylabel('Flux in FN')

subplot(338) %%% a_BD
x = arrayfun(@(s) -s.P1(2,4), recs);
y = arrayfun(@(s) -s.P2(2,4), recs);
[counts, ~] = hexbin_plot(x,y,30,[0 1 0 1]);
maxCounts(3,2) = max(counts);
hold on
plot([0 1],[0 1],'k--')
axis([0 1 0 1])
set(gca,'xtick',0:.2:1,'ytick',0:.2:1)
set(gca,'box','on')
title('Flux partitioning a_{BD}','FontWeight','normal', 'FontSize', titleFS)
xlabel('Flux in DS')

subplot(339) %%% a_GD
x = arrayfun(@(s) -s.P1(3,4), recs);
y = arrayfun(@(s) -s.P2(3,4), recs);
[counts, ~] = hexbin_plot(x,y,30,[0 1 0 1]);
maxCounts(3,3) = max(counts);
hold on
plot([0 1],[0 1],'k--')
axis([0 1 0 1])
set(gca,'xtick',0:.2:1,'ytick',0:.2:1)
set(gca,'box','on')
title('Flux partitioning a_{GD}','FontWeight','normal', 'FontSize', titleFS)
xlabel('Flux in DS')

maxCounts % maximum hexagon density observed for each subplot

for cx=1:3
    for cy=1:3
        axis equal
        set(gca,'PositionConstraint','innerposition')
        if  ~(cy==2 && cx==3)
            subplot(3,3,(cy-1)*3 + cx)
            text(0.03, 0.98, sprintf('n_{max}=%d', maxCounts(cy, cx)), ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top', ...
            'fontsize',11);
        end 
    end
end

% print figure

set(gcf,'color','w')
aux=get(gcf,'position');
set(gcf,'position',[aux(1:2) 800 745])

%% Save figure

% exportgraphics(gcf, 'fig_5_remove2.png', 'Resolution',600)