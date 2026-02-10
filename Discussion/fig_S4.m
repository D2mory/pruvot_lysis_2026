%%% FIGURE S4: Distribution of experimentally estimated viral lysis rates

% read excel files from Mojica & Brussaard
% Modified Dilution assay
T1=readtable('MDA_data.xlsx');
aux = T1.ViralLysisRate__d_1_;

if iscell(aux) || isstring(aux)
    aux = str2double(aux);
end

aux = aux(~isnan(aux));
dataMD = aux;

% Virus Reduction assay
T2 = readtable('VP_data.xlsx', ...
    'VariableNamingRule','preserve', ...
    'HeaderLines',1);
T2a = T2(:,[30 21 27]);   % VP, bacteria, burst size
T2a.Properties.VariableNames = ...
    {'VP_ml_1D_1_','Bacteria_ml_1_','BS_virusesCell_1_'};

aux = rmmissing(T2a);
aux = table2array(aux);
aux = aux(:,1)./aux(:,2)./aux(:,3);
dataVR = aux;


%% publication figure

% Subplot 1
col=[.9 .2 .2];
subplot(121) %%%
histogram(100*dataMD,0:5:160,'normalization','probability', ...
    'facecolor',col,'edgecolor','w','facealpha',.8)
xlabel('Process-based lysis percentage (%/day)')
ylabel('Probability')
ylim([0 .4])
set(gca,'ytick',0:.1:.4)
text(50,.3333,{'N=450 measurements', ...
             '25th percentile = 0.7 %/day', ...
             '75th percentile = 25 %/day'})
title('Modified dilution (MD) assays','fontweight','normal')

% Subplot 2
subplot(122) %%%
histogram(100*dataVR,0:20:680,'normalization','probability', ...
    'facecolor',col,'edgecolor','w','facealpha',.8)
xlabel('Process-based lysis percentage (%/day)')
ylabel('Probability')
ylim([0 .18])
set(gca,'ytick',0:.05:.15)
text(220,.15,{'N=86 measurements', ...
             '25th percentile = 35 %/day', ...
             '75th percentile = 140 %/day'})
title('Virus reduction (VR) assays','fontweight','normal')
% [length(dataMD) prctile(100*dataMD,[25 75])]
% [length(dataVR) prctile(100*dataVR,[25 75])]

aux=get(gcf,'position');
set(gcf,'position',[aux(1:2) 900 325])
set(gcf,'color','w')

%% Save figure

exportgraphics(gcf, 'figS4.png', 'Resolution',600)
