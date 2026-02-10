function [F_mean, F] = WEI()
% function [F_mean, F] = WEI()
% Computes flux matrix for different the model of Weitz et al. 2015
%
% Inputs: none
%
% Outputs: 
% . F : flux matrices for all 1000 selected parameter sets
% . F_mean : mean flux matrix of the 1000 selected flux matrices

%% Fix parameter ranges (from Weitz et al. 2015)
pr=struct('muH',[.5 2], ...
          'Ko',[.25 1], ...
          'eH',[.05 .2], ...
          'muC',[.5 2], ...
          'KiC',[.05 1], ...
          'muE',[.2 2], ...
          'KiE',[.5 10], ...
          'phH',[1e-13 1e-10], ...
          'phC',[1e-13 1e-10], ...
          'phE',[1e-12 1e-10], ...
          'beH',[12.5 50], ...
          'beC',[12.5 100], ...
          'beE',[125 500], ...
          'mVH',[.05 5], ...
          'mVC',[.05 5], ...
          'mVE',[.05 5], ...
          'psH',[1e-6 1e-4], ...
          'psC',[1e-6 1e-4], ...
          'psE',[1e-6 1e-4], ...
          'pg',[.4 .4], ...
          'po',[.3 .3], ...
          'pi',[.3 .3], ...
          'mZ',[.025 .1], ...
          'mZP',[1e-8 1e-4], ...
          'omg',[.005 .02], ...
          'xsb',[2.5 10],...
          'qH',[5e-10 4e-9], ...
          'qC',[5e-10 4e-9], ...
          'qE',[5e-8 4e-7], ...
          'qZ',[5e-5 4e-4], ...
          'qV',[.5e-12 20e-12], ...
          'miH',[.001 .1], ...
          'miC',[.001 .1], ...
          'miE',[.001 .1], ...
          'moH',[.005 .1], ...
          'moC',[.005 .1], ...
          'moE',[.005 .1]);

%% Fix variable ranges (from Weitz et al. 2015)
%target
tar = [2e8 ... %H
       2e8 ... %C
       2e6 ... %E
       5e4 ... %Z
       2e9 ... %VH
       2e9 ... %VC
       2e7 ... %VE
       5 ... %xo
       .1 ... %xi
       10 ... viros-to-host ratio
       .5 ... %fractional mortality due to VH
       .25 ... %fractional mortality due to VC
       .1]; %fractional mortality due to VE

%% Find 10000 parameter sets compatible with non-negative population densities
nsim = 1e4; % desired number of parameter sets
ps(nsim) = pr; % structure to store accepted parameter sets

c1 = 1;
c2 = 1;
while c2 <= nsim
    pn = wei_smp(pr); % uniform sampling of parameter sets within the ranges
    xeq = wei_eqs(pn); % variables at equilibrium
    if all(xeq>0) % if all variables are non negative
        if mod(c2,10) == 0
            fprintf('Found parameter set %d at iteration %d\n', c2, c1)
        end
        ps(c2) = pn; % selected parameter set
        c2 = c2+1;
    end
    c1 = c1+1;
end

%% Selection of 1000 parameter sets closest to the desired variable ranges

% Out of the 10000 selected realisations, we calculate the distance to the
% targets (densities and fraction of viral mortality) 

dev = zeros(nsim,1);
for cc = 1:nsim 
    sol = wei_eqs(ps(cc));
    [~,ph] = wei_rhs(sol, ps(cc));

    aux = [ ...
        -(ph(1,5)+ph(1,10))/ph(1,1), ... %fraction of mortality to VH
        -(ph(2,6)+ph(2,10))/ph(2,2), ... % '' VC
        -(ph(3,7)+ph(3,10))/ph(3,3)]; % '' VE

    cmp = [sol(:)' sum(sol(5:7))/sum(sol(1:3)) aux]; %all criteria of selection
    dev(cc) = sum((log(cmp) - log(tar)).^2); %distance to the targets
end

% We then select the 1000 realisations the closest to the targets. 
nsel = 1000;
[~,idx] = mink(dev,nsel); % 1000 selected realisations

% Store selected solutions
res(nsel) = struct('pars',[],'sol',[]);
for k = 1:nsel
    i = idx(k);
    res(k).pars = ps(i);
    res(k).sol  = wei_eqs(ps(i));
end

%% Construct flux matrix and solve fluxes

% For the 1000 selected realisations closest to the targets, we extract the
% flux matrix from the original format of the Weitz et al. model (F_ori)
% and it is transformed into our common format: 5x5 (P,B,G,D,DV) 

F = struct('mat', cell(1,1000));   % 1×1000 struct
for c1 = 1:1e3
    pn = ps(c1); % parameter set values
    xeq = wei_eqs(pn); % equilibrium values
    [~,F_ori] = wei_rhs(xeq,pn); % original flux matrix from Weitz et al. 2015

    % Weitz et al. 2015: 10x10 (H,C,E,Z,VH,VC,VE,xo,xi,DV)
    % Our common format: 5x5 (P,B,G,D,DV)
    F_trf = zeros(5);

    % Phytoplankton
    F_trf(1,1) = F_ori(2,2)+F_ori(3,3); %primary production
    F_trf(1,3) = sum(F_ori(2:3,4)); %grazing
    F_trf(1,4) = sum(F_ori(2:3,8)); %DOC release
    F_trf(1,5) = sum(F_ori(2:3,10))+F_ori(2,6)+F_ori(3,7); %viral lysis

    % Heterotrophic bacteria
    F_trf(2,2) = F_ori(1,1); %bacterial growth
    F_trf(2,3) = F_ori(1,4); %grazing
    F_trf(2,4) = F_ori(1,8); %DOC release
    F_trf(2,5) = F_ori(1,5)+F_ori(1,10); %viral lysis

    % Grazers
    F_trf(3,3) = F_ori(4,4); %grazers growth (grazing on phytoplankton and bacteria)
    F_trf(3,4) = F_ori(4,8); % DOC release

    % DOC
    F_trf(4,2) = F_ori(8,1); %DOC bacterial consumption
    F_trf(4,4) = F_ori(8,8); %total of incoming DOC

    % DOC_V
    F_trf(5,4) = sum(F_trf(1:2,5)); %total of viral-induced DOC
    F_trf(5,5) = -sum(F_trf(1:2,5)); %inclusion in the total DOC
    
    F_trf = F_trf./F_trf(1,1)*100;

    F(c1).mat = F_trf;
end

%% Calculating the mean of all flux matrices

F_all = cat(3, F.mat); % 5×5×1000
F_mean = mean(F_all, 3); % 5×5

end