function [F_mean, F_real] = model_WEI_flx()
% function [F_mean, F_real] = model_WEI_flx()
%
% Generates flux matrixes in the format used in Pruvôt et al. (2026)
% associated to equilibrium states of the marine microbial ecosystem model
% described in Weitz et al. (2015)
% To account for the variability of the model of Weitz et al. (2015),
% we generate 1000 model realisations
%
% Inputs: none
%
% Outputs: 
% . F_real: flux matrices for the 1000 model realisations
% . F_mean: mean flux matrix of the 1000 model realisations

%% Set parameter ranges
% parameter ranges are taken from Weitz et al. (2015)

pr = struct('muH',[.5 2], ...
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

%% Set target values for dynamical variables
% target values are taken from Weitz et al. (2015)

tar = [2e8 ...  % abundance of heterotrophs
       2e8 ...  % abundance of cyanobacteria
       2e6 ...  % abundance of eukaryotic autotrophs
       5e4 ...  % abundance of zooplankton
       2e9 ...  % abundance of viruses of heterotrophs
       2e9 ...  % abundance of viruses of cyanobacteria
       2e7 ...  % abundance of viruses of eukaryotic autotrophs
       5 ...    % concentration of organic nitrogen
       .1 ...   % concentration of inorganic nitrogen
       10 ...   % virus-to-host abundance ratio
       .5 ...   % mortality fraction due to viruses of heterotrophs
       .25 ...  % mortality fraction due to viruses of cyanobacteria
       .1];     % mortality fraction due to virues of eukaryotic autotrophs

%% Construct 10000 feasible realisations
% randomly sample parameter values from parameter ranges and require
% that all dynamical variables have positive equilibrium values
% here target values for dynamical variables are not yet considered

nSim = 1e4;     % desired number of parameter sets
ps(nSim) = pr;  % structure to store accepted parameter sets

fprintf('Starting simulations of Weitz et al. (2015)\n');

c1 = 1;
c2 = 1;
while c2 <= nSim
    pn = model_WEI_smp(pr);   % sample parameter values uniformly from ranges
    xeq = model_WEI_eqs(pn);  % compute dynamical variables at equilibrium
    if all(xeq>0)             % all variables should be positive
        if mod(c2,100) == 0
            fprintf('. found parameter set %d (out of 10000) at simulation %d\n', c2, c1)
        end
        ps(c2) = pn;          % store accepted parameter set
        c2 = c2+1;
    end
    c1 = c1+1;
end

%% Select 1000 realisations closest to target values

% calculate distance to targets for 10000 simulated realisations
dev = zeros(nSim,1);
for cc = 1:nSim 
    sol = model_WEI_eqs(ps(cc));
    [~,ph] = model_WEI_rhs(sol, ps(cc));
    aux = [ ...
      -(ph(1,5)+ph(1,10))/ph(1,1), ... % mortality fraction due to VH
      -(ph(2,6)+ph(2,10))/ph(2,2), ... % mortality fraction due to VC
      -(ph(3,7)+ph(3,10))/ph(3,3)];    % mortality fraction due to VE
    % construct vector with all values to be compared to targets
    cmp = [sol(:)' sum(sol(5:7))/sum(sol(1:3)) aux];
    % compute Euclidean distance to targets
    dev(cc) = sum((log(cmp) - log(tar)).^2);
end

% select the 1000 realisations that are closest to targets
nReal = 1000;
[~,idx] = mink(dev,nReal);

% store the selected realisations
res(nReal) = struct('pars',[],'sol',[]);
for k = 1:nReal
    i = idx(k);
    res(k).pars = ps(i);
    res(k).sol  = model_WEI_eqs(ps(i));
end

%% Construct flux matrix
% for the 1000 selected realisations, extract the flux matrix
% in the format used in Weitz et al. (2015) and transform
% into the format used in Pruvôt et al. (2026)

% format used in Weitz et al. (2015): matrix of dimensions 10x10
% -> variables in the order H,C,E,Z,VH,VC,VE,xo,xi,DV
% format used in Pruvôt et al. (2026): matrix of dimensions 5x5
% -> variables in the order P,B,G,D,DV

F_real = struct('mat', cell(1,1000));
for c1 = 1:1e3
    pn = ps(c1);                        % current parameter set
    xeq = model_WEI_eqs(pn);            % corresponding equilibrium
    [~,F_ori] = model_WEI_rhs(xeq,pn);  % flux matrix from Weitz et al.

    F_trf = zeros(5);

    % Phytoplankton
    F_trf(1,1) = F_ori(2,2)+F_ori(3,3);                    % primary production
    F_trf(1,3) = sum(F_ori(2:3,4));                        % grazing
    F_trf(1,4) = sum(F_ori(2:3,8));                        % DOC release
    F_trf(1,5) = sum(F_ori(2:3,10))+F_ori(2,6)+F_ori(3,7); % viral lysis

    % Heterotrophic bacteria
    F_trf(2,2) = F_ori(1,1);                               % bacterial growth
    F_trf(2,3) = F_ori(1,4);                               % grazing
    F_trf(2,4) = F_ori(1,8);                               % DOC release
    F_trf(2,5) = F_ori(1,5)+F_ori(1,10);                   % viral lysis

    % Grazers
    F_trf(3,3) = F_ori(4,4);                               % growth grazers
    F_trf(3,4) = F_ori(4,8);                               % DOC release

    % DOC
    F_trf(4,2) = F_ori(8,1);                               % DOC consumption
    F_trf(4,4) = F_ori(8,8);                               % incoming DOC

    % DOC_V
    F_trf(5,4) = sum(F_trf(1:2,5));                        % virus-driven DOC
    F_trf(5,5) = -sum(F_trf(1:2,5));                       % to total DOC
    
    F_trf = F_trf./F_trf(1,1)*100;
    F_real(c1).mat = F_trf;
end

%% Compute mean of all flux matrices

F_all = cat(3, F_real.mat);     % matrix of dimensions 5×5×1000
F_mean = mean(F_all, 3);   % by averaging over dimension 3, we get matrix of dimensions 5×5

end
