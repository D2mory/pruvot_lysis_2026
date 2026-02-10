% data_genetation_virus_removal.m
%
% This code generates the data used in section 5 on virus removal
%  . 10000 plausible realisations of a dynamical model are generated
%  . Viruses are removed and the model is computed
%  . Simulations that reach steady states without viruses are kept for analysis (~8700)
%  . Flux matrices are extracted for the selected realisations, with and  without viruses

rng(1); %fixed seed

%% Fix parameter ranges
pr=struct('qP',[3e-7 3e-6], ...
          'qB',[2e-9 2e-8], ...
          'qG',[2e-4 2e-3], ...
          'qVP',[2e-12 6e-11], ...
          'qVB',[2e-12 6e-11], ...
          'muP',[.2 2], ...
          'KP',[1e7 1e9], ...
          'm1P',[.005 .1], ...
          'r1P',[.001 .1], ...
          'psP',[1e-6 2e-5], ...
          'muB',[.5 2], ...
          'KB',[30 300], ...
          'eB',[.05 .2], ...
          'm1B',[.005 .1], ...
          'r1B',[.001 .1], ...
          'psB',[1e-6 3e-5], ...
          'm2G',[1e-7 1e-5],...
          'r1G',[1e-3 1e-2], ...
          'pi',[.3 .5], ...
          'pe',[.2 .4], ...
          'phP',[1e-9 1e-8], ...
          'beP',[50 500], ...
          'mVP',[.05 5], ...
          'phB',[1e-11 1e-9], ...
          'beB',[10 100], ...
          'mVB',[.05 5]);

% Fix variable ranges for acceptable densities
vr=[2e5 2e7;   % P
    2e7 2e9;   % B
    2e6 2e8;   % VP
    2e8 2e10;  % VB
    2e3 2e5;   % G
    30 100];   % D


%% Generate 10000 parameter sets to find plausible realisations

nSim = 1e4;
pars(nSim) = pr; % structure to store accepted parameter sets
c1 = 1;
c2 = 1;
while c2 <= nSim
    pn = vr_smp(pr); %uniform sampling of parameters
    xeq = vr_eqs(pn); %variables steady state
    if vr_chk(xeq,vr) %check if variables' equilibria are within the desired ranges
        if mod(c2,10) == 0
            fprintf('Found parameter set %d at iteration %d\n', c2, c1)
        end
        pars(c2) = pn; % All selected parameter sets
        c2 = c2+1;
    end
    c1 = c1+1;
end

%% Virus removal

% Structure "sols" to store all solutions
solref=struct('with',[],'withC',[],'wout',[],'woutC',[]);
% 'with' : solution with viruses, abundances
% 'withC' : solution with viruses, carbon content
% 'wout' : solution without viruses, abundances
% 'woutC' : solution without viruses, carbon content
sols(nSim)=solref;

nrms=zeros(1,nSim);
for c1=1:nSim
    if mod(c1,10) == 0
        fprintf('Starting parameter set %d\n', c1)
    end

    %%% With viruses %%%
    pn = pars(c1); %parameters for abundances
    pc = vr_trf(pn,3,pn); %parameters from abundance to carbon
    xeq = vr_eqs(pn); %variables in abundance
    sols(c1).with=xeq;
    xeqC = vr_trf(xeq,1,pn); %variables from abundance to carbon
    sols(c1).withC=xeqC;

    %%% Without viruses %%%
    % The model is simulated, the initial conditions are the steady states
    % of the model obtained in the presence of viruses from which we remove
    % viruses.
    x0C = xeqC;
    x0C(3:4) = 0; %virus removal in the initial conditions
    opts = odeset('abstol',1e-8,'reltol',1e-8,'nonnegative',1:6); 
    [TC,XC] = ode45(@(t,x) vr_rhsC(x,pc),[0 1e4],x0C,opts); %simulation
    xenC = XC(end,:); 
    sols(c1).woutC=xenC; %carbon "equilibrium" without viruses
    sols(c1).wout=vr_trf(xenC,2,pn); %from carbon to abundances without viruses
    nrms(c1)=norm(vr_rhsC(xenC,pc)); %how far is the system from steady state at xenC
end

%% Construct flux matrix and solve fluxes
% For the simulations that reached a steady state, the flux matrix is computed

% indices of cases at steady state
idx = find(nrms < 1e-6);
nsel = numel(idx);

% preallocate only selected cases
flxref = struct('with', [], 'wout', []);
flxs = repmat(flxref, nsel, 1);

xeqref = struct('with', [], 'wout', []);
xeqsC = repmat(xeqref, nsel, 1);

for k = 1:nsel
    c1 = idx(k);

    pn = pars(c1);
    pc = vr_trf(pn, 3, pn); % parameters in carbon

    xeqC = sols(c1).withC; % steady state with viruses
    [~, ph_with] = vr_rhsC(xeqC, pc);
    ph_with = vr_flx_trf(ph_with(:,1:7));

    xenC = sols(c1).woutC; % steady state without viruses
    [~, ph_wout] = vr_rhsC(xenC, pc);
    ph_wout = vr_flx_trf(ph_wout(:,1:7));

    flxs(k).with = ph_with;
    flxs(k).wout = ph_wout;

    xeqsC(k).with = xeqC;
    xeqsC(k).wout = xenC';
end

psel = pars(idx);

%% Saving data

save('data_virus_removal.mat', 'flxs', 'psel', 'xeqsC');